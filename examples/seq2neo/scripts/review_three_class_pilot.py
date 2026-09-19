#!/usr/bin/env python3
"""Read-only pilot review, with scratch/reports in a fresh local directory.

Independently replay the allele/class union (not the production adapter), check
all recorded file hashes, and show observed stage timings. This does not grant
training approval or execute a cohort. Large allele sets are disk-backed.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import sqlite3
import time

import pysam

import run_refined_cohort as cohort


def key(record):
    return (record.contig, record.pos, record.ref, ','.join(record.alts or ()))


def label(record):
    values = list(record.filter)
    if len(values) != 1:
        raise ValueError('Require exactly one biological FILTER')
    return values[0]


def review_stage(baseline, native, final, scratch):
    """Exact full-allele union, no truth lookup and no in-memory record sets."""
    counts = Counter()
    with sqlite3.connect(scratch) as db:
        db.execute('PRAGMA cache_size=-16384')
        db.execute('CREATE TABLE alleles(c TEXT,p INTEGER,r TEXT,a TEXT,b TEXT,n TEXT,'
                   'seen INTEGER DEFAULT 0,PRIMARY KEY(c,p,r,a)) WITHOUT ROWID')
        for path, column in ((baseline, 'b'), (native, 'n')):
            with pysam.VariantFile(str(path)) as reader:
                for record in reader:
                    value = label(record)
                    if column == 'n' and value in {'Germline', 'Reference'}:
                        trace = str(record.info.get('CLASSIFICATION_RATIONALE', '')).split('|')
                        if 'three_class_policy:native_three_class_v1' not in trace:
                            raise ValueError('Unproven native negative nomination')
                    cursor = db.execute(
                        f'INSERT INTO alleles(c,p,r,a,{column}) VALUES (?,?,?,?,?) '
                        f'ON CONFLICT(c,p,r,a) DO UPDATE SET {column}=excluded.{column} '
                        f'WHERE alleles.{column} IS NULL', (*key(record), value))
                    if cursor.rowcount != 1:
                        raise ValueError('Duplicate source allele')
            db.commit()
        with pysam.VariantFile(str(final)) as reader:
            for record in reader:
                row = db.execute('SELECT b,n,seen FROM alleles WHERE c=? AND p=? AND r=? AND a=?',
                                 key(record)).fetchone()
                if row is None or row[2]:
                    raise ValueError('Unexpected or duplicate final allele')
                b, n, _ = row
                expected = 'Somatic' if b == 'Somatic' else n if n in {'Germline', 'Reference'} else 'NoConsensus'
                if label(record) != expected:
                    raise ValueError(f'Class/membership mismatch: {key(record)}')
                if (record.info.get('TRAINING_ELIGIBLE') != 'NO'
                        or record.info.get('UNIFIED_FILTER') != expected
                        or record.info.get('THREE_CLASS_POLICY') != cohort.THREE_CLASS_POLICY
                        or record.info.get('THREE_CLASS_BASELINE_FILTER') != (b or 'MISSING')
                        or record.info.get('THREE_CLASS_NATIVE_FILTER') != (n or 'MISSING')):
                    raise ValueError('Eligibility/class provenance mismatch')
                counts[expected] += 1
                db.execute('UPDATE alleles SET seen=1 WHERE c=? AND p=? AND r=? AND a=?', key(record))
        if db.execute('SELECT count(*) FROM alleles WHERE seen=0').fetchone()[0]:
            raise ValueError('Missing final alleles')
    return dict(status='exact_all_class_union_pass', counts=dict(counts), records=sum(counts.values()))


def timings(dest, state):
    # Logs without contents retain creation/start mtime; successful command
    # outputs/reports bound the preceding step's finish. These are estimates,
    # not CPU profiling or guaranteed per-command duration telemetry.
    boundaries = [
        ('normalization', (dest/'index_strelka.log').stat().st_mtime),
        ('baseline_consensus', (dest/'index_consensus.log').stat().st_mtime),
        ('baseline_rescue', (dest/'refined.rescue.vcf.gz').stat().st_mtime),
        ('native_nomination', (dest/'index_native.log').stat().st_mtime),
        ('three_class_consensus', (dest/'three_class.consensus.report.json').stat().st_mtime),
        ('three_class_rescue', (dest/'three_class.rescue.report.json').stat().st_mtime),
        ('audit_and_integrity', state['completed']),
    ]
    prior = state['started']
    result = {}
    for name, stamp in boundaries:
        if stamp < prior or stamp > state['completed'] + 1:
            return dict(status='unavailable_modified_timestamps')
        result[name] = round((stamp - prior) / 60, 3)
        prior = stamp
    return dict(status='estimated_from_local_artifact_mtimes', minutes=result,
                total_minutes=(state['completed'] - state['started']) / 60)


def review(cfg, outdir):
    root = Path(cfg['output_root']).resolve(strict=True)
    outdir = outdir.resolve()
    if cohort.below(outdir, root) or cohort.below(root, outdir):
        raise ValueError('Review destination must be disjoint from cohort output')
    frozen = json.loads((root/'run_identity.json').read_text())
    gate = cohort.validation_gate(cfg)
    if (frozen['policy'] != cohort.THREE_CLASS_POLICY or frozen['code'] != cohort.code_hashes()
            or frozen['validation_gate'] != gate['sha256']
            or frozen['manifest'] != cohort.digest(cfg['manifest'])
            or frozen['reference'] != cohort.digest(cfg['fasta'])
            or frozen['reference_fai'] != cohort.digest(cfg['fasta'] + '.fai')):
        raise ValueError('Pilot identity no longer matches code/inputs/validation')
    outdir.mkdir(parents=True, exist_ok=False)
    result = dict(status='running', training_approved=False, policy=cohort.THREE_CLASS_POLICY,
                  root=str(root), identity=frozen, samples={}, started=time.time())
    try:
        for sid in cfg['pilot_samples']:
            state_path = root/sid/'state.json'
            state = json.loads(state_path.read_text())
            source_hashes = {p:cohort.digest(p) for p in state['source_hashes']}
            if not cohort.completed(root/sid, frozen, source_hashes):
                raise ValueError('Incomplete/stale/corrupt pilot: ' + sid)
            dest = Path(state['output'])
            completion = dest/'completion.json'
            if json.loads(completion.read_text()) != state:
                raise ValueError('State/completion disagreement')
            sample = dict(completion=str(completion), completion_sha256=cohort.digest(completion),
                          timing=timings(dest, state), stages={})
            print(sid + ': hashes verified; reviewing full allele unions', flush=True)
            for stage, base in (('consensus', 'refined.vcf.gz'), ('rescue', 'refined.rescue.vcf.gz')):
                final = Path(state['final_artifacts'][stage])
                outcome = review_stage(dest/base, dest/'native.vcf.gz', final, outdir/f'{sid}.{stage}.sqlite')
                report_path = dest/f'three_class.{stage}.report.json'
                report = json.loads(report_path.read_text())
                audit_path = dest/f'{stage}.audit.json'
                audit = json.loads(audit_path.read_text())
                if (report['somatic_membership_mismatches'] != 0 or not report['sources_unchanged']
                        or audit['issues'] or not audit['sources_unchanged']
                        or audit['sha256'] != state['outputs'][str(final)]
                        or any(audit['counts'].get(k) != v for k,v in outcome['counts'].items())):
                    raise ValueError('Pilot report/audit mismatch')
                outcome.update(final=str(final), sha256=state['outputs'][str(final)],
                               report_sha256=cohort.digest(report_path), audit_sha256=cohort.digest(audit_path))
                sample['stages'][stage] = outcome
            # Detect changes during the independent replay, not just before it.
            if (any(cohort.digest(p) != h for p,h in source_hashes.items())
                    or not cohort.completed(root/sid, frozen, source_hashes)
                    or json.loads(state_path.read_text()) != state):
                raise ValueError('Pilot inputs/outputs/state changed during review')
            result['samples'][sid] = sample
            print(sid + ': exact three-class union PASS', flush=True)
        if cohort.code_hashes() != frozen['code'] or cohort.validation_gate(cfg)['sha256'] != frozen['validation_gate']:
            raise ValueError('Policy code/validation changed during review')
        result.update(status='candidate_execution_review_pass', completed=time.time())
    except Exception as exc:
        result.update(status='failed', error=str(exc))
        raise
    finally:
        (outdir/'review.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    review(json.loads(args.config.read_text()), args.outdir)


if __name__ == '__main__':
    main()
