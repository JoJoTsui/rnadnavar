#!/usr/bin/env python3
"""Audit immutable v2 exports and prepare quarantined weak-label review subsets.

No mapping/calling, no new absence-of-ALT thresholds, no training approval.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pysam

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'bin'))
from vcf_utils.three_class_policy import native_negative
from build_refined_rerun_artifacts import sha

EXCLUDED = {'PRJNA298330_4032', 'PRJNA298376_4081', 'PRJNA298376_4255'}
THREE = {'PRJNA298376_4007', 'PRJNA298376_4060', 'PRJNA298376_4072'}
POLICY = 'separated_three_class_v2'


def selection(table):
    allowed = pc.invert(pc.is_in(table['sample_id'], value_set=pa.array(sorted(EXCLUDED))))
    # Quarantine explicit conflicts, not weak labels merely awaiting approval.
    reason = pc.fill_null(table['review_reason'], '')
    conflict = pc.and_(pc.equal(table['FILTER'], 'Somatic'),
                       pc.not_equal(reason, 'somatic_candidate_requires_label_qc'))
    return pc.and_(allowed, pc.invert(conflict))


def counts(table):
    return Counter({r['sample_id']+'|'+r['FILTER']: r['POS_count'] for r in
        table.group_by(['sample_id', 'FILTER'], use_threads=False).aggregate([('POS','count')]).to_pylist()})


def run(root, out):
    summary_path, manifest_path = root/'exports/summary.json', root/'exports/manifest.tsv'
    verification_path = root/'review_20260921/parquet_verification.json'
    summary, verification = (json.loads(p.read_text()) for p in (summary_path, verification_path))
    source = Path(summary['parquet'])
    bindings = {str(p): sha(p) for p in (summary_path, manifest_path, verification_path, source,
                                       Path(__file__), ROOT/'bin/vcf_utils/three_class_policy.py')}
    if (verification['status'] != 'all_parquet_rows_verified'
            or verification['parquet_sha256'] != bindings[str(source)]
            or verification['summary_sha256'] != bindings[str(summary_path)]
            or verification['manifest_sha256'] != bindings[str(manifest_path)]):
        raise ValueError('Invalid source export binding')
    with manifest_path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    manifest = {r['sample_id']: r for r in rows}
    if len(rows) != 66 or len(manifest) != 66 or not EXCLUDED|THREE <= set(manifest):
        raise ValueError('Unexpected cohort identity')
    for sid, row in manifest.items():
        state = json.loads((root/sid/'state.json').read_text())
        if (state['status'] != 'candidate_complete_not_training_approved'
                or state['identity']['policy'] != POLICY
                or row['candidate_policy'] != POLICY
                or row['rescue_vcf_path'] != state['final_artifacts']['rescue']
                or row['new_rescue_sha256'] != state['outputs'][row['rescue_vcf_path']]):
            raise ValueError('Wrong sample/version binding: '+sid)
    out.mkdir(parents=True, exist_ok=False)
    report = dict(status='running', training_approved=False, sources=bindings,
                  excluded_samples=sorted(EXCLUDED), three_sample_stage=sorted(THREE),
                  single_sample_stage_skipped=True, policy=POLICY,
                  limitations=['Weak workflow supervision, not independent truth',
                               'Negative native-caller spot checks are bounded, not class accuracy',
                               'Full original BAM/sample biological identity not certified',
                               'No 299-read or zero-ALT training gate imposed',
                               'No global training split or approval assigned'])
    all_counts, kept_counts, three_counts, dropped = Counter(), Counter(), Counter(), Counter()
    probes = {}
    try:
        reader = pq.ParquetFile(source)
        with pq.ParquetWriter(out/'cohort63.review.parquet', reader.schema_arrow, compression='zstd') as writer, \
             pq.ParquetWriter(out/'three_samples.review.parquet', reader.schema_arrow, compression='zstd') as small:
            for batch in reader.iter_batches(batch_size=32768, use_threads=False):
                table = pa.Table.from_batches([batch])
                bad = pc.or_(pc.not_equal(table['THREE_CLASS_POLICY'], POLICY), table['training_eligible'])
                if pc.any(pc.fill_null(bad, True)).as_py():
                    raise ValueError('Wrong policy or unexpected source approval')
                negative = pc.is_in(table['FILTER'], value_set=pa.array(['Germline','Reference']))
                valid_native = pc.and_(pc.equal(table['FILTER'], table['THREE_CLASS_NATIVE_FILTER']),
                    pc.match_substring(table['THREE_CLASS_NATIVE_RATIONALE'], 'three_class_policy:native_three_class_v1'))
                if pc.any(pc.and_(negative, pc.invert(pc.fill_null(valid_native, False)))).as_py():
                    raise ValueError('Negative class without native nomination provenance')
                if pc.any(pc.not_equal(table['TRAINING_ELIGIBLE'], 'NO')).as_py():
                    raise ValueError('Unexpected source training flag')
                all_counts.update(counts(table))
                mask = selection(table)
                rejected = table.filter(pc.invert(mask))
                dropped.update(counts(rejected))
                retained = table.filter(mask)
                kept_counts.update(counts(retained))
                writer.write_table(retained)
                training_stage = retained.filter(pc.is_in(retained['sample_id'], value_set=pa.array(sorted(THREE))))
                three_counts.update(counts(training_stage))
                small.write_table(training_stage)
                check = table.filter(pc.and_(negative, pc.is_in(table['sample_id'], value_set=pa.array(sorted(EXCLUDED|THREE)))))
                # Reproducible prefix per sample/class: a plumbing check, not prevalence estimation.
                for item in check.select(['sample_id','FILTER','CHROM','POS','REF','ALT']).to_pylist():
                    key = (item['sample_id'], item['FILTER'])
                    if len(probes.setdefault(key, [])) < 4:
                        probes[key].append(item)
        expected = Counter({s['sample_id']+'|'+label: n for s in summary['samples_detail']
                            for label,n in s['exported_counts'].items()})
        if all_counts != expected or set(k.split('|')[0] for k in kept_counts) != set(manifest)-EXCLUDED:
            raise ValueError('Cohort class coverage mismatch')
        report['native_evidence_checks'] = []
        for sid in sorted(EXCLUDED|THREE):
            path = Path(manifest[sid]['caller_dna_deepsomatic'])
            before = (path.stat().st_size, path.stat().st_mtime_ns)
            with pysam.VariantFile(str(path)) as vcf:
                names = list(vcf.header.samples)
                if len(names) != 1 or names[0] not in {sid+'DT',manifest[sid]['vcf_prefix']+'DT'}:
                    raise ValueError('Unexpected native caller sample header: '+sid+str(names))
                for label in ('Germline','Reference'):
                    for item in probes.get((sid,label), []):
                        matches = [r for r in vcf.fetch(item['CHROM'],item['POS']-1,item['POS'])
                                   if r.pos == item['POS'] and r.ref == item['REF'] and r.alts == (item['ALT'],)]
                        if len(matches) != 1:
                            raise ValueError('Native source allele not uniquely found: '+str(item))
                        rec = matches[0]
                        raw = rec.samples[names[0]]
                        ev = {'tumor_'+k:raw.get(k) for k in ('GQ','DP','AD','PL')}
                        observed, reason = native_negative(dict(REF=rec.ref, ALT=rec.alts[0], callers=['deepsomatic'],
                            filters_original=[';'.join(rec.filter)], native_evidence={'deepsomatic':ev}))
                        if observed != label:
                            raise ValueError('Native source evidence does not reproduce nomination: '+str(item))
                        report['native_evidence_checks'].append(dict(**item, evidence=ev, reason=reason,
                            caller_vcf=str(path), caller_sample=names[0]))
            if before != (path.stat().st_size,path.stat().st_mtime_ns):
                raise ValueError('Native caller changed')
            print(sid+': current native source checks passed; '+('quarantined' if sid in EXCLUDED else 'three-sample stage'), flush=True)
        for name, ids, parquet in [('cohort63.review.tsv',set(manifest)-EXCLUDED,'cohort63.review.parquet'),
                                  ('three_samples.review.tsv',THREE,'three_samples.review.parquet')]:
            with (out/name).open('x') as handle:
                fields = list(rows[0])+['handoff_status','source_variant_parquet']
                csvwriter = csv.DictWriter(handle,fieldnames=fields,delimiter='\t')
                csvwriter.writeheader()
                for sid in sorted(ids):
                    row = dict(manifest[sid], source_variant_parquet=str(source),
                               variant_parquet_path=str(out/parquet), handoff_status='review_ready_not_training_approved')
                    row['training_label_vcf'] = ''
                    csvwriter.writerow(row)
        if any(sha(Path(p)) != h for p,h in bindings.items()):
            raise ValueError('Source export/code changed during audit')
        report.update(status='review_ready_not_training_approved', source_counts=dict(all_counts),
                      selected_counts=dict(kept_counts), three_sample_counts=dict(three_counts),
                      quarantined_counts=dict(dropped), outputs={str(p):sha(p) for p in out.iterdir() if p.is_file()})
    except Exception as exc:
        report.update(status='failed',error=str(exc))
        raise
    finally:
        with (out/'audit.json').open('x') as handle:
            json.dump(report,handle,indent=2)
            handle.write('\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    args = parser.parse_args()
    print(run(args.source_root,args.outdir)['status'])
