#!/usr/bin/env python3
"""Read-only evidence review of frozen baseline annotation conflicts.

Never changes rules or labels. Replays scoring only to retain partitions and
requires exact metric parity. Missing partition matches are unresolved, not FP.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import subprocess
import zlib

import pysam

from aggregate_benchmark import parse_metrics_json
from audit_refined_rescue_conflicts import file_stat, paired_evidence
from apply_refined_rescue import digest, info_dict
from vcf_utils.refined_rescue_policy import biological_veto
from validate_refined_native_integration import alleles


def af_comparison(info, matches):
    raw = next((v for k, v in info.items() if k.casefold() == 'gnomad_af'), None)
    if raw in (None, '', '.'):
        return {'source_af': raw, 'af_status': 'unavailable'}
    try:
        value = float(raw)
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError('invalid AF')
    except (ValueError, TypeError):
        return {'source_af': raw, 'af_status': 'invalid_or_non_scalar'}
    confirmed = any(m['af'] is not None and math.isclose(value, m['af'], rel_tol=1e-5, abs_tol=1e-8)
                    for m in matches)
    return {'source_af': value, 'af_status': 'exact_allele_confirmed' if confirmed else 'not_confirmed'}


def partition_label(key, partitions):
    found = [label for label, keys in partitions.items() if key in keys]
    return found or ['unscored_or_representation_unresolved']


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--validation', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    ap.add_argument('--gnomad', type=Path, default=Path('/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/gnomAD/exomes'))
    args = ap.parse_args()
    source = json.loads(args.validation.read_text())
    if source['status'] != 'complete_experimental_not_training_ready':
        raise ValueError('Require a completed frozen evaluation')
    root = args.validation.parent.resolve()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    report = {'status': 'running', 'scope': __doc__, 'sources': {}, 'database_stats': {},
              'sqlite_stats': {}, 'rounds': {}, 'benchmarks': {}, 'commands': [],
              'code': {str(Path(__file__).resolve()): digest(Path(__file__))}}

    def save():
        (out / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')

    def track(path):
        path = Path(path).resolve()
        report['sources'][str(path)] = digest(path)

    try:
        track(args.validation)
        baseline = Path(source['outputs']['refined_consensus']['vcf'])
        track(baseline)
        if digest(baseline) != source['outputs']['refined_consensus']['sha256']:
            raise ValueError('Baseline changed since evaluation')
        positive = alleles(baseline, True)
        for round_name in ('first', 'realignment'):
            dbpath = root / round_name / 'union.sqlite'
            report['sqlite_stats'][str(dbpath)] = file_stat(dbpath)
            adapter_report = root / round_name / 'report.json'
            track(adapter_report)
            expected = json.loads(adapter_report.read_text())['counts'].get('baseline_annotation_conflicts', 0)
            entries = []
            with sqlite3.connect(dbpath.as_uri() + '?mode=ro', uri=True) as db:
                for key in sorted(positive):
                    row = db.execute('SELECT rescue FROM variants WHERE chrom=? AND pos=? AND ref=? AND alt=?', key).fetchone()
                    if not row or row[0] is None:
                        continue
                    raw = zlib.decompress(row[0]).decode()
                    fields = raw.split('\t')
                    info = info_dict(fields[7])
                    reason = biological_veto(info)
                    if reason:
                        entries.append({'allele': key, 'reason': reason, 'source_label': fields[6],
                                        'source_info': info, 'source_record_sha256': hashlib.sha256(raw.encode()).hexdigest()})
            if len(entries) != expected:
                raise ValueError(f'{round_name}: conflict count {len(entries)} != recorded {expected}')
            report['rounds'][round_name] = {'expected_conflicts': expected, 'conflicts': entries}
        save()

        # Reuse the frozen full query, never select only the flagged sites.
        query = root / 'refined_consensus.query.vcf.gz'
        track(query)
        for region in ('ukb', 'medexome'):
            prefix = root / f'refined_consensus.{region}'
            original = next(c for c in source['commands'] if c[-1] == str(prefix))
            cmd = list(original)
            cmd[-1] = str(out / region)
            scratch = out / f'{region}_scratch'
            cmd += ['--scratch-prefix', str(scratch)]
            report['commands'].append(cmd)
            save()
            with (out / f'{region}.log').open('w') as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
            metrics = parse_metrics_json(out / f'{region}.metrics.json')
            if metrics != source['metrics'][f'refined_consensus/{region}']:
                raise ValueError(f'{region}: benchmark replay parity failure')
            partitions = {label: alleles(scratch / path) for label, path in
                          {'TP': 'tpfn/0003.vcf.gz', 'FP': 'fp.vcf.gz',
                           'ambiguous': 'ambi.vcf.gz', 'unscored': 'unk.vcf.gz'}.items()}
            report['benchmarks'][region] = {'metrics': metrics, 'frozen_parity': True}
            for rd in report['rounds'].values():
                for row in rd['conflicts']:
                    row.setdefault('partitions', {})[region] = partition_label(tuple(row['allele']), partitions)
            save()

        handles = {}
        try:
            for rd in report['rounds'].values():
                for row in rd['conflicts']:
                    chrom, pos, ref, alt = row['allele']
                    if chrom not in handles:
                        path = args.gnomad / f'gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz'
                        report['database_stats'][str(path)] = file_stat(path)
                        handles[chrom] = pysam.VariantFile(str(path))
                    matches = []
                    for record in handles[chrom].fetch(chrom, pos-1, pos):
                        if record.pos != pos or record.ref != ref:
                            continue
                        for i, allele in enumerate(record.alts or []):
                            if allele == alt:
                                af = record.info.get('AF')
                                value = af[i] if isinstance(af, tuple) else af
                                matches.append({'af': value, 'filters': list(record.filter),
                                                'record_sha256': hashlib.sha256(str(record).encode()).hexdigest()})
                    row['database_matches'] = matches
                    row.update(af_comparison(row['source_info'], matches))
        finally:
            for reader in handles.values():
                reader.close()
        for caller, path in source['manifest']['dna'].items():
            track(path)
            with pysam.VariantFile(path) as reader:
                for rd in report['rounds'].values():
                    for row in rd['conflicts']:
                        row.setdefault('paired_evidence', {})[caller] = paired_evidence(reader, caller, row['allele'])
        for rd in report['rounds'].values():
            rd['reasons'] = dict(Counter(r['reason'] for r in rd['conflicts']))
            rd['af_status'] = dict(Counter(r['af_status'] for r in rd['conflicts']))
            rd['partition_counts'] = {region: dict(Counter(label for r in rd['conflicts'] for label in r['partitions'][region]))
                                      for region in ('ukb', 'medexome')}
        report['sources_unchanged'] = all(digest(Path(p)) == h for p, h in report['sources'].items())
        report['large_files_stats_unchanged'] = all(file_stat(Path(p)) == s for group in ('database_stats', 'sqlite_stats')
                                                   for p, s in report[group].items())
        report['status'] = ('complete_read_only_not_training_approval' if report['sources_unchanged'] and
                            report['large_files_stats_unchanged'] else 'integrity_failure')
        save()
        print(json.dumps({k: {f: v[f] for f in ('expected_conflicts', 'reasons', 'af_status', 'partition_counts')}
                          for k, v in report['rounds'].items()}, indent=2), flush=True)
    except Exception as exc:
        report.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        save()
        raise


if __name__ == '__main__':
    main()
