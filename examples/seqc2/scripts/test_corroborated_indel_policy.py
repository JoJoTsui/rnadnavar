#!/usr/bin/env python3
"""Bounded SEQC2 development assay; no production changes or HG008 inputs.

Retain current consensus and add indels with DeepSomatic PASS Q>=20, tumor
AD>=3, Mutect2 tumor AD>=2, normal DP>=10/ADalt=0, GERMQ>=20 and TLOD>0.
Only PASS or contamination/weak_evidence Mutect2 filters are eligible.
Callers observe the same reads; their counts are never summed.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess

from aggregate_benchmark import parse_metrics_json
from audit_current_native_policy import digest
from attribute_native_indel_errors import rows
from replay_historical_native_gate import write_query


def qualifies(ds, m2):
    try:
        quality = float(ds['qual'])
        if ds['filter'] not in {'PASS', '.'} or not math.isfinite(quality) or quality < 20:
            return False
        filters = set(m2['filter'].split(';')) - {'PASS', '.'}
        if not filters <= {'contamination', 'weak_evidence'}:
            return False
        dt = [v for k, v in ds['samples'].items() if k.endswith('_T_1')]
        mt = [v for k, v in m2['samples'].items() if k.endswith('_T_1')]
        mn = [v for k, v in m2['samples'].items() if k.endswith('_N_1')]
        if not all(len(v) == 1 for v in (dt, mt, mn)):
            return False
        info = dict(v.split('=', 1) for v in m2['info'].split(';') if '=' in v)
        germq, tlod = float(info['GERMQ']), float(info['TLOD'])
        return (int(dt[0]['AD'].split(',')[1]) >= 3 and
                int(mt[0]['AD'].split(',')[1]) >= 2 and
                int(mn[0]['AD'].split(',')[1]) == 0 and int(mn[0]['DP']) >= 10 and
                math.isfinite(germq) and germq >= 20 and math.isfinite(tlod) and tlod > 0)
    except (KeyError, ValueError, IndexError):
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--audit', required=True, type=Path)
    ap.add_argument('--evidence', required=True, type=Path)
    ap.add_argument('--dataset', required=True, choices=['wes_ll', 'wgs_il'])
    ap.add_argument('--outdir', required=True, type=Path)
    ap.add_argument('--strict-snv-evidence', action='store_true',
                    help='Reject non-DeepSomatic-PASS SNV additions with other Mutect2 rejection filters')
    args = ap.parse_args()
    audit = json.loads(args.audit.read_text())
    args.outdir = args.outdir.resolve()
    args.outdir.mkdir(parents=True, exist_ok=False)
    results = {'status': 'exploratory_not_promoted', 'rule': __doc__, 'cells': {},
               'strict_snv_evidence': args.strict_snv_evidence,
               'source_audit': str(args.audit.resolve()), 'source_audit_sha256': digest(args.audit),
               'script_sha256': digest(Path(__file__)), 'source_sha256': {}}
    for region in ('ukb', 'medexome'):
        source = args.evidence / region
        for path in (source / 'current_native_scratch/normalized_query.vcf.gz',
                     source / 'deepsomatic.normalized.vcf.gz', source / 'mutect2.normalized.vcf.gz'):
            results['source_sha256'][str(path.resolve())] = digest(path)
        # Selection reads caller evidence and current query only, never truth.
        current = dict(rows(source / 'current_native_scratch/normalized_query.vcf.gz'))
        ds = dict(rows(source / 'deepsomatic.normalized.vcf.gz'))
        m2 = dict(rows(source / 'mutect2.normalized.vcf.gz'))
        removed_snvs = set()
        if args.strict_snv_evidence:
            for key in current:
                if len(key[2]) != 1 or len(key[3]) != 1:
                    continue
                if ds.get(key, {}).get('filter') in {'PASS', '.'}:
                    continue
                record = m2.get(key)
                if record is None or not (set(record['filter'].split(';')) - {'PASS', '.'}) <= {'contamination', 'weak_evidence'}:
                    removed_snvs.add(key)
        additions = {k for k in ds.keys() & m2.keys() if k not in current
                     and ',' not in k[3] and len(k[2]) != len(k[3])
                     and qualifies(ds[k], m2[k])}
        dest = args.outdir / region
        dest.mkdir()
        query = write_query(dest / 'query.vcf', (set(current) - removed_snvs) | additions)
        cmd = next(list(c) for c in audit['commands']
                   if c[-1].endswith(f'/{args.dataset}/{region}/current_native'))
        cmd[6] = query
        cmd[-1] = str(dest / 'benchmark')
        with (dest / 'benchmark.log').open('w') as handle:
            subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=True)
        measured = parse_metrics_json(dest / 'benchmark.metrics.json')
        results['cells'][region] = {'additions': sorted(additions), 'removed_snvs': sorted(removed_snvs),
                                    'metrics': measured, 'command': cmd}
        (args.outdir / 'evaluation.json').write_text(json.dumps(results, indent=2) + '\n')
        print(args.dataset, region, len(additions), json.dumps(measured), flush=True)


if __name__ == '__main__':
    main()
