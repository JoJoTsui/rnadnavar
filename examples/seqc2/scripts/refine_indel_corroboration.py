#!/usr/bin/env python3
"""SEQC2-informed refinement; requires later frozen HG008 evaluation.

Test two limited additions to the fixed candidate: reciprocal Mutect2 PASS
support only for ECNT=1 events; or DeepSomatic PASS QUAL>=10 with >=3 tumor
alternate reads, >=2 Mutect2 alternate reads spanning both strands, and clean
Mutect2 normal evidence (DP>=10, alt=0, GERMQ>=20). Mutect2 may only carry
contamination/weak_evidence filters and must have positive TLOD.
"""
import argparse
import json
from pathlib import Path
import subprocess

from aggregate_benchmark import parse_metrics_json
from audit_current_native_policy import digest
from attribute_native_indel_errors import rows
from explore_indel_corroboration import gates, number, sample
from replay_historical_native_gate import write_query


def refined_gates(ds, m2):
    result = set()
    try:
        info = dict(v.split('=', 1) for v in m2['info'].split(';') if '=' in v)
        if 'mutect_pass' in gates(ds, m2, {}) and int(info['ECNT']) == 1:
            result.add('single_event')
        if ds['filter'] not in {'PASS', '.'} or number(ds['qual']) < 10:
            return result
        if (set(m2['filter'].split(';')) - {'PASS', '.'}) - {'contamination', 'weak_evidence'}:
            return result
        dt, mt, mn = sample(ds, '_T_1'), sample(m2, '_T_1'), sample(m2, '_N_1')
        strand = [int(x) for x in mt['SB'].split(',')]
        if (int(dt['AD'].split(',')[1]) >= 3 and int(mt['AD'].split(',')[1]) >= 2
                and int(mn['DP']) >= 10 and int(mn['AD'].split(',')[1]) == 0
                and number(info['GERMQ']) >= 20 and number(info['TLOD']) > 0
                and len(strand) == 4 and min(strand[2:]) >= 1):
            result.add('moderate_strand')
    except (KeyError, ValueError, IndexError, TypeError):
        pass
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--dataset', choices=['wes', 'wgs'], required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    args = ap.parse_args()
    args.outdir = args.outdir.resolve()
    args.outdir.mkdir(parents=True, exist_ok=False)
    previous = args.root / (args.dataset + '_combined_candidate/evaluation.json')
    source = json.loads(previous.read_text())
    evidence = args.root / ('wes_indel_attribution' if args.dataset == 'wes' else 'wgs_attribution')
    report = {'status': 'development_not_promoted', 'rules': __doc__,
              'script_sha256': digest(Path(__file__)),
              'helper_sha256': digest(Path(__file__).with_name('explore_indel_corroboration.py')),
              'sources': {str(previous.resolve()): digest(previous)}, 'cells': {}}
    for region in ('ukb', 'medexome'):
        basepath = Path(source['cells'][region]['command'][6])
        base = dict(rows(basepath))
        paths = [evidence / region / (c + '.normalized.vcf.gz') for c in ('deepsomatic', 'mutect2')]
        report['sources'].update({str(p.resolve()): digest(p) for p in [basepath, *paths]})
        ds, m2 = [dict(rows(p)) for p in paths]
        selected = {name: set() for name in ('single_event', 'moderate_strand', 'combined')}
        for key in (ds.keys() & m2.keys()) - base.keys():
            if ',' in key[3] or len(key[2]) == len(key[3]) or key[3].startswith('<') or key[3] == '*':
                continue
            for name in refined_gates(ds[key], m2[key]):
                selected[name].add(key)
                selected['combined'].add(key)
        for name, additions in selected.items():
            dest = args.outdir / region / name
            dest.mkdir(parents=True)
            query = write_query(dest / 'query.vcf', set(base) | additions)
            cmd = list(source['cells'][region]['command'])
            cmd[6], cmd[-1] = query, str(dest / 'benchmark')
            with (dest / 'benchmark.log').open('w') as handle:
                subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=True)
            metrics = parse_metrics_json(dest / 'benchmark.metrics.json')
            report['cells'][region + '/' + name] = {'additions': sorted(additions), 'metrics': metrics, 'command': cmd}
            (args.outdir / 'evaluation.json').write_text(json.dumps(report, indent=2) + '\n')
            print(args.dataset, region, name, len(additions), json.dumps(metrics['indel']), flush=True)
    report['sources_unchanged'] = all(digest(Path(p)) == sha for p, sha in report['sources'].items())
    report['status'] = 'complete_not_promoted' if report['sources_unchanged'] else 'source_changed'
    (args.outdir / 'evaluation.json').write_text(json.dumps(report, indent=2) + '\n')
    if not report['sources_unchanged']:
        raise RuntimeError('Source integrity changed')


if __name__ == '__main__':
    main()
