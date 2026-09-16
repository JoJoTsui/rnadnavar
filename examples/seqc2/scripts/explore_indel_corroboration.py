#!/usr/bin/env python3
"""Bounded SEQC2 indel ablations on the fixed 20260916 consensus candidate.

No truth is read during selection. Three additions are tested separately and
jointly: Mutect2 PASS corroborated by DeepSomatic tumor reads; Strelka PASS
corroborated by Mutect2; lower-QUAL DeepSomatic PASS with Mutect2 bidirectional
alternate strand support. All require adequate clean Mutect2 normal evidence.
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


def sample(record, suffix):
    found = [v for k, v in record.get('samples', {}).items() if k.endswith(suffix)]
    return found[0] if len(found) == 1 else {}


def number(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('Non-finite evidence')
    return value


def gates(ds, m2, st):
    """Return eligible hypotheses from sample-resolved evidence; missing fails."""
    passing = {'PASS', '.'}
    result = set()
    try:
        info = dict(v.split('=', 1) for v in m2['info'].split(';') if '=' in v)
        tumor, normal = sample(m2, '_T_1'), sample(m2, '_N_1')
        ma = int(tumor['AD'].split(',')[1])
        if (int(normal['DP']) < 10 or int(normal['AD'].split(',')[1]) != 0
                or number(info['GERMQ']) < 20 or number(info['TLOD']) <= 0):
            return result
        soft = (set(m2['filter'].split(';')) - passing) <= {'contamination', 'weak_evidence'}
        if not soft or ma < 3:
            return result
    except (KeyError, ValueError, IndexError, TypeError):
        return result
    try:
        da = int(sample(ds, '_T_1')['AD'].split(',')[1])
        quality = number(ds['qual'])
        if ds['filter'] in passing | {'RefCall'} and quality > 0 and da >= 3:
            if m2['filter'] in passing:
                result.add('mutect_pass')
            strand = [int(x) for x in tumor['SB'].split(',')]
            if ds['filter'] in passing and len(strand) == 4 and min(strand[2:]) >= 1:
                result.add('ds_strand')
    except (KeyError, ValueError, IndexError, TypeError):
        pass
    try:
        if (st['filter'] in passing and int(st['samples']['TUMOR']['TIR'].split(',')[0]) >= 3
                and int(st['samples']['NORMAL']['TIR'].split(',')[0]) == 0
                and int(st['samples']['NORMAL']['DP']) >= 10):
            result.add('strelka_pass')
    except (KeyError, ValueError, IndexError, TypeError):
        pass
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, required=True, help='Completed current_native_audit root')
    ap.add_argument('--dataset', choices=['wes', 'wgs'], required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    args = ap.parse_args()
    args.outdir = args.outdir.resolve()
    args.outdir.mkdir(parents=True, exist_ok=False)
    previous = args.root / (args.dataset + '_combined_candidate/evaluation.json')
    source = json.loads(previous.read_text())
    evidence = args.root / ('wes_indel_attribution' if args.dataset == 'wes' else 'wgs_attribution')
    report = {'status': 'development_not_promoted', 'rules': __doc__,
              'script_sha256': digest(Path(__file__)), 'sources': {str(previous.resolve()): digest(previous)},
              'cells': {}}
    for region in ('ukb', 'medexome'):
        basepath = Path(source['cells'][region]['command'][6])
        base = dict(rows(basepath))
        paths = [evidence / region / (c + '.normalized.vcf.gz') for c in ('deepsomatic', 'mutect2', 'strelka')]
        report['sources'].update({str(p.resolve()): digest(p) for p in [basepath, *paths]})
        ds, m2, st = [dict(rows(p)) for p in paths]
        selected = {name: set() for name in ('mutect_pass', 'strelka_pass', 'ds_strand', 'combined')}
        for key in (ds.keys() | m2.keys() | st.keys()) - base.keys():
            if ',' in key[3] or len(key[2]) == len(key[3]) or key[3].startswith('<') or key[3] == '*':
                continue
            for name in gates(ds.get(key, {}), m2.get(key, {}), st.get(key, {})):
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
