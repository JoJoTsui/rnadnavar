#!/usr/bin/env python3
"""Export normalized consensus/DeepSomatic indel disagreements and caller evidence.

Diagnostic only: truth labels never enter candidate selection. Uses existing
audit commands and retains som.py scratch for independently inspectable labels.
"""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import subprocess

from aggregate_benchmark import parse_metrics_json


def rows(path):
    with gzip.open(path, 'rt') as handle:
        samples = []
        for line in handle:
            if line.startswith('#CHROM'):
                samples = line.rstrip().split('\t')[9:]
            if line.startswith('#'):
                continue
            row = line.rstrip().split('\t')
            key = (row[0], int(row[1]), row[3], row[4])
            fmt = row[8].split(':') if len(row) > 8 else []
            yield key, {'filter': row[6], 'qual': row[5], 'info': row[7],
                        'samples': {s: dict(zip(fmt, v.split(':')))
                                    for s, v in zip(samples, row[9:])}}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--audit', type=Path, required=True)
    ap.add_argument('--dataset', choices=['wes_ll', 'wgs_il'], required=True)
    ap.add_argument('--variant-type', choices=['indel', 'snp', 'all'], default='indel')
    ap.add_argument('--outdir', type=Path, required=True)
    args = ap.parse_args()
    doc = json.loads(args.audit.read_text())
    args.outdir = args.outdir.resolve()
    args.outdir.mkdir(parents=True, exist_ok=False)
    dataset = args.dataset
    results = {}
    commands = []
    def run(cmd, log):
        commands.append(cmd)
        (args.outdir / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
        with log.open('w') as handle:
            subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=True)
    for region in ('ukb', 'medexome'):
        cell = args.outdir / region
        cell.mkdir()
        queries = {}
        for method in ('current_native', 'deepsomatic'):
            suffix = f'/{dataset}/{region}/{method}'
            cmd = next(list(c) for c in doc['commands'] if c[-1].endswith(suffix))
            original_output = cmd[-1]
            scratch = cell / (method + '_scratch')
            cmd[-1] = str(cell / method)
            cmd += ['--scratch-prefix', str(scratch)]
            run(cmd, cell / (method + '.log'))
            measured = parse_metrics_json(cell / (method + '.metrics.json'))
            if measured != doc['metrics'][f'{dataset}/{region}/{method}']:
                raise RuntimeError(f'Benchmark parity failed for {original_output}')
            queries[method] = dict(rows(scratch / 'normalized_query.vcf.gz'))
        truth = set(k for k, _ in rows(cell / 'current_native_scratch/normalized_truth.vcf.gz'))
        current, baseline = queries['current_native'], queries['deepsomatic']
        changed = {k for k in set(current) ^ set(baseline)
                   if args.variant_type == 'all'
                   or (args.variant_type == 'indel' and len(k[2]) != len(k[3]))
                   or (args.variant_type == 'snp' and len(k[2]) == len(k[3]) == 1)}
        evidence = {k: {'key': k, 'direction': 'added' if k in current else 'lost',
                        'truth': 'TP' if k in truth else 'FP', 'callers': {}} for k in changed}
        repo = Path(__file__).resolve().parents[3]
        source_dir = repo / 'examples/seqc2/verified/20260914' / dataset / 'dna_callers'
        for caller in ('deepsomatic', 'mutect2', 'strelka'):
            selected = cell / (caller + '.selected.vcf.gz')
            normalized = cell / (caller + '.normalized.vcf.gz')
            run(['bcftools', 'view', str(source_dir / (caller + '.vcf.gz')),
                 '-R', cmd[cmd.index('-R') + 1], '-T', cmd[cmd.index('-T') + 1],
                 '-Oz', '-o', str(selected)], cell / (caller + '.select.log'))
            run(['bcftools', 'norm', '-f', doc['reference'], '-c', 'x', '-d', 'exact',
                 str(selected), '-Oz', '-o', str(normalized)], cell / (caller + '.norm.log'))
            for key, value in rows(normalized):
                if key in evidence:
                    evidence[key]['callers'].setdefault(caller, []).append(value)
        records = [evidence[k] for k in sorted(evidence)]
        counts = dict(Counter(r['direction'] + '_' + r['truth'] for r in records))
        results[region] = {'counts': counts, 'records': records}
        print(dataset, region, counts, flush=True)
    (args.outdir / (args.variant_type + '_attribution.json')).write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
