#!/usr/bin/env python3
"""Extract region-restricted rescue additions and their original caller evidence.

All additions here must be biallelic SNVs. Exact allele truth attribution is
cross-checked against the verified som.py count deltas; it is not an indel or
haplotype comparison. Missing records/fields remain missing, never zero support.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

from replay_historical_native_gate import records


def scan(path, wanted):
    found = {}
    samples = []
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                samples = line.rstrip().split('\t')[9:]
            if line.startswith('#'):
                continue
            f = line.rstrip().split('\t')
            for index, alt in enumerate(f[4].split(','), 1):
                key = (f[0], int(f[1]), f[3], alt)
                if key not in wanted:
                    continue
                row = {'filter': f[6], 'qual': f[5], 'alt_index': index,
                       'info': dict(x.split('=', 1) if '=' in x else (x, True)
                                    for x in f[7].split(';') if x != '.'),
                       'samples': {sample: dict(zip(f[8].split(':'), values.split(':')))
                                   for sample, values in zip(samples, f[9:])}}
                found.setdefault(key, []).append(row)
    return found


def regional(path, hc, target):
    output = subprocess.check_output(['bcftools', 'view', '-H', '-R', str(hc),
                                      '-T', str(target), str(path)], text=True)
    return {(f[0], int(f[1]), f[3], f[4]) for line in output.splitlines()
            if len(f := line.split('\t')) >= 8}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[3]
    bundle = repo / 'examples/seqc2/verified/20260914'
    comparison = bundle / 'comparison'
    provenance = json.loads((comparison / 'provenance.json').read_text())
    hc = next(Path(p) for p in provenance['sources'] if 'High-Confidence_Regions' in p)
    targets = {'ukb': next(Path(p) for p in provenance['sources'] if 'ukb.pad' in p),
               'medexome': next(Path(p) for p in provenance['sources'] if 'MedExome' in p)}
    truth_path = next(Path(p) for p in provenance['sources'] if p.endswith('benchmark_truth.vcf.gz'))
    truth = set(records(truth_path))
    results = []
    sources = {str(truth_path): hashlib.sha256(truth_path.read_bytes()).hexdigest()}
    for ds in ('wes_ll', 'wgs_il'):
        root = bundle / ds
        domains = {}
        for domain, bed in targets.items():
            before = regional(root / 'native_consensus.vcf.gz', hc, bed)
            after = regional(root / 'native_gated_rescue.vcf.gz', hc, bed)
            if before - after:
                raise ValueError(f'Unexpected baseline loss: {ds}/{domain}')
            domains[domain] = after - before
        wanted = set.union(*domains.values())
        if any(len(k[2]) != 1 or len(k[3]) != 1 for k in wanted):
            raise ValueError('This evidence extractor is restricted to SNV additions')
        panel = {}
        for modality in ('dna', 'rna', 'rna_realign'):
            for caller in ('deepsomatic', 'mutect2', 'strelka'):
                path = (root / f'{modality}_callers/{caller}.vcf.gz').resolve(strict=True)
                sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
                panel[f'{modality}/{caller}'] = scan(path, wanted)
        rescue_root = root / 'original_workflow_rescue'
        suffix = '*.filtered.vcf.gz' if ds == 'wes_ll' else '*.rescue.filtered.stripped.vep.vcf.gz'
        paths = list(rescue_root.glob('*/' + suffix))
        if len(paths) != 1:
            raise ValueError(f'Ambiguous source rescue: {paths}')
        rescue = scan(paths[0], wanted)
        sources[str(paths[0].resolve())] = hashlib.sha256(paths[0].read_bytes()).hexdigest()
        for key in sorted(wanted):
            results.append({'dataset': ds, 'chrom': key[0], 'pos': key[1], 'ref': key[2], 'alt': key[3],
                            'truth_status': 'TP' if key in truth else 'FP',
                            'domains': [name for name, keys in domains.items() if key in keys],
                            'rescue': rescue.get(key, []),
                            'callers': {name: data.get(key, []) for name, data in panel.items()}})
    summary = {}
    for ds in ('wes_ll', 'wgs_il'):
        for domain in targets:
            rows = [r for r in results if r['dataset'] == ds and domain in r['domains']]
            summary[f'{ds}/{domain}'] = {label: sum(r['truth_status'] == label for r in rows) for label in ('TP', 'FP')}
    expected = {'wes_ll/ukb': {'TP': 11, 'FP': 4}, 'wes_ll/medexome': {'TP': 6, 'FP': 2},
                'wgs_il/ukb': {'TP': 0, 'FP': 3}, 'wgs_il/medexome': {'TP': 0, 'FP': 1}}
    if summary != expected:
        raise ValueError(f'Allele attribution disagrees with verified som.py deltas: {summary}')
    args.outdir.mkdir(parents=True, exist_ok=False)
    (args.outdir / 'evidence.json').write_text(json.dumps({'summary': summary, 'sources': sources, 'sites': results}, indent=2) + '\n')
    with (args.outdir / 'sites.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['dataset','chrom','pos','ref','alt','truth_status','domains'])
        for r in results:
            writer.writerow([r[k] for k in ('dataset','chrom','pos','ref','alt','truth_status')] + [','.join(r['domains'])])
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
