#!/usr/bin/env python3
"""Validate experimental output adapter on existing SEQC2 rescue rounds.

The DNA baseline is benchmark-domain scoped. These are validation artifacts,
not whole-genome training-label releases. Never launches Nextflow or HG008.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from aggregate_benchmark import parse_metrics_json
from audit_current_native_policy import digest
from replay_historical_native_gate import write_query
from validate_refined_native_integration import alleles


def rescue_pattern(dataset, round_name):
    # WGS first rescue is published as filtered.vcf.gz, not the realignment
    # VEP basename. Both include the population/editing annotations used here.
    return ('*/*.rescue.filtered.stripped.vep.vcf.gz'
            if dataset == 'wgs' and round_name == 'realignment'
            else '*/*.filtered.vcf.gz')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--dataset', choices=['wes', 'wgs'], required=True)
    ap.add_argument('--region', choices=['ukb', 'medexome'], required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    args = ap.parse_args()
    repo = Path(__file__).resolve().parents[3]
    dataset = 'wes_ll' if args.dataset == 'wes' else 'wgs_il'
    bundle = repo / 'examples/seqc2/verified/20260914' / dataset
    baseline = args.root / 'driver_parity_v1' / args.dataset / args.region / 'refined.vcf.gz'
    frozen_path = args.root / f'{args.dataset}_refined_gate_replay/evaluation.json'
    frozen = json.loads(frozen_path.read_text())['cells'][args.region]
    args.outdir.mkdir(parents=True, exist_ok=False)
    report = {'status': 'running', 'scope': __doc__, 'dataset': dataset, 'region': args.region,
              'sources': {str(p.resolve()): digest(p) for p in (baseline, frozen_path)}, 'rounds': {}}
    for round_name in ('first', 'realignment'):
        if round_name == 'realignment':
            source_dir = bundle / 'original_workflow_rescue'
            rna_dir = bundle / 'rna_realign_callers'
        else:
            source_dir = bundle / ('workflow_dna_source' if args.dataset == 'wes' else 'workflow_source') / 'rescue'
            rna_dir = bundle / 'rna_callers'
        pattern = rescue_pattern(args.dataset, round_name)
        matches = list(source_dir.glob(pattern))
        if len(matches) != 1:
            raise ValueError(f'Ambiguous annotated {round_name} source: {matches}')
        source = matches[0].resolve(strict=True)
        report['sources'][str(source)] = digest(source)
        dest = args.outdir / round_name
        cmd = [sys.executable, str(repo / 'bin/apply_refined_rescue.py'), '--dna-consensus', str(baseline),
               '--annotated-rescue', str(source), '--alignment-round', round_name, '--outdir', str(dest)]
        for modality, folder in [('dna', bundle / 'dna_callers'), ('rna', rna_dir)]:
            for caller in ('deepsomatic', 'mutect2', 'strelka'):
                path = (folder / f'{caller}.vcf.gz').resolve(strict=True)
                report['sources'][str(path)] = digest(path)
                cmd += [f'--{modality}-vcf', f'{caller}={path}']
        with (args.outdir / f'{round_name}.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        adapter = json.loads((dest / 'report.json').read_text())
        query = write_query(dest / 'benchmark.query.vcf', alleles(dest / 'refined.rescue.vcf.gz', True))
        benchmark = list(frozen['command'])
        benchmark[6], benchmark[-1] = query, str(dest / 'benchmark')
        with (dest / 'benchmark.log').open('w') as log:
            subprocess.run(benchmark, stdout=log, stderr=subprocess.STDOUT, check=True)
        metrics = parse_metrics_json(dest / 'benchmark.metrics.json')
        report['rounds'][round_name] = {'command': cmd, 'benchmark_command': benchmark, 'metrics': metrics,
                                      'adapter_counts': adapter['counts'], 'sources_unchanged': adapter['sources_unchanged'],
                                      'matches_frozen_realignment': metrics == frozen['metrics'] if round_name == 'realignment' else None}
        (args.outdir / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
        print(dataset, args.region, round_name, metrics['records'], flush=True)
    report['sources_unchanged'] = all(digest(Path(p)) == sha for p, sha in report['sources'].items())
    report['status'] = 'complete_experimental_not_training_ready' if report['sources_unchanged'] else 'integrity_failure'
    (args.outdir / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
