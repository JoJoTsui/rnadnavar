#!/usr/bin/env python3
"""Transfer preserved SEQC2 realignment-gate additions to a new consensus query.

This is a historical gate replay, not validation of a portable rescue driver.
The frozen candidate scope and annotations are inherited and disclosed.
"""
import argparse
import json
from pathlib import Path
import subprocess

from aggregate_benchmark import parse_metrics_json
from audit_current_native_policy import alleles, digest
from replay_historical_native_gate import write_query


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--candidate', type=Path, required=True)
    ap.add_argument('--dataset', choices=['wes_ll', 'wgs_il'], required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    ap.add_argument('--policy', help='Select region/POLICY cells from an ablation report')
    args = ap.parse_args()
    repo = Path(__file__).resolve().parents[3]
    root = repo / 'examples/seqc2'
    historical = root / 'verified/20260914' / args.dataset / 'native_consensus.vcf.gz'
    gated = root / 'comparison/rescue_fp_investigation_20260914/gate_tests' / args.dataset / 'nomination_biological/query.vcf.gz'
    accepted = {'Somatic', 'PASS', '.'}
    additions = alleles(gated, accepted) - alleles(historical, accepted)
    if any(len(k[2]) != 1 or len(k[3]) != 1 for k in additions):
        raise ValueError('Frozen rescue additions must be SNPs')
    candidate = json.loads(args.candidate.read_text())
    cells = ({region: candidate['cells'][region + '/' + args.policy] for region in ('ukb', 'medexome')}
             if args.policy else candidate['cells'])
    if set(cells) != {'ukb', 'medexome'}:
        raise ValueError('Expected exactly UKB and MedExome cells')
    args.outdir = args.outdir.resolve()
    args.outdir.mkdir(parents=True, exist_ok=False)
    report = {'status': 'historical_gate_transfer_not_promoted', 'scope': __doc__,
              'script_sha256': digest(Path(__file__)),
              'sources': {str(p.resolve()): digest(p) for p in (historical, gated, args.candidate)},
              'frozen_additions': sorted(additions), 'cells': {}}
    for region, cell in cells.items():
        cmd = list(cell['command'])
        source = Path(cmd[6])
        report['sources'][str(source)] = digest(source)
        native = alleles(source, accepted)
        dest = args.outdir / region
        dest.mkdir()
        cmd[6] = write_query(dest / 'query.vcf', native | additions)
        cmd[-1] = str(dest / 'benchmark')
        with (dest / 'benchmark.log').open('w') as handle:
            subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=True)
        metrics = parse_metrics_json(dest / 'benchmark.metrics.json')
        report['cells'][region] = {'command': cmd, 'metrics': metrics}
        (args.outdir / 'evaluation.json').write_text(json.dumps(report, indent=2) + '\n')
        print(args.dataset, region, json.dumps(metrics), flush=True)


if __name__ == '__main__':
    main()
