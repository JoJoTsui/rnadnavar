#!/usr/bin/env python3
"""Audit fixed rescue-gate exclusions by truth confidence and target domain.

No candidate selection, workflow execution, or production-label modification.
Outside-HC records are explicitly unassessed, never labelled false positive.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess

from investigate_verified_rescue import scan, regional
from replay_historical_native_gate import records


def hc_records(path, hc):
    text = subprocess.check_output(['bcftools', 'view', '-H', '-R', str(hc), str(path)], text=True)
    return {(f[0], int(f[1]), f[3], f[4]) for line in text.splitlines()
            if len(f := line.split('\t')) >= 8}


def reasons(rows):
    result = set()
    for row in rows:
        info = row['info']
        af = info.get('GNOMAD_AF')
        if af not in (None, '.', '') and max(map(float, str(af).split(','))) > 0.001:
            result.add('common_af')
        if (info.get('REDI_ACCESSION') not in (None, '.', '')
                and info.get('REDI_CANONICAL') == 'YES'
                and info.get('N_DNA_CALLERS_SOMATIC') == '0'):
            result.add('editing_without_dna_somatic')
    return result


def domain(key, hc, ukb, medexome):
    if key not in hc:
        return 'outside_hc_unassessed'
    return ('hc_ukb' if key in ukb else 'hc_outside_ukb') + (
        '_medexome' if key in medexome else '_outside_medexome')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[3]
    bundle = repo / 'examples/seqc2/verified/20260914'
    experiment = repo / 'examples/seqc2/comparison/rescue_fp_investigation_20260914/gate_tests'
    evaluation = json.loads((experiment / 'evaluation.json').read_text())
    provenance = json.loads((bundle / 'comparison/provenance.json').read_text())
    command = provenance['commands'][0]
    hc = Path(command[command.index('-R') + 1])
    truth_path = Path(command[5])
    truth = set(records(truth_path))
    report = {'status': 'diagnostic_not_independent_validation', 'sources': {}, 'summary': {}, 'sites': []}
    for path in (hc, truth_path, experiment / 'evaluation.json'):
        report['sources'][str(path.resolve())] = hashlib.sha256(path.read_bytes()).hexdigest()
    for ds in ('wes_ll', 'wgs_il'):
        root = bundle / ds
        gate = root / 'native_gated_rescue.vcf.gz'
        base = set(records(root / 'native_consensus.vcf.gz'))
        additions = set(records(gate)) - base
        hc_set = hc_records(gate, hc)
        targets = {}
        for target in ('ukb', 'medexome'):
            cmd = next(c for c in provenance['commands'] if c[-1].endswith(f'{ds}/{target}/deepsomatic'))
            bed = Path(cmd[cmd.index('-T') + 1])
            report['sources'][str(bed.resolve())] = hashlib.sha256(bed.read_bytes()).hexdigest()
            targets[target] = regional(gate, hc, bed)
        pattern = '*.filtered.vcf.gz' if ds == 'wes_ll' else '*.rescue.filtered.stripped.vep.vcf.gz'
        paths = list((root / 'original_workflow_rescue').glob('*/' + pattern))
        assert len(paths) == 1
        rescue = scan(paths[0], additions)
        assert set(rescue) == additions
        for path in (gate, root / 'native_consensus.vcf.gz', paths[0]):
            report['sources'][str(path.resolve())] = hashlib.sha256(path.read_bytes()).hexdigest()
        nomination_removed = {tuple(k) for k in evaluation['selection'][f'{ds}/nomination']['removed_keys']}
        removed = {tuple(k) for k in evaluation['selection'][f'{ds}/nomination_biological']['removed_keys']}
        counts = Counter()
        for key in sorted(additions):
            why = reasons(rescue[key])
            if key in nomination_removed:
                why.add('no_dna_nomination')
            assert bool(why) == (key in removed), key
            region = domain(key, hc_set, targets['ukb'], targets['medexome'])
            truth_status = ('truth_present' if key in truth else 'truth_absent') if key in hc_set else 'unassessed'
            decision = 'excluded' if key in removed else 'retained'
            signature = '+'.join(sorted(why)) or 'none'
            counts[(region, decision, truth_status, signature)] += 1
            report['sites'].append({'dataset': ds, 'allele': key, 'domain': region,
                                    'decision': decision, 'truth_status': truth_status,
                                    'reasons': sorted(why)})
        report['summary'][ds] = [dict(zip(('domain', 'decision', 'truth_status', 'reasons'), key), count=count)
                                 for key, count in sorted(counts.items())]
    args.outdir.mkdir(parents=True, exist_ok=False)
    (args.outdir / 'domain_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['summary'], indent=2))


if __name__ == '__main__':
    main()
