#!/usr/bin/env python3
"""VCF-only exploratory gate tests; do not change production labels or defaults.

Keep the entire verified native baseline. For historical gated-rescue additions,
require a DNA variant nomination with positive tumor alternate-read evidence,
not a DeepSomatic RefCall alone. An optional biological veto excludes common
gnomAD alleles and annotated canonical RNA editing without DNA Somatic votes.
Truth is supplied only to som.py after candidate selection.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess

from investigate_verified_rescue import scan
from replay_historical_native_gate import records, write_query
from aggregate_benchmark import parse_metrics_json


def tumor_alt(row, caller, alt):
    samples = [value for name, value in row.get('samples', {}).items()
               if name != 'NORMAL' and not name.endswith('_N_1')]
    if len(samples) != 1:
        return None
    value = samples[0]
    try:
        if caller == 'strelka':
            return int(value[alt + 'U'].split(',')[0])
        return int(value['AD'].split(',')[row['alt_index']])
    except (ValueError, KeyError, IndexError):
        return None


def nominated(panel, key):
    for caller in ('deepsomatic', 'mutect2', 'strelka'):
        for row in panel[caller].get(key, []):
            if caller == 'deepsomatic' and row['filter'] != 'PASS':
                continue
            count = tumor_alt(row, caller, key[3])
            if count is not None and count > 0:
                return True
    return False


def biological_veto(rows):
    for row in rows:
        info = row['info']
        af = info.get('GNOMAD_AF')
        if af not in (None, '.', ''):
            try:
                if max(float(v) for v in str(af).split(',')) > 0.001:
                    return True
            except ValueError:
                raise ValueError(f'Malformed GNOMAD_AF: {af}')
        editing = info.get('REDI_ACCESSION') not in (None, '.', '') and info.get('REDI_CANONICAL') == 'YES'
        if editing and info.get('N_DNA_CALLERS_SOMATIC') == '0':
            return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=False)
    repo = Path(__file__).resolve().parents[3]
    bundle = repo / 'examples/seqc2/verified/20260914'
    provenance = json.loads((bundle / 'comparison/provenance.json').read_text())
    report = {'status': 'exploratory_not_promoted', 'sources': {}, 'selection': {}, 'commands': [], 'metrics': {}}
    jobs = []
    def track(path):
        p = path.resolve(strict=True)
        report['sources'][str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
        return p
    for ds in ('wes_ll','wgs_il'):
        root = bundle / ds
        native = set(records(track(root / 'native_consensus.vcf.gz')))
        original_gate = set(records(track(root / 'native_gated_rescue.vcf.gz')))
        assert native <= original_gate
        additions = original_gate - native
        assert all(len(k[2]) == len(k[3]) == 1 for k in additions)
        panel = {caller: scan(track(root / f'dna_callers/{caller}.vcf.gz'), additions)
                 for caller in ('deepsomatic','mutect2','strelka')}
        pattern = '*.filtered.vcf.gz' if ds == 'wes_ll' else '*.rescue.filtered.stripped.vep.vcf.gz'
        rescue_paths = list((root / 'original_workflow_rescue').glob('*/'+pattern))
        assert len(rescue_paths) == 1
        rescue = scan(track(rescue_paths[0]), additions)
        assert set(rescue) == additions
        for policy in ('nomination', 'nomination_biological'):
            kept = {k for k in additions if nominated(panel, k)
                    and (policy == 'nomination' or not biological_veto(rescue[k]))}
            folder = args.outdir / ds / policy
            folder.mkdir(parents=True)
            query = write_query(folder / 'query.vcf', native | kept)
            report['selection'][f'{ds}/{policy}'] = {
                'baseline_count': len(native), 'baseline_removed': 0,
                'historical_additions': len(additions), 'kept_additions': len(kept),
                'removed_keys': [list(k) for k in sorted(additions-kept)]}
            for domain in ('ukb','medexome'):
                output = folder / domain
                output.mkdir()
                cmd = list(next(c for c in provenance['commands'] if c[-1].endswith(f'{ds}/{domain}/deepsomatic')))
                cmd[6] = query
                cmd[-1] = str(output / 'benchmark')
                report['commands'].append(cmd)
                jobs.append((ds,policy,domain,cmd,output))
    dest = args.outdir / 'evaluation.json'
    dest.write_text(json.dumps(report,indent=2)+'\n')
    def run(job):
        ds,policy,domain,cmd,output = job
        with (output/'benchmark.log').open('w') as log:
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        metrics = parse_metrics_json(output/'benchmark.metrics.json')
        key=f'{ds}/{policy}/{domain}'
        print(key, metrics['records'],flush=True)
        return key,metrics
    with ThreadPoolExecutor(max_workers=2) as pool:
        report['metrics'] = dict(pool.map(run,jobs))
    dest.write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    main()
