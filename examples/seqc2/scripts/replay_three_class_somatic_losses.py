#!/usr/bin/env python3
"""Replay observed Somatic label losses through the actual CLI, without truth fitting."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import pysam
from validate_refined_native_integration import digest

ROOT = Path(__file__).resolve().parents[3]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--validation', required=True, type=Path)
    ap.add_argument('--transitions', required=True, type=Path)
    ap.add_argument('--outdir', required=True, type=Path)
    args = ap.parse_args()
    validation = json.loads(args.validation.read_text())
    frozen = json.loads(Path(validation['frozen_validation']).read_text())
    original = next(c for c in frozen['commands'] if '--input_dir' in c)
    inputs = Path(original[original.index('--input_dir')+1])
    losses = json.loads(args.transitions.read_text())['somatic_label_losses']
    sites = {tuple(r['allele'][:3])+(tuple(r['allele'][3]),) for r in losses}
    args.outdir.mkdir(parents=True, exist_ok=False)
    staged = args.outdir / 'inputs'
    staged.mkdir()
    regions = args.outdir / 'losses.bed'
    regions.write_text(''.join(f'{c}\t{p-1}\t{p}\n' for c,p in sorted({(s[0],s[1]) for s in sites})))
    report = {'scope':__doc__, 'sources':{}, 'commands':[], 'results':[],
              'classification_sha256':digest(ROOT/'bin/vcf_utils/classification.py')}
    for caller in ('deepsomatic','mutect2','strelka'):
        source = inputs/f'sample.{caller}.vcf.gz'
        report['sources'][str(source)] = digest(source)
        dest = staged/source.name
        cmd = ['bcftools','view','-R',str(regions),str(source),'-Oz','-o',str(dest)]
        report['commands'].append(cmd)
        subprocess.run(cmd,check=True)
        subprocess.run(['bcftools','index','-t',str(dest)],check=True)
    cmd = [sys.executable,str(ROOT/'bin/run_consensus_vcf.py'),'--input_dir',str(staged),
           '--expected_callers','deepsomatic,mutect2,strelka','--experimental-refined-native',
           '--experimental-three-class','--out_prefix',str(args.outdir/'replayed')]
    report['commands'].append(cmd)
    with (args.outdir/'consensus.log').open('x') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    observed = set()
    with pysam.VariantFile(str(args.outdir/'replayed.vcf.gz')) as reader:
        for r in reader:
            key = r.contig,r.pos,r.ref,r.alts
            if key in sites:
                observed.add(key)
                report['results'].append(dict(allele=key,label=next(iter(r.filter),''),rationale=r.info.get('CLASSIFICATION_RATIONALE')))
    report['restored_somatic'] = sum(r['label']=='Somatic' for r in report['results'])
    report['expected_losses'] = len(sites)
    report['missing'] = sorted(sites-observed)
    report['sources_unchanged'] = all(digest(Path(p))==sha for p,sha in report['sources'].items())
    (args.outdir/'replay.json').write_text(json.dumps(report,indent=2)+'\n')
    print({k:report[k] for k in ('restored_somatic','expected_losses','missing','sources_unchanged')})


if __name__ == '__main__':
    main()
