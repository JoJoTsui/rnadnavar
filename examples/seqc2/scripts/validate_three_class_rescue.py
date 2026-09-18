#!/usr/bin/env python3
"""Standalone first/realignment rescue validation from existing VCFs only.

Preparation is the default; --execute runs the two adapters and benchmarks.
Never launches Nextflow, alignment, variant calling or a cohort rerun.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from aggregate_benchmark import parse_metrics_json
from validate_refined_native_integration import digest
from validate_three_class_policy import class_queries

ROOT = Path(__file__).resolve().parents[3]


def plan(validation, out, truth=None):
    if validation.get('status') != 'complete_somatic_truth_screen_not_training_approved':
        raise ValueError('Require completed corrected consensus validation')
    if not validation.get('sources_unchanged') or not validation.get('code_unchanged'):
        raise ValueError('Consensus validation integrity must pass')
    config = validation['manifest']
    consensus_cmd = next(c for c in validation['commands'] if '--out_prefix' in c)
    dna = Path(consensus_cmd[consensus_cmd.index('--out_prefix')+1]+'.vcf.gz')
    jobs = []
    for alignment_round in ('first','realignment'):
        destination = out/alignment_round
        cmd = [sys.executable,str(ROOT/'bin/apply_refined_rescue.py'),
               '--dna-consensus',str(dna),'--annotated-rescue',config['rescues'][alignment_round],
               '--alignment-round',alignment_round,'--experimental-three-class','--outdir',str(destination)]
        for modality,panel in (('dna',config['dna']),('rna',config['rna_'+alignment_round])):
            if set(panel) != {'deepsomatic','mutect2','strelka'}:
                raise ValueError('Incomplete caller panel')
            for caller in ('deepsomatic','mutect2','strelka'):
                cmd += ['--'+modality+'-vcf',caller+'='+panel[caller]]
        jobs.append(dict(alignment_round=alignment_round,command=cmd,outdir=str(destination)))
    return jobs, str(truth or config['truth']), dna


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--consensus-validation',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--truth',type=Path,help='Explicit truth override; recommended tumorvariants for HG008')
    ap.add_argument('--execute',action='store_true')
    args = ap.parse_args()
    validation = json.loads(args.consensus_validation.read_text())
    jobs, truth, dna = plan(validation,args.outdir.resolve(),args.truth)
    config = validation['manifest']
    paths = [dna,Path(truth),Path(config['fasta']),Path(config['hc']),args.consensus_validation]
    for key in ('dna','rna_first','rna_realignment','rescues','targets'):
        paths += [Path(p) for p in config[key].values()]
    for p in paths:
        if not p.is_file(): raise ValueError('Missing source: '+str(p))
    if digest(dna) != validation['output_sha256']:
        raise ValueError('Consensus does not match completed validation')
    # Reject pre-fix consensus provenance even if its structural screen passed.
    classifier = str(ROOT/'bin/vcf_utils/classification.py')
    if validation['code'].get(classifier) != digest(Path(classifier)):
        raise ValueError('Consensus classifier differs from current code; validate current consensus first')
    args.outdir.mkdir(parents=True,exist_ok=False)
    report = {'status':'prepared_not_executed','scope':__doc__,'jobs':jobs,'truth':truth,
              'training_approved':False,'sources':{},'rounds':{},'commands':[],
              'consensus_validation':str(args.consensus_validation.resolve())}

    def save():
        (args.outdir/'validation.json').write_text(json.dumps(report,indent=2)+'\n')

    save()
    if not args.execute:
        print(json.dumps(report,indent=2))
        return
    try:
        report['sources'] = {str(p.resolve()):digest(p) for p in set(paths)}
        report['status'] = 'running'
        save()
        for job in jobs:
            round_name = job['alignment_round']
            report['commands'].append(job['command'])
            save()
            with (args.outdir/(round_name+'.log')).open('x') as log:
                subprocess.run(job['command'],stdout=log,stderr=subprocess.STDOUT,check=True)
            dest = Path(job['outdir'])
            adapter = json.loads((dest/'report.json').read_text())
            if not adapter['sources_unchanged']:
                raise ValueError('Rescue input integrity failed')
            counts,reasons = class_queries(dest/'refined.rescue.vcf.gz',dest)
            metrics = {}
            for label in ('Somatic','Germline','Reference'):
                for domain,target in config['targets'].items():
                    prefix = dest/f'{label}.{domain}'
                    cmd = ['micromamba','run','-n','happy','som.py',truth,str(dest/f'{label}.query.vcf.gz'),
                           '-R',config['hc'],'-T',target,'-r',config['fasta'],'-N','-o',str(prefix)]
                    report['commands'].append(cmd)
                    save()
                    with Path(str(prefix)+'.log').open('x') as log:
                        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
                    metrics[label+'/'+domain] = parse_metrics_json(Path(str(prefix)+'.metrics.json'))
            report['rounds'][round_name] = dict(class_counts=counts,decisions=reasons,metrics=metrics,
                                               adapter=adapter,negative_metric_interpretation='TP means known-somatic collision; not negative-class accuracy')
            save()
        report['sources_unchanged'] = all(digest(Path(p))==sha for p,sha in report['sources'].items())
        if not report['sources_unchanged']: raise ValueError('Source integrity failed')
        report['status'] = 'complete_not_training_approved'
        save()
    except Exception as error:
        report.update(status='failed',error=str(error))
        save()
        raise


if __name__=='__main__':main()
