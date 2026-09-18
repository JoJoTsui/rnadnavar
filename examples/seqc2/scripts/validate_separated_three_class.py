#!/usr/bin/env python3
"""Validate separated candidate labels; never grant biological training approval."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'bin'))
from apply_three_class_labels import run as label, digest
from validate_three_class_policy import class_queries
from validate_frozen_hybrid_policy import pass_query
from aggregate_benchmark import parse_metrics_json


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--native-validation',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--truth',type=Path)
    ap.add_argument('--samplesheet',type=Path,required=True)
    ap.add_argument('--per-class',type=int,default=128)
    args=ap.parse_args()
    native=json.loads(args.native_validation.read_text())
    if native['status']!='complete_somatic_truth_screen_not_training_approved' or not native['sources_unchanged']:
        raise ValueError('Require completed native candidate validation')
    frozen=json.loads(Path(native['frozen_validation']).read_text())
    config=native['manifest']
    truth=str(args.truth.resolve()) if args.truth else config['truth']
    candidates=args.native_validation.parent/'three_class.vcf.gz'
    args.outdir.mkdir(parents=True,exist_ok=False)
    report=dict(status='running',training_approved=False,policy='separated_three_class_v2',
                native_validation=str(args.native_validation.resolve()),truth=truth,
                stages={},commands=[],sources={},code={})
    for p in (Path(truth),Path(config['fasta']),Path(config['hc']),args.samplesheet,
              args.native_validation,Path(native['frozen_validation']),
              *[Path(p) for p in config['targets'].values()]):
        report['sources'][str(p.resolve())]=digest(p)
    for p in (Path(__file__),ROOT/'bin/apply_three_class_labels.py',ROOT/'bin/assess_negative_label_evidence.py',
              ROOT/'examples/seqc2/scripts/pilot_three_class_bam_evidence.py'):
        report['code'][str(p.resolve())]=digest(p)

    def save():
        (args.outdir/'validation.json').write_text(json.dumps(report,indent=2)+'\n')

    def execute(cmd,log):
        report['commands'].append(list(map(str,cmd)));save()
        with log.open('x') as handle:
            subprocess.run(list(map(str,cmd)),stdout=handle,stderr=subprocess.STDOUT,check=True)

    def benchmark(query,prefix,target):
        execute(['micromamba','run','-n','happy','som.py',truth,query,'-R',config['hc'],
                 '-T',target,'-r',config['fasta'],'-N','-o',prefix],Path(str(prefix)+'.log'))
        return parse_metrics_json(Path(str(prefix)+'.metrics.json'))

    try:
        for stage,name in (('consensus','refined_consensus'),('first','refined_first_rescue'),('realignment','refined_realignment_rescue')):
            source=frozen['outputs'][name]
            destination=args.outdir/stage
            outcome=label(Path(source['vcf']),candidates,destination,source['sha256'],native['output_sha256'],stage)
            counts,reasons=class_queries(Path(outcome['output']),destination)
            report['stages'][stage]=dict(adapter=outcome,class_counts=counts,metrics={},baseline_metrics={})
            current=report['stages'][stage]
            baseline_query=destination/'baseline.query.vcf.gz'
            pass_query(Path(source['vcf']),baseline_query,True)
            for domain,target in config['targets'].items():
                original=benchmark(baseline_query,destination/('baseline.'+domain),target)
                current['baseline_metrics'][domain]=original
                for cls in ('Somatic','Germline','Reference'):
                    measured=benchmark(destination/(cls+'.query.vcf.gz'),destination/(cls+'.'+domain),target)
                    current['metrics'][cls+'/'+domain]=measured
                    if cls=='Somatic' and measured!=original:
                        raise ValueError('Somatic metric parity failed: '+stage+'/'+domain)
                    save()
            execute([sys.executable,ROOT/'examples/seqc2/scripts/pilot_three_class_bam_evidence.py',
                     '--vcf',outcome['output'],'--samplesheet',args.samplesheet,'--fasta',config['fasta'],
                     '--per-class',args.per_class,'--out',destination/'bam_pilot.json'],destination/'bam_pilot.log')
            execute([sys.executable,ROOT/'bin/assess_negative_label_evidence.py','--vcf',outcome['output'],
                     '--bam-evidence',destination/'bam_pilot.json','--outdir',destination/'evidence_gate'],destination/'evidence_gate.log')
            current['negative_evidence']=json.loads((destination/'evidence_gate/report.json').read_text())
            current['somatic_parity']=True
            save()
            print(stage,counts,flush=True)
        report['sources_unchanged']=all(digest(Path(p))==sha for p,sha in report['sources'].items())
        report['code_unchanged']=all(digest(Path(p))==sha for p,sha in report['code'].items())
        if not report['sources_unchanged'] or not report['code_unchanged']:
            raise ValueError('Validation integrity failed')
        report['status']='complete_candidate_validation_not_training_approved'
    except Exception as exc:
        report.update(status='failed',error=str(exc));raise
    finally:
        save()


if __name__=='__main__':main()
