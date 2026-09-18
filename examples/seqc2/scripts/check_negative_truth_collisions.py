#!/usr/bin/env python3
"""Challenge negative evidence at known Somatic exact-allele collisions.

This targeted regression panel is truth-selected, unlike the random-hash pilot.
It must never be used as an accuracy estimate or to fit thresholds. Indel
negative evidence remains withheld pending haplotype-aware adjudication.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

import pysam

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'bin'))
from assess_negative_label_evidence import assess, digest
from pilot_three_class_bam_evidence import evidence, FLAG_FILTER


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--vcf',type=Path,required=True)
    ap.add_argument('--truth',type=Path,required=True)
    ap.add_argument('--samplesheet',type=Path,required=True)
    ap.add_argument('--fasta',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():raise ValueError('Refuse overwrite')
    paths=(args.vcf,args.truth,args.samplesheet,args.fasta)
    sources={str(p.resolve()):digest(p) for p in paths}
    with pysam.VariantFile(str(args.truth)) as f:
        truth={(r.contig,r.pos,r.ref,a) for r in f for a in r.alts or []}
    with args.samplesheet.open() as f: samples=list(csv.DictReader(f))
    pair={}
    for status in ('0','1'):
        matches=[s for s in samples if s['status']==status]
        if len(matches)!=1:raise ValueError('Require one DNA pair')
        pair[status]=matches[0]
    report=dict(scope=__doc__,sources=sources,results=[],training_approved=False,
                bam_paths={k:v['bam'] for k,v in pair.items()},flag_filter=FLAG_FILTER)
    with pysam.AlignmentFile(pair['0']['bam'],index_filename=pair['0']['bai']) as normal, \
         pysam.AlignmentFile(pair['1']['bam'],index_filename=pair['1']['bai']) as tumor, \
         pysam.FastaFile(str(args.fasta)) as fasta, pysam.VariantFile(str(args.vcf)) as reader:
        for record in reader:
            label=next(iter(record.filter),'')
            key=(record.contig,record.pos,record.ref,','.join(record.alts or []))
            if label not in {'Germline','Reference'} or key not in truth:continue
            n=t=None
            if len(key[2])==len(key[3])==1 and not set(key[2]+key[3])-set('ACGT'):
                n,t=evidence(normal,fasta,key),evidence(tumor,fasta,key)
            status,reason=assess(label,key[2],key[3],n,t)
            report['results'].append(dict(site=key,label=label,normal=n,tumor=t,status=status,reason=reason))
    report['counts']=dict(Counter(r['label']+':'+r['status']+':'+r['reason'] for r in report['results']))
    report['supported_known_somatic_collisions']=sum(r['status']=='SUPPORTED' for r in report['results'])
    report['sources_unchanged']=all(digest(Path(p))==sha for p,sha in sources.items())
    report['status']='pass_targeted_screen_not_accuracy' if report['supported_known_somatic_collisions']==0 and report['sources_unchanged'] else 'failed'
    with args.out.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(report['counts']))
    if report['status']=='failed':raise SystemExit(1)


if __name__=='__main__':main()
