#!/usr/bin/env python3
"""Verify RNA nominated alleles against caller-ready DNA tumor/normal BAMs.

The command is deliberately evidence-only: it emits confirmed, rejected, or
inconclusive outcomes and never writes a Somatic FILTER.  Rules are explicit
and can be swept by changing the JSON configuration.
"""
import argparse, json, os, shutil, subprocess, sys
from pathlib import Path

def parse_sites(path):
    opener = open
    if str(path).endswith('.gz'):
        import gzip; opener = gzip.open
    rows=[]
    with opener(path, 'rt') as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'): continue
            f=line.rstrip().split('	')
            if len(f) >= 5 and f[1].isdigit(): rows.append((f[0],int(f[1]),f[3].upper(),f[4].upper()))
    return rows

def pileup(samtools,bam,ref,site,min_mapq,min_baseq):
    c,pos,r,a=site
    cmd=[samtools,'mpileup','-f',ref,'-r',f'{c}:{pos}-{pos}','-q',str(min_mapq),'-Q',str(min_baseq),bam]
    p=subprocess.run(cmd,text=True,capture_output=True)
    if p.returncode: raise RuntimeError(p.stderr.strip() or f'mpileup failed ({p.returncode})')
    for line in p.stdout.splitlines():
        f=line.split('	')
        if len(f)>=5:
            return int(f[3]), f[4]
    return 0,''

def counts(depth,bases,ref,alt):
    if not bases or depth<=0: return 0,depth
    i=0; n=0
    while i<len(bases):
        c=bases[i]
        if c=='^': i+=2; continue
        if c=='$': i+=1; continue
        if c in '+-':
            j=i+1
            while j<len(bases) and bases[j].isdigit(): j+=1
            ln=int(bases[i+1:j] or 0); seq=bases[j:j+ln].upper(); i=j+ln
            if c=='+' and alt.startswith(ref) and alt[len(ref):]==seq: n+=1
            continue
        if c.upper()==alt: n+=1
        i+=1
    return n,depth

def main():
    ap=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--candidates',required=True,help='TSV/VCF with CHROM POS REF ALT columns')
    ap.add_argument('--reference',required=True); ap.add_argument('--tumor-bam',required=True); ap.add_argument('--normal-bam',required=True)
    ap.add_argument('--out',required=True); ap.add_argument('--min-tumor-alt',type=int,default=3); ap.add_argument('--max-normal-alt',type=int,default=0); ap.add_argument('--min-dp',type=int,default=10); ap.add_argument('--samtools',default='samtools'); ap.add_argument('--min-mapq',type=int,default=20); ap.add_argument('--min-baseq',type=int,default=20)
    args=ap.parse_args(); sam=shutil.which(args.samtools) or (args.samtools if os.path.isfile(args.samtools) else None)
    if not sam: ap.error(f'samtools not found: {args.samtools}')
    for x in (args.reference,args.tumor_bam,args.normal_bam,args.candidates):
        if not os.path.isfile(x): ap.error(f'input not found: {x}')
    rules={'min_tumor_alt':args.min_tumor_alt,'max_normal_alt':args.max_normal_alt,'min_dp':args.min_dp,'min_mapq':args.min_mapq,'min_baseq':args.min_baseq}
    if any(',' in site[3] or site[3] in {'.','*'} for site in parse_sites(args.candidates)):
        raise SystemExit('ambiguous or symbolic ALT is not supported; split candidates before verification')
    out=[]
    for site in parse_sites(args.candidates):
        row={'chrom':site[0],'pos':site[1],'ref':site[2],'alt':site[3],'status':'inconclusive','tumor_dp':None,'tumor_alt':None,'normal_dp':None,'normal_alt':None,'reason':''}
        try:
            td,tb=pileup(sam,args.tumor_bam,args.reference,site,args.min_mapq,args.min_baseq); nd,nb=pileup(sam,args.normal_bam,args.reference,site,args.min_mapq,args.min_baseq)
            ta,_=counts(td,tb,site[2],site[3]); na,_=counts(nd,nb,site[2],site[3]); row.update(tumor_dp=td,tumor_alt=ta,normal_dp=nd,normal_alt=na)
            if td < args.min_dp or nd < args.min_dp: row['reason']='insufficient_tumor_or_normal_coverage'
            elif ta >= args.min_tumor_alt and na <= args.max_normal_alt: row.update(status='confirmed',reason='tumor_support_and_normal_clean')
            elif na > args.max_normal_alt: row.update(status='rejected',reason='matched_normal_support')
            else: row['reason']='insufficient_tumor_alt_support'
        except Exception as e: row['reason']=f'verification_error:{e}'
        out.append(row)
    Path(args.out).write_text(json.dumps({'rules':rules,'results':out},indent=2)+'\n')
    print(json.dumps({'confirmed':sum(r['status']=='confirmed' for r in out),'rejected':sum(r['status']=='rejected' for r in out),'inconclusive':sum(r['status']=='inconclusive' for r in out)}))
if __name__=='__main__': main()
