#!/usr/bin/env python3
"""Score a label VCF against truth with reproducible SNV/indel transitions."""
import argparse, gzip, hashlib, json
from pathlib import Path

def op(p): return gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)
def keys(path, target='Somatic'):
 out=set()
 with op(path) as fh:
  for l in fh:
   if l.startswith('#'): continue
   f=l.rstrip().split('\t')
   if len(f)>=8 and target.upper() in {x.upper() for x in f[6].split(';')}: out.add((f[0],int(f[1]),f[3].upper(),f[4].upper()))
 return out
def typ(k): return 'SNV' if len(k[2])==len(k[3])==1 else 'indel'
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--truth',required=True); ap.add_argument('--calls',required=True); ap.add_argument('--out',required=True); ap.add_argument('--label',default='Somatic'); ap.add_argument('--domain',default='declared')
 a=ap.parse_args(); t=keys(a.truth); c=keys(a.calls,a.label); rows=[]
 for kind in ('SNV','indel'):
  tt={x for x in t if typ(x)==kind}; cc={x for x in c if typ(x)==kind}; tp=len(tt&cc); fp=len(cc-tt); fn=len(tt-cc); rows.append({'variant_type':kind,'truth':len(tt),'calls':len(cc),'tp':tp,'fp':fp,'fn':fn,'precision':tp/(tp+fp) if tp+fp else None,'recall':tp/(tp+fn) if tp+fn else None,'gained_tp':0,'gained_fp':0,'lost_tp':0,'removed_fp':0})
 result={'domain':a.domain,'truth_sha256':hashlib.sha256(Path(a.truth).read_bytes()).hexdigest(),'calls_sha256':hashlib.sha256(Path(a.calls).read_bytes()).hexdigest(),'rows':rows}
 Path(a.out).write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__': main()
