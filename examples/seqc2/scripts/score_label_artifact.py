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
   if len(f)>=8 and (target is None or target.upper() in {x.upper() for x in f[6].split(';')}): out.add((f[0],int(f[1]),f[3].upper(),f[4].upper()))
 return out
def typ(k): return 'SNV' if len(k[2])==len(k[3])==1 else 'indel'
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--truth',required=True); ap.add_argument('--calls',required=True); ap.add_argument('--out',required=True); ap.add_argument('--label',default='Somatic'); ap.add_argument('--baseline-label',default='PASS'); ap.add_argument('--domain',default='declared'); ap.add_argument('--stage',default='unknown'); ap.add_argument('--baseline'); ap.add_argument('--provenance')
 a=ap.parse_args(); t=keys(a.truth,None); provenance=None
 if a.stage != 'unknown' and not a.provenance: ap.error('--provenance is required when --stage is declared')
 if a.provenance:
  provenance=json.loads(Path(a.provenance).read_text())
  if not isinstance(provenance,dict): ap.error('--provenance must contain a JSON object')
  artifact=provenance.get('artifact')
  if not isinstance(artifact,dict) or not isinstance(artifact.get('sha256'),str): ap.error('--provenance artifact.sha256 is required')
  expected=artifact['sha256']
  actual=hashlib.sha256(Path(a.calls).read_bytes()).hexdigest()
  if expected != actual: ap.error('--provenance artifact sha256 does not match --calls')
 c=keys(a.calls,a.label); baseline=keys(a.baseline,a.baseline_label) if a.baseline else set(); rows=[]
 for kind in ('SNV','indel'):
  tt={x for x in t if typ(x)==kind}; cc={x for x in c if typ(x)==kind}; tp=len(tt&cc); fp=len(cc-tt); fn=len(tt-cc)
  bt={x for x in baseline if typ(x)==kind}; gained=cc-bt; removed=bt-cc
  rows.append({'variant_type':kind,'truth':len(tt),'calls':len(cc),'baseline_calls':len(bt),'tp':tp,'fp':fp,'fn':fn,'precision':tp/(tp+fp) if tp+fp else None,'recall':tp/(tp+fn) if tp+fn else None,'gained_tp':len(gained&tt),'gained_fp':len(gained-tt),'lost_tp':len((bt-cc)&tt),'removed_fp':len((bt-cc)-tt)})
 result={'domain':a.domain,'stage':a.stage,'selectors':{'truth':'all_records','calls':a.label,'baseline':a.baseline_label if a.baseline else None},'truth_sha256':hashlib.sha256(Path(a.truth).read_bytes()).hexdigest(),'calls_sha256':hashlib.sha256(Path(a.calls).read_bytes()).hexdigest(),'baseline_sha256':hashlib.sha256(Path(a.baseline).read_bytes()).hexdigest() if a.baseline else None,'baseline_calls_sha256':hashlib.sha256(Path(a.baseline).read_bytes()).hexdigest() if a.baseline else None,'provenance_sha256':hashlib.sha256(Path(a.provenance).read_bytes()).hexdigest() if a.provenance else None,'provenance':provenance,'rows':rows}
 Path(a.out).write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__': main()
