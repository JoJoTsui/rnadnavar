#!/usr/bin/env python3
"""Select and freeze a label policy from declared development metrics."""
import argparse, json
from pathlib import Path

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--metrics',required=True); ap.add_argument('--out',required=True); ap.add_argument('--min-delta',type=float,default=0.0)
 a=ap.parse_args(); data=json.loads(Path(a.metrics).read_text()); baseline=data.get('baseline',{}); candidates=data.get('candidates',[]); decisions=[]
 for cand in candidates:
  rows=cand.get('rows',[]); by={(r['pair'],r['variant_type']):r for r in rows}; failures=[]
  for key, b in baseline.items():
   pair, kind=key.split(':',1); r=by.get((pair,kind))
   if not r: failures.append({'key':key,'reason':'missing_candidate_slice'}); continue
   if r.get('precision') is None or r.get('recall') is None: failures.append({'key':key,'reason':'incomplete_metrics'}); continue
   if r['precision'] < b['precision'] + a.min_delta: failures.append({'key':key,'reason':'precision_gate','candidate':r['precision'],'baseline':b['precision']})
   if r['recall'] < b['recall'] + a.min_delta: failures.append({'key':key,'reason':'recall_gate','candidate':r['recall'],'baseline':b['recall']})
  decisions.append({'policy':cand.get('policy'),'qualifies':not failures,'failures':failures})
 qualified=[d for d in decisions if d['qualifies']]
 result={'status':'qualified' if qualified else 'no_qualifying_policy','development_only':True,'min_delta':a.min_delta,'decisions':decisions,'frozen_policy':qualified[0]['policy'] if qualified else None}
 Path(a.out).write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__': main()
