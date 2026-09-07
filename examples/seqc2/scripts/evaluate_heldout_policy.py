#!/usr/bin/env python3
"""Evaluate a frozen policy on held-out slices without retuning."""
import argparse, json
from pathlib import Path

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--frozen',required=True); ap.add_argument('--metrics',required=True); ap.add_argument('--out',required=True); a=ap.parse_args()
 frozen=json.loads(Path(a.frozen).read_text()); held=json.loads(Path(a.metrics).read_text());
 if frozen.get('status')!='qualified' or not frozen.get('frozen_policy'): status='incomplete_no_frozen_policy'; rows=[]
 else:
  rows=held.get('rows',[]); status='complete' if rows else 'incomplete_missing_metrics'
  for r in rows:
   r['accepted']=bool(r.get('precision') is not None and r.get('recall') is not None and r.get('precision_gate',True) and r.get('recall_gate',True))
 result={'status':status,'rules_source':'frozen_development_policy','retuning_allowed':False,'rows':rows,'limitations':held.get('limitations',[])}
 Path(a.out).write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__': main()
