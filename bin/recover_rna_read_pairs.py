#!/usr/bin/env python3
"""Recover candidate RNA pairs from original FASTQ files without cross-library mixing."""
import argparse, gzip, json
from pathlib import Path

def op(path): return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)
def rid(name): return name.split()[0].rstrip('/1').rstrip('/2')
def read_fastq(path):
 with op(path) as fh:
  while True:
   h=fh.readline()
   if not h: break
   seq,plus,qual=fh.readline(),fh.readline(),fh.readline()
   if not qual: break
   yield rid(h[1:].strip()), h+seq+plus+qual

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--read-ids',required=True); ap.add_argument('--fastq-1',required=True); ap.add_argument('--fastq-2',required=True); ap.add_argument('--library',required=True); ap.add_argument('--out-prefix',required=True); ap.add_argument('--stats',required=True); a=ap.parse_args()
 ids={x.strip().split()[0] for x in Path(a.read_ids).read_text().splitlines() if x.strip()}; r1={k:v for k,v in read_fastq(a.fastq_1) if k in ids}; r2={k:v for k,v in read_fastq(a.fastq_2) if k in ids}; pairs=sorted(set(r1)&set(r2)); singles1=sorted(set(r1)-set(r2)); singles2=sorted(set(r2)-set(r1))
 for out,items in [(a.out_prefix+'_R1.fastq',[(k,r1[k]) for k in pairs+singles1]),(a.out_prefix+'_R2.fastq',[(k,r2[k]) for k in pairs+singles2])]: Path(out).write_text(''.join(v for _,v in items))
 Path(a.stats).write_text(json.dumps({'library':a.library,'requested':len(ids),'recovered_r1':len(r1),'recovered_r2':len(r2),'paired':len(pairs),'singleton_r1':len(singles1),'singleton_r2':len(singles2),'missing':len(ids-set(r1)-set(r2))},indent=2)+'\n')
if __name__=='__main__': main()
