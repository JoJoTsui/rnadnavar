#!/usr/bin/env python3
"""Validate the explicit, opt-in native/gated policy profile."""
import argparse, json
from pathlib import Path

EXPECTED={"policy_profile":"native_gated_experimental","native_evidence_snv":True,"rescue_promotion_enabled":True,"rescue_min_dna_callers":1,"rescue_min_rna_callers":2,"rescue_veto":"dna","indel_policy":"threshold_consensus"}
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--profile',type=Path,required=True); ap.add_argument('--out',type=Path,required=True); a=ap.parse_args()
 values={}
 for line in a.profile.read_text().splitlines():
  line=line.split('#',1)[0].strip()
  if not line or ':' not in line: continue
  k,v=(x.strip() for x in line.split(':',1)); values[k]=({'true':True,'false':False}.get(v.lower(), int(v) if v.isdigit() else v.strip('"\'')))
 failures=[k for k,v in EXPECTED.items() if values.get(k)!=v]
 result={"status":"pass" if not failures else "blocked","default_loaded":False,"failures":failures,"values":values}
 a.out.parent.mkdir(parents=True,exist_ok=True); a.out.write_text(json.dumps(result,indent=2)+"\n"); print(json.dumps(result)); return 0 if not failures else 1
if __name__=='__main__': raise SystemExit(main())
