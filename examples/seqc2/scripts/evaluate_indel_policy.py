#!/usr/bin/env python3
"""Evaluate indel policy candidates from score-label JSON artifacts.

The script compares supplied candidates against truth and reports an evidence
summary; it does not select or enable a production policy.
"""
import argparse, json
from pathlib import Path

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--scores', type=Path, nargs='+', required=True,
                    help='score_label_artifact JSON files')
    ap.add_argument('--out', type=Path, required=True)
    a=ap.parse_args(); rows=[]
    for path in a.scores:
        doc=json.loads(path.read_text())
        row=next((r for r in doc.get('rows',[]) if r.get('variant_type','').lower()=='indel'), None)
        if row is None: raise SystemExit(f'indel row missing: {path}')
        rows.append({"source":str(path.resolve()),"source_sha256":__import__('hashlib').sha256(path.read_bytes()).hexdigest(),**row})
    result={"status":"descriptive_not_policy_selection",
            "decision":"requires_independent_validation",
            "criteria":["precision","recall","f1","fp","fn"],"rows":rows}
    a.out.parent.mkdir(parents=True,exist_ok=True); a.out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"rows":len(rows)}))
if __name__=='__main__': main()
