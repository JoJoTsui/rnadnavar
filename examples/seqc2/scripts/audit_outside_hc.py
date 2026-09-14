#!/usr/bin/env python3
"""Publish outside-HC rescue additions as explicitly unassessed.

No truth comparison is made outside the declared high-confidence region.
"""
import argparse, hashlib, json
from pathlib import Path

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--domain-audit', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a=ap.parse_args(); doc=json.loads(a.domain_audit.read_text())
    rows=[]
    for dataset, entries in doc.get('summary', {}).items():
        for row in entries:
            if row.get('domain') != 'outside_hc_unassessed': continue
            rows.append({"dataset":dataset, **row, "truth_status":"unassessed"})
    result={"status":"pass", "truth_policy":"outside_hc_is_unassessed",
            "source_sha256":hashlib.sha256(a.domain_audit.read_bytes()).hexdigest(),
            "rows":rows, "total":sum(int(r.get('count',0)) for r in rows)}
    a.out.parent.mkdir(parents=True, exist_ok=True); a.out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"total":result["total"]}))
if __name__=='__main__': main()
