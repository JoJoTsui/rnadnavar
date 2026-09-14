#!/usr/bin/env python3
"""Validate that rerun label paths identify the final realignment-rescue artifact.

This is a read-only manifest/config check. It rejects intermediate ``*.rescued``
artifacts and first-round RNA paths when realignment rescue is required.
"""
import argparse
import csv
import json
from pathlib import Path

FINAL_SUFFIXES = (".filtered.vcf.gz", ".filtered.vcf.stripped.vcf.gz")

def validate_rows(rows, *, require_realign=True):
    results = []
    for row in rows:
        sample = row.get("sample_id") or row.get("sample") or "<unknown>"
        path = str(row.get("rescue_vcf_path") or row.get("training_label_vcf") or "")
        reasons = []
        if not path:
            reasons.append("missing rescue VCF path")
        if path.endswith(".rescued.vcf.gz") or "/rescued/" in path and ".filtered." not in path:
            reasons.append("intermediate rescue artifact")
        if path and not path.endswith(FINAL_SUFFIXES):
            reasons.append("not a filtered final VCF")
        if require_realign and path and "vcf_realignment" not in path and "RT_realign" not in path:
            reasons.append("realignment lineage is not identifiable")
        results.append({"sample_id": sample, "path": path,
                        "status": "pass" if not reasons else "error",
                        "reasons": reasons})
    return results

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--allow-first-round", action="store_true",
                    help="allow first-round rescue paths; still rejects intermediates")
    args=ap.parse_args()
    with args.manifest.open(newline="") as handle:
        rows=list(csv.DictReader(handle, delimiter="\t"))
    results=validate_rows(rows, require_realign=not args.allow_first_round)
    report={"status":"pass" if all(r["status"]=="pass" for r in results) else "blocked",
            "require_realign":not args.allow_first_round, "rows":results}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status":report["status"],"rows":len(results),
                      "failures":sum(r["status"]!="pass" for r in results)}))
    return 0 if report["status"]=="pass" else 1

if __name__ == "__main__":
    raise SystemExit(main())
