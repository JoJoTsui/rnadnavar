#!/usr/bin/env python3
"""Compare two PASS/Somatic VCFs and report paired allele transitions.

The parser is intentionally stdlib-only and uses exact normalized VCF keys
(CHROM, POS, REF, ALT). It reports added/removed records and truth status when
an optional truth VCF is supplied. Inputs are never modified.
"""
import argparse
import gzip
import json
from pathlib import Path

def rows(path, include_all=False):
    opener = gzip.open if str(path).endswith(".gz") else open
    out = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            if not include_all and cols[6] not in {"PASS", "Somatic"}:
                continue
            for alt in cols[4].split(","):
                out[(cols[0], int(cols[1]), cols[3], alt)] = {"filter": cols[6], "info": cols[7]}
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--before", required=True, type=Path)
    ap.add_argument("--after", required=True, type=Path)
    ap.add_argument("--truth", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    before, after = rows(a.before), rows(a.after)
    truth = rows(a.truth, include_all=True) if a.truth else {}
    added = set(after) - set(before); removed = set(before) - set(after)
    result = {"before": str(a.before), "after": str(a.after),
              "counts": {"before": len(before), "after": len(after),
                         "added": len(added), "removed": len(removed)},
              "transitions": {"added_true": sum(k in truth for k in added),
                              "added_false": sum(k not in truth for k in added),
                              "removed_true": sum(k in truth for k in removed),
                              "removed_false": sum(k not in truth for k in removed)}}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
if __name__ == "__main__": main()
