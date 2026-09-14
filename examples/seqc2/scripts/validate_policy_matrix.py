#!/usr/bin/env python3
"""Enforce shared-policy benchmark gates across all matrix cells.

The candidate must beat every declared DNA comparator on F1 and match or beat
DeepSomatic precision for each cell and variant class. CSVs are produced by
aggregate_benchmark.py and are never modified.
"""
import argparse, csv, json
from pathlib import Path

def load(path):
    with open(path, newline="") as fh:
        return {(r["query"], r["variant_type"]): r for r in csv.DictReader(fh)}

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True, type=Path, help="Root containing dataset/domain subdirectories")
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--comparators", nargs="+", default=["mutect2", "strelka", "deepsomatic"])
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(); decisions = []; failures = []
    cells = sorted(p for p in a.matrix.glob("*/*/*.benchmark_comparison.csv"))
    if len(cells) != 4:
        raise SystemExit(f"expected four dataset/domain CSVs, found {len(cells)}")
    for path in cells:
        rows = load(path); cell = str(path.parent.relative_to(a.matrix))
        for kind in ("snp", "indel", "records"):
            key = (a.candidate, kind)
            if key not in rows: failures.append({"cell": cell, "type": kind, "reason": "candidate_missing"}); continue
            cand = rows[key]; p = float(cand["precision"]); f1 = float(cand["f1"])
            ds = rows.get(("deepsomatic", kind))
            if ds is None or p < float(ds["precision"]): failures.append({"cell": cell, "type": kind, "reason": "deep_somatic_precision"})
            for comparator in a.comparators:
                base = rows.get((comparator, kind))
                if base is None or f1 <= float(base["f1"]): failures.append({"cell": cell, "type": kind, "reason": f"f1_not_above_{comparator}"})
            decisions.append({"cell": cell, "variant_type": kind, "candidate": cand})
    result = {"status": "qualified" if not failures else "no_qualifying_policy", "candidate": a.candidate, "comparators": a.comparators, "cells": [str(p.parent.relative_to(a.matrix)) for p in cells], "decisions": decisions, "failures": failures}
    a.out.parent.mkdir(parents=True, exist_ok=True); a.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True)); raise SystemExit(0 if not failures else 2)
if __name__ == "__main__": main()
