#!/usr/bin/env python3
"""Compare a candidate indel policy with a baseline across frozen metric cells.

The comparator is deliberately conservative: a shared rule qualifies only if it
never loses a baseline TP, never adds a baseline FP, and improves at least one
cell. It emits a no-qualifying-rule outcome instead of selecting a strongest
caller or tuning thresholds per dataset. Inputs are existing som.py metrics; no
workflow, mapping, or caller process is launched.
"""
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def indel_row(path):
    doc = json.loads(path.read_text())
    tables = doc.get("metrics", [])
    if not tables or not tables[0].get("data"):
        raise ValueError(f"metrics table missing: {path}")
    values = {row["id"]: row.get("values", []) for row in tables[0]["data"]}
    try:
        idx = next(i for i, value in enumerate(values["type"])
                   if str(value).lower() == "indels")
        return {key: values[key][idx] for key in ("tp", "fp", "fn", "precision", "recall")}
    except (KeyError, StopIteration, IndexError) as exc:
        raise ValueError(f"indel row missing or malformed: {path}") from exc


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidate", action="append", required=True, metavar="CELL=METRICS",
                    help="Candidate som.py metrics JSON; repeat for each frozen cell")
    ap.add_argument("--baseline", action="append", required=True, metavar="CELL=METRICS",
                    help="Baseline som.py metrics JSON matching every candidate cell")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    def parse(items):
        result = {}
        for item in items:
            if "=" not in item:
                ap.error(f"expected CELL=METRICS, got {item!r}")
            cell, raw = item.split("=", 1)
            if not cell or cell in result:
                ap.error(f"duplicate/empty cell: {cell!r}")
            path = Path(raw).resolve()
            if not path.is_file():
                ap.error(f"metrics file not found: {path}")
            result[cell] = path
        return result

    candidates = parse(args.candidate)
    baselines = parse(args.baseline)
    if set(candidates) != set(baselines):
        ap.error(f"candidate/baseline cells differ: {sorted(set(candidates) ^ set(baselines))}")

    rows = []
    for cell in sorted(candidates):
        candidate = indel_row(candidates[cell])
        baseline = indel_row(baselines[cell])
        non_regression = (candidate["tp"] >= baseline["tp"] and
                          candidate["fp"] <= baseline["fp"])
        strict_gain = (candidate["tp"] > baseline["tp"] or
                       candidate["fp"] < baseline["fp"] or
                       (candidate.get("f1") is not None and baseline.get("f1") is not None and
                        candidate["f1"] > baseline["f1"]))
        rows.append({"cell": cell, "candidate": candidate, "baseline": baseline,
                     "non_regression": non_regression, "strict_gain": strict_gain,
                     "candidate_source": str(candidates[cell]),
                     "candidate_sha256": digest(candidates[cell]),
                     "baseline_source": str(baselines[cell]),
                     "baseline_sha256": digest(baselines[cell])})

    qualifies = bool(rows) and all(row["non_regression"] for row in rows) and any(row["strict_gain"] for row in rows)
    result = {
        "status": "pass" if qualifies else "no_qualifying_shared_rule",
        "decision": "select_candidate" if qualifies else "retain_existing_threshold_consensus",
        "selection_criteria": {"never_lose_tp": True, "never_add_fp": True, "strict_gain_required": True},
        "candidate_policy": "shared_candidate",
        "baseline_policy": "DeepSomatic_reference",
        "rows": rows,
        "independent_validation": "not_run",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "cells": len(rows), "decision": result["decision"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
