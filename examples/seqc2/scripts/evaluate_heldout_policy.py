#!/usr/bin/env python3
"""Evaluate a frozen policy on held-out slices without retuning."""
import argparse
import json
import math
from pathlib import Path


def metric(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{label} must be a finite number between 0 and 1")
    return float(value)


def validate_bound(frozen, held):
    required = frozen.get("required_slices")
    if not isinstance(required, list) or not required:
        raise ValueError("frozen policy has no required slices")
    expected_partition = frozen.get("partitions", {}).get("holdout")
    if held.get("partition") != expected_partition:
        raise ValueError("held-out partition does not match frozen holdout identity")
    if held.get("required_slices") != required:
        raise ValueError("held-out required slices do not match frozen policy")
    baseline = frozen.get("baseline")
    if not isinstance(baseline, dict) or set(baseline) != set(required):
        raise ValueError("frozen policy has no complete baseline binding")
    rows = held.get("rows")
    if not isinstance(rows, list):
        raise ValueError("held-out metrics must contain rows")
    by = {}
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get("pair"), str) or not isinstance(row.get("variant_type"), str):
            raise ValueError("held-out rows must declare pair and variant_type")
        key = f"{row['pair']}:{row['variant_type']}"
        if key in by:
            raise ValueError(f"held-out metrics contain duplicate slice {key}")
        by[key] = row
        metric(row.get("precision"), f"held-out {key} precision")
        metric(row.get("recall"), f"held-out {key} recall")
    extra = [key for key in by if key not in required]
    if extra:
        raise ValueError("held-out metrics contain undeclared extra slices: " + ", ".join(sorted(extra)))
    missing = [key for key in required if key not in by]
    if missing:
        raise ValueError("held-out metrics are missing required slices: " + ", ".join(missing))
    return required, baseline, by


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    frozen = json.loads(Path(args.frozen).read_text())
    held = json.loads(Path(args.metrics).read_text())
    if frozen.get("status") != "qualified" or not frozen.get("frozen_policy"):
        result = {"status": "incomplete_no_frozen_policy", "rules_source": "frozen_development_policy", "retuning_allowed": False, "rows": [], "limitations": held.get("limitations", [])}
        Path(args.out).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        return
    try:
        required, baseline, by = validate_bound(frozen, held)
        min_delta = frozen.get("min_delta", 0.0)
        metric(min_delta, "frozen min_delta")
        rows = []
        for key in required:
            row = dict(by[key])
            base = baseline[key]
            metric(base.get("precision"), f"baseline {key} precision")
            metric(base.get("recall"), f"baseline {key} recall")
            row["precision_gate"] = row["precision"] >= base["precision"] + min_delta
            row["recall_gate"] = row["recall"] >= base["recall"] + min_delta
            row["accepted"] = row["precision_gate"] and row["recall_gate"]
            rows.append(row)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.error(str(error))
    result = {"status": "complete", "rules_source": "frozen_development_policy", "retuning_allowed": False, "required_slices": required, "rows": rows, "limitations": held.get("limitations", [])}
    Path(args.out).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
