#!/usr/bin/env python3
"""Select and freeze a label policy from declared development metrics."""
import argparse
import json
from pathlib import Path
from policy_metrics import f1, metric


def _schema(data):
    if not isinstance(data, dict):
        raise ValueError("metrics must be a JSON object")
    required = data.get("required_slices")
    if not isinstance(required, list) or not required or any(not isinstance(x, str) or ":" not in x for x in required):
        raise ValueError("required_slices must be a non-empty list of pair:variant_type keys")
    if len(required) != len(set(required)):
        raise ValueError("required_slices contains duplicates")
    partitions = data.get("partitions")
    if not isinstance(partitions, dict) or not partitions.get("development") or not partitions.get("holdout"):
        raise ValueError("partitions must declare development and holdout identities")
    baseline = data.get("baseline")
    if not isinstance(baseline, dict) or not baseline:
        raise ValueError("baseline must contain required comparison slices")
    if set(baseline) != set(required):
        raise ValueError("baseline keys must exactly match required_slices")
    for key, row in baseline.items():
        if not isinstance(row, dict):
            raise ValueError(f"baseline slice {key} must be an object")
        metric(row.get("precision"), f"baseline {key} precision")
        metric(row.get("recall"), f"baseline {key} recall")
        row["f1"] = f1(row, f"baseline {key}")
    candidates = data.get("candidates")
    if not isinstance(candidates, list):
        raise ValueError("candidates must be a list")
    return required, baseline, candidates


def _candidate_rows(candidate, required):
    if not isinstance(candidate, dict) or not isinstance(candidate.get("policy"), dict):
        raise ValueError("each candidate must contain a policy object")
    rows = candidate.get("rows")
    if not isinstance(rows, list):
        raise ValueError("each candidate must contain a rows list")
    by = {}
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get("pair"), str) or not isinstance(row.get("variant_type"), str):
            raise ValueError("candidate rows must declare pair and variant_type")
        key = f"{row['pair']}:{row['variant_type']}"
        if key in by:
            raise ValueError(f"candidate contains duplicate slice {key}")
        by[key] = row
        metric(row.get("precision"), f"candidate {key} precision")
        metric(row.get("recall"), f"candidate {key} recall")
        row["f1"] = f1(row, f"candidate {key}")
    extra = [key for key in by if key not in required]
    if extra:
        raise ValueError("candidate contains undeclared extra slices: " + ", ".join(sorted(extra)))
    missing = [key for key in required if key not in by]
    if missing:
        raise ValueError("candidate is missing required slices: " + ", ".join(missing))
    return by


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-delta", type=float, default=0.0)
    parser.add_argument("--min-f1-delta", type=float, default=0.0,
                        help="Minimum candidate F1 delta versus the declared baseline")
    parser.add_argument("--min-precision-delta", type=float, default=0.0,
                        help="Minimum candidate precision delta versus the declared baseline")
    args = parser.parse_args()
    try:
        data = json.loads(Path(args.metrics).read_text())
        required, baseline, candidates = _schema(data)
        min_delta = metric(args.min_delta, "min-delta")
        min_f1_delta = metric(args.min_f1_delta, "min-f1-delta")
        min_precision_delta = metric(args.min_precision_delta, "min-precision-delta")
        decisions = []
        for candidate in candidates:
            by = _candidate_rows(candidate, required)
            failures = []
            for key in required:
                base = baseline[key]
                row = by[key]
                if row["precision"] < base["precision"] + max(min_delta, min_precision_delta):
                    failures.append({"key": key, "reason": "precision_gate", "candidate": row["precision"], "baseline": base["precision"]})
                if row["recall"] < base["recall"] + min_delta:
                    failures.append({"key": key, "reason": "recall_gate", "candidate": row["recall"], "baseline": base["recall"]})
                if row["f1"] < base["f1"] + min_f1_delta:
                    failures.append({"key": key, "reason": "f1_gate", "candidate": row["f1"], "baseline": base["f1"]})
            decisions.append({"policy": candidate["policy"], "qualifies": not failures, "failures": failures})
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.error(str(error))
    qualified = sorted((decision for decision in decisions if decision["qualifies"]), key=lambda decision: json.dumps(decision["policy"], sort_keys=True, separators=(",", ":")))
    result = {
        "status": "qualified" if qualified else "no_qualifying_policy",
        "development_only": True,
        "required_slices": required,
        "partitions": data["partitions"],
        "baseline": baseline,
        "min_delta": min_delta,
        "min_f1_delta": min_f1_delta,
        "min_precision_delta": min_precision_delta,
        "decisions": decisions,
        "frozen_policy": qualified[0]["policy"] if qualified else None,
    }
    Path(args.out).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
