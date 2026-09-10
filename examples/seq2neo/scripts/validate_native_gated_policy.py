#!/usr/bin/env python3
"""Evaluate native-consensus + gated-rescue release criteria.

The script consumes precomputed benchmark slices and does not invoke Nextflow,
hap.py, mapping, or variant callers.  Input JSON is a list of slices with
candidate and deepsomatic metrics (precision, f1, tp, fp, fn) plus optional
valid_input_failures.  It writes a machine-readable gate decision.
"""
import argparse
import json
import math
from pathlib import Path

def _metric(obj, name):
    value = obj.get(name)
    if value is None or not math.isfinite(float(value)):
        raise ValueError(f"missing/non-finite metric: {name}")
    return float(value)


def evaluate(doc):
    slices = doc.get("slices")
    if not isinstance(slices, list) or not slices:
        raise ValueError("comparison JSON must contain a non-empty 'slices' list")
    results = []
    reasons = []
    for index, row in enumerate(slices):
        label = f"{row.get('cohort', 'slice')}:{row.get('variant_type', 'records')}"
        candidate = row.get("candidate") or {}
        baseline = row.get("deepsomatic") or {}
        try:
            cp = _metric(candidate, "precision")
            cf = _metric(candidate, "f1")
            bp = _metric(baseline, "precision")
            bf = _metric(baseline, "f1")
            cfp = int(candidate["fp"])
            bfp = int(baseline["fp"])
            ctp = int(candidate["tp"])
            btp = int(baseline["tp"])
            cfn = int(candidate["fn"])
            bfn = int(baseline["fn"])
        except (KeyError, TypeError, ValueError) as exc:
            reasons.append(f"{label}: invalid metrics ({exc})")
            results.append({"slice": label, "accepted": False, "reason": str(exc)})
            continue
        failures = int(row.get("valid_input_failures", 0))
        checks = {
            "f1_not_lower": cf + 1e-12 >= bf,
            "precision_not_lower": cp + 1e-12 >= bp,
            "tp_not_lower": ctp >= btp,
            "fp_not_higher": cfp <= bfp,
            "fn_not_higher": cfn <= bfn,
            "no_valid_input_failures": failures == 0,
        }
        accepted = all(checks.values())
        if not accepted:
            failed = [name for name, ok in checks.items() if not ok]
            reasons.append(f"{label}: failed {', '.join(failed)}")
        results.append({"slice": label, "accepted": accepted, "checks": checks,
                        "candidate": candidate, "deepsomatic": baseline,
                        "valid_input_failures": failures})
    return {"status": "passed" if not reasons else "failed",
            "slices": results, "reasons": reasons}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-json", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        result = evaluate(json.loads(args.comparison_json.read_text()))
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        parser.error(str(exc))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "slices": len(result["slices"]), "reasons": result["reasons"]}))
    return 0 if result["status"] == "passed" else 1

if __name__ == "__main__":
    raise SystemExit(main())
