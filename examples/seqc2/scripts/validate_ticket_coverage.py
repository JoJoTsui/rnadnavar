#!/usr/bin/env python3
"""Summarize implementation coverage for updated-policy tickets 01-07.

This consumes a completed validation report only; it never launches workflow
processes or modifies source, cache, or pipeline outputs.
"""
import argparse
import json
from pathlib import Path


def _status(report, names):
    checks = {row["name"]: row["status"] for row in report.get("checks", [])}
    values = [checks.get(name, "not_run") for name in names]
    if any(value in {"error", "blocked"} for value in values):
        return "fail"
    if any(value in {"inconclusive", "not_run"} for value in values):
        return "inconclusive"
    return "pass"


def build_matrix(report):
    return [
        {"ticket": "01", "status": _status(report, [
            "benchmark:wes_ll/ukb", "benchmark:wes_ll/medexome",
            "benchmark:wgs_il/ukb", "benchmark:wgs_il/medexome",
        ]), "evidence": "benchmark cells and source hashes"},
        {"ticket": "02", "status": _status(report, ["read-evidence:available"]),
         "evidence": "four-stream read/depth/strand metrics"},
        {"ticket": "03", "status": _status(report, [
            "default:native-flag", "policy:experimental-not-promoted",
            "policy:no-production-nomination-wiring",
        ]), "evidence": "default wiring and policy-difference checks"},
        {"ticket": "04", "status": _status(report, [
            "read-evidence:available", "process-reachability:available",
        ]), "evidence": "RNA and realignment route evidence; annotation parity remains bounded"},
        {"ticket": "05", "status": _status(report, ["focused-tests"]),
         "evidence": "adversarial and label-contract regression tests"},
        {"ticket": "06", "status": _status(report, [
            "ingress-contract:available", "process-reachability:available",
            "rerun:checksum-output", "cohort-rerun:checksum-output",
        ]), "evidence": "ingress, route, namespace, and cache-isolation checks"},
        {"ticket": "07", "status": "pass_with_known_gates" if report.get("status") == "pass_with_known_gates" else report.get("status", "not_run"),
         "evidence": "final verdict and production-fix backlog"},
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = json.loads(args.report.read_text())
    matrix = {"source_report": str(args.report.resolve()), "tickets": build_matrix(report)}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(matrix, indent=2) + "\n")
    print(json.dumps(matrix, indent=2))
    return 0 if all(row["status"] in {"pass", "pass_with_known_gates"} for row in matrix["tickets"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
