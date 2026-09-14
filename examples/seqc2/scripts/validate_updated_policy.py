#!/usr/bin/env python3
"""Run the bounded, read-only validation for the frozen rescue experiment.

This validator consumes completed benchmark/evidence artifacts and inspects
effective configuration wiring. It never launches Nextflow, mapping, calling,
or annotation; no source VCF, cache, or workflow output is modified.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys


EXPECTED = {
    "wes_ll/ukb": (1062, 38, 1238),
    "wes_ll/medexome": (570, 19, 259),
    "wgs_il/ukb": (2169, 19, 131),
    "wgs_il/medexome": (702, 6, 127),
}


def check(name, passed, detail, *, severity="error"):
    return {"name": name, "status": "pass" if passed else severity,
            "detail": detail}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def metric_tuple(metrics):
    rows = metrics["records"]
    return tuple(rows[key] for key in ("tp", "fp", "fn"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[3]
    bundle = repo / "examples/seqc2/verified/20260914"
    gate_root = repo / "examples/seqc2/comparison/rescue_fp_investigation_20260914/gate_tests"
    domain_path = repo / "examples/seqc2/comparison/rescue_fp_investigation_20260914/domain_audit/domain_audit.json"
    checks = []
    evidence = {}

    evaluation_path = gate_root / "evaluation.json"
    evaluation = json.loads(evaluation_path.read_text())
    evidence["evaluation_sha256"] = digest(evaluation_path)
    metrics = evaluation.get("metrics", {})
    for cell, expected in EXPECTED.items():
        key = f"{cell.split('/')[0]}/nomination_biological/{cell.split('/')[1]}"
        actual = metric_tuple(metrics[key])
        checks.append(check(f"benchmark:{cell}", actual == expected,
                            f"expected {expected}, observed {actual}"))

    domain = json.loads(domain_path.read_text())
    evidence["domain_audit_sha256"] = digest(domain_path)
    for dataset, rows in domain["summary"].items():
        total = sum(row["count"] for row in rows)
        checks.append(check(f"domain-accounting:{dataset}",
                            total == {"wes_ll": 303, "wgs_il": 610}[dataset],
                            f"accounted for {total} rescue additions"))
        truth_excluded = [row for row in rows
                          if row["truth_status"] == "truth_present" and row["decision"] == "excluded"]
        checks.append(check(f"domain-retention:{dataset}", not truth_excluded,
                            f"truth-present exclusions: {truth_excluded}"))

    nextflow = (repo / "nextflow.config").read_text()
    checks.extend([
        check("default:native-flag", "native_evidence_snv             = true" in nextflow,
              "root native_evidence_snv default is true"),
        check("default:rescue-promotion", "rescue_promotion_enabled        = true" in nextflow,
              "root rescue promotion default is true"),
        check("default:rescue-floors", "rescue_min_dna_callers          = 1" in nextflow and
              "rescue_min_rna_callers          = 2" in nextflow,
              "root rescue floors are DNA=1 and RNA=2"),
    ])
    experimental = (repo / "examples/seqc2/scripts/test_rescue_nomination_gate.py").read_text()
    checks.append(check("policy:experimental-not-promoted",
                        "exploratory_not_promoted" in experimental,
                        "nomination+biological gate remains an isolated experiment"))
    tracked = subprocess.check_output(["git", "ls-files", "bin", "modules", "subworkflows", "conf"], cwd=repo, text=True).splitlines()
    tracked_code = [repo / path for path in tracked if Path(path).suffix in {".py", ".nf", ".config"}]
    unexpected = [str(path.relative_to(repo)) for path in tracked_code
                  if "nomination_biological" in path.read_text(errors="replace")]
    checks.append(check("policy:no-production-nomination-wiring", not unexpected,
                        f"production references: {unexpected}"))

    rerun = (repo / "examples/seq2neo/config/rerun.yaml").read_text()
    cohort = (repo / "examples/seq2neo/config/rerun_cohort.yaml").read_text()
    for label, config in (("rerun", rerun), ("cohort-rerun", cohort)):
        checks.append(check(f"{label}:consensus-rescue-tools",
                            "tools: consensus,rescue,filtering,vep" in config,
                            "consensus/rescue-only tool set present"))
        checks.append(check(f"{label}:no-caller-tools",
                            not re.search(r"tools:.*(?:mutect2|strelka|deepsomatic)", config),
                            "no variant caller is requested"))
        checks.append(check(f"{label}:consensus-step",
                            bool(re.search(r"step:\s+consensus\b", config)),
                            "step is consensus"))

        checks.append(check(f"{label}:checksum-output", "checksum_dir:" in config, "rerun has a dedicated checksum output namespace"))

    tests = [
        [sys.executable, "-m", "pytest", "tests/test_rescue_nomination_experiment.py",
         "tests/test_rescue_gate_domain_audit.py", "tests/test_historical_native_gate_replay.py", "-q"],
    ]
    test_logs = []
    for command in tests:
        result = subprocess.run(command, cwd=repo, text=True, capture_output=True)
        test_logs.append({"command": command, "returncode": result.returncode,
                          "stdout": result.stdout, "stderr": result.stderr})
    checks.append(check("focused-tests", all(item["returncode"] == 0 for item in test_logs),
                        "12+ focused policy tests pass"))

    report = {
        "status": "pass_with_known_gates" if all(item["status"] == "pass" for item in checks)
                  else "blocked_or_inconclusive",
        "policy": "frozen_nomination_biological_gate",
        "independent_validation": "not_run",
        "production_defaults_changed": False,
        "mapping_or_calling_run": False,
        "checks": checks,
        "evidence": evidence,
        "test_logs": test_logs,
        "known_gates": [
            "HG008 independent validation is not complete",
            "production workflow is not equivalent to the frozen experimental gate",
            "indel policy remains provisional",
            "outside-HC additions are unassessed",
        ],
    }
    args.outdir.mkdir(parents=True, exist_ok=False)
    (args.outdir / "validation_report.json").write_text(json.dumps(report, indent=2) + "\n")
    failures = [item for item in checks if item["status"] != "pass"]
    print(json.dumps({"status": report["status"], "checks": len(checks),
                      "failures": failures}, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
