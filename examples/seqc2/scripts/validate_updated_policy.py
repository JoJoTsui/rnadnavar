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
    parser.add_argument("--read-evidence", type=Path, required=True)
    parser.add_argument("--ingress-report", type=Path, required=True)
    parser.add_argument("--reachability-report", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error(f"output directory already exists: {args.outdir}; choose a fresh namespace")
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

    read_evidence_path = args.read_evidence.resolve()
    if read_evidence_path.exists():
        read_evidence = json.loads(read_evidence_path.read_text())
        evidence["read_evidence_sha256"] = digest(read_evidence_path)
        sites = read_evidence.get("sites", [])
        expected_sites = {(row["dataset"], tuple(row[k] for k in ("chrom", "pos", "ref", "alt"))) for row in json.loads((repo / "examples/seqc2/comparison/rescue_fp_investigation_20260914/evidence.json").read_text())["sites"]}
        actual_sites = {(row.get("dataset"), tuple(row.get("allele", []))) for row in sites}
        streams = ("dna_tumor", "dna_normal", "rna_tumor", "rna_realign")
        stream_counts = {stream: {status: sum(1 for row in sites if row.get(stream, {}).get("status") == status) for status in ("observed", "inconclusive")} for stream in streams}
        scored = [row for row in sites if row.get("role") == "scored_rescue"]
        complete = len(scored) == len(expected_sites) and {(row.get("dataset"), tuple(row.get("allele", []))) for row in scored} == expected_sites and len(sites) == 26 and all(sum(counts.values()) == len(sites) for counts in stream_counts.values())
        source_stats = read_evidence.get("source_stats", {})
        expected_alignments = {str((bundle / dataset / "alignments" / filename).resolve()) for dataset in ("wes_ll", "wgs_il") for filename in ("dna_tumor.bam", "dna_normal.bam", "rna_tumor.cram", "rna_realign.cram")}
        sources_ok = expected_alignments <= set(source_stats) and str(Path("/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta")) in source_stats
        checks.append(check("read-evidence:available", complete and sources_ok, f"stream counts {stream_counts}; roles={sorted(set(row.get('role') for row in sites))}"))
    else:
        checks.append(check("read-evidence:available", False, "bounded DNA pileup evidence not run", severity="inconclusive"))

    ingress_path = args.ingress_report.resolve()
    if ingress_path.exists():
        ingress = json.loads(ingress_path.read_text())
        evidence["ingress_report_sha256"] = digest(ingress_path)
        checks.append(check("ingress-contract:available", ingress.get("status") == "pass", ingress.get("status", "missing")))
    else:
        checks.append(check("ingress-contract:available", False, "ingress contract report not run", severity="inconclusive"))

    reachability_path = args.reachability_report.resolve()
    if reachability_path.exists():
        reachability = json.loads(reachability_path.read_text())
        evidence["reachability_report_sha256"] = digest(reachability_path)
        checks.append(check("process-reachability:available", reachability.get("status") == "pass", reachability.get("status", "missing")))
    else:
        checks.append(check("process-reachability:available", False, "process reachability report not run", severity="inconclusive"))

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
                            ((label == "rerun" and "tools: consensus,rescue,filtering,vep" in config) or (label == "cohort-rerun" and "tools: consensus,rescue,filtering" in config and "tools: consensus,rescue,filtering,vep" not in config)),
                            "consensus/rescue-only tool set present"))
        checks.append(check(f"{label}:no-caller-tools",
                            not re.search(r"tools:.*(?:mutect2|strelka|deepsomatic)", config),
                            "no variant caller is requested"))
        checks.append(check(f"{label}:consensus-step",
                            bool(re.search(r"step:\s+consensus\b", config)),
                            "step is consensus"))

        checks.append(check(f"{label}:checksum-output", "checksum_dir:" in config, "rerun has a dedicated checksum output namespace"))

    for label, config in (("rerun", rerun), ("cohort-rerun", cohort)):
        values = {}
        for key in ("outdir_base", "csv_dir", "state_file", "checksum_dir"):
            match = re.search(rf"^{key}:\s*(\S+)", config, re.MULTILINE)
            values[key] = match.group(1) if match else ""
        checks.append(check(f"{label}:isolated-namespaces", len(set(values.values())) == len(values) and all(values.values()), str(values)))
        checks.append(check(f"{label}:new-output-namespace", "output_reconsensus" in values["outdir_base"], values["outdir_base"]))

    mode_configs = [repo / "examples/seqc2/seqc2.shared.config", repo / "examples/seqc2/seqc2.hybrid.shared.config", repo / "examples/seqc2/hybrid/seqc2.hybrid.wgs.config"]
    for path in mode_configs:
        text = path.read_text()
        checks.append(check(f"input-mode:{path.name}", "consensus" in text and "rescue" in text, "consensus/rescue tools are declared"))

    tests = [
        [sys.executable, "-m", "pytest", "tests/test_rescue_nomination_experiment.py",
         "tests/test_rescue_gate_domain_audit.py", "tests/test_historical_native_gate_replay.py", "tests/test_rescue_read_evidence.py", "-q"],
    ]
    test_logs = []
    for command in tests:
        result = subprocess.run(command, cwd=repo, text=True, capture_output=True)
        test_logs.append({"command": command, "returncode": result.returncode,
                          "stdout": result.stdout, "stderr": result.stderr})
    checks.append(check("focused-tests", all(item["returncode"] == 0 for item in test_logs),
                        "19 focused policy tests pass"))

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
