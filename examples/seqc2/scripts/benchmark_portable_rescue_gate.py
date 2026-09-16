#!/usr/bin/env python3
"""Benchmark caller-reconstructed rescue atop actual refined driver outputs.

Uses only SEQC2 data. Allele-only PASS queries are benchmark copies, never
training labels. Candidate decisions are completed before truth is read.
"""
import argparse
import json
from pathlib import Path
import subprocess

from aggregate_benchmark import parse_metrics_json
from audit_current_native_policy import digest
from replay_historical_native_gate import write_query
from validate_refined_native_integration import alleles


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--audit", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    audit = json.loads(args.audit.read_text())
    if audit["status"] != "complete_not_promoted" or not audit["sources_unchanged"]:
        raise ValueError("Require completed integrity-checked audit")
    args.outdir.mkdir(parents=True, exist_ok=False)
    report = {"status": "running", "sources": {str(args.audit.resolve()): digest(args.audit)},
              "script_sha256": digest(Path(__file__)), "cells": {}}
    for short, dataset in (("wes", "wes_ll"), ("wgs", "wgs_il")):
        candidate = audit["datasets"][dataset]
        if candidate["eligible_extra"] or candidate["eligible_missing"]:
            raise ValueError("Eligible rescue has unresolved historical parity differences")
        accepted = {tuple(k) for k in candidate["eligible_all_accepted"]}
        frozen_path = args.root / f"{short}_refined_gate_replay/evaluation.json"
        frozen = json.loads(frozen_path.read_text())
        report["sources"][str(frozen_path.resolve())] = digest(frozen_path)
        for region in ("ukb", "medexome"):
            source = args.root / "driver_parity_v1" / short / region / "refined.vcf.gz"
            report["sources"][str(source.resolve())] = digest(source)
            native = alleles(source, True)
            dest = args.outdir / short / region
            dest.mkdir(parents=True)
            query = write_query(dest / "query.vcf", native | accepted)
            cmd = list(frozen["cells"][region]["command"])
            previous = alleles(Path(cmd[6]))
            report["sources"][str(Path(cmd[6]).resolve())] = digest(Path(cmd[6]))
            cmd[6], cmd[-1] = query, str(dest / "benchmark")
            with (dest / "benchmark.log").open("w") as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
            metrics = parse_metrics_json(dest / "benchmark.metrics.json")
            report["cells"][short + "/" + region] = {
                "command": cmd, "metrics": metrics, "previous_metrics": frozen["cells"][region]["metrics"],
                "extra_vs_replay": sorted((native | accepted) - previous),
                "missing_vs_replay": sorted(previous - (native | accepted)),
                "native_records": len(native), "accepted_rescue_candidates": len(accepted)}
            (args.outdir / "evaluation.json").write_text(json.dumps(report, indent=2) + "\n")
            print(short, region, metrics["records"], flush=True)
    report["sources_unchanged"] = all(digest(Path(p)) == sha for p, sha in report["sources"].items())
    report["status"] = "complete_not_promoted" if report["sources_unchanged"] else "integrity_failure"
    (args.outdir / "evaluation.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
