#!/usr/bin/env python3
"""Standalone, read-only-source parity check against frozen SEQC2 assays.

Runs only the consensus Python driver on existing normalized caller VCFs.
Does not launch Nextflow, mapping, calling, rescue, or HG008 evaluation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pysam


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def alleles(path, somatic=False):
    with pysam.VariantFile(str(path)) as handle:
        return {(r.contig, r.pos, r.ref, ",".join(r.alts or [])) for r in handle
                if not somatic or "Somatic" in r.filter}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir = args.outdir.resolve()
    args.outdir.mkdir(parents=True, exist_ok=False)
    repo = Path(__file__).resolve().parents[3]
    report = {"status": "running", "sources": {}, "cells": {}, "code": {}}
    for p in [repo / "bin/run_consensus_vcf.py", *sorted((repo / "bin/vcf_utils").glob("*.py"))]:
        report["code"][str(p)] = digest(p)
    for dataset in ("wes", "wgs"):
        source = args.root / ("wes_indel_attribution" if dataset == "wes" else "wgs_attribution")
        for region in ("ukb", "medexome"):
            dest = args.outdir / dataset / region
            inputs = dest / "inputs"
            inputs.mkdir(parents=True)
            for caller in ("deepsomatic", "mutect2", "strelka"):
                p = (source / region / f"{caller}.normalized.vcf.gz").resolve()
                report["sources"][str(p)] = digest(p)
                staged = inputs / f"sample.{caller}.vcf.gz"
                staged.symlink_to(p)
                # Make a local index, never modify a source or its index.
                subprocess.run(["bcftools", "index", "-t", "-o", str(staged) + ".tbi", str(staged)], check=True)
            expected = (args.root / f"{dataset}_refined_indel" / region / "combined/query.vcf.gz").resolve()
            report["sources"][str(expected)] = digest(expected)
            cmd = [sys.executable, str(repo / "bin/run_consensus_vcf.py"),
                   "--input_dir", str(inputs), "--expected_callers", "deepsomatic,mutect2,strelka",
                   "--out_prefix", str(dest / "refined"), "--experimental-refined-native"]
            with (dest / "consensus.log").open("w") as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
            outputs = list(dest.glob("refined*.vcf.gz"))
            if len(outputs) != 1:
                raise RuntimeError(f"Ambiguous consensus output: {outputs}")
            actual, frozen = alleles(outputs[0], True), alleles(expected)
            cell = {"command": cmd, "output": str(outputs[0]), "expected": str(expected),
                    "actual_count": len(actual), "expected_count": len(frozen),
                    "extra": sorted(actual - frozen), "missing": sorted(frozen - actual)}
            report["cells"][dataset + "/" + region] = cell
            (args.outdir / "parity.json").write_text(json.dumps(report, indent=2) + "\n")
            print(dataset, region, "extra", len(cell["extra"]), "missing", len(cell["missing"]), flush=True)
    report["sources_unchanged"] = all(digest(Path(p)) == sha for p, sha in report["sources"].items())
    report["code_unchanged"] = all(digest(Path(p)) == sha for p, sha in report["code"].items())
    report["status"] = ("parity_passed" if report["sources_unchanged"] and report["code_unchanged"]
                        and all(not c["extra"] and not c["missing"] for c in report["cells"].values())
                        else "parity_failed")
    (args.outdir / "parity.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["status"] == "parity_passed" else 1


if __name__ == "__main__":
    sys.exit(main())
