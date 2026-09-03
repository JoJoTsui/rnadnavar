#!/usr/bin/env python3
"""Validate a bounded seq2neo execution trace and final VCF contracts."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path


REQUIRED_PROCESS_GROUPS = {
    "FASTQ trimming": ("FASTP",),
    "DNA alignment": ("BWAMEM", "BWA_MEM", "FASTQ_ALIGN"),
    "RNA alignment": ("STAR_ALIGN", "FASTQ_ALIGN_STAR"),
    "Mutect2": ("MUTECT2",),
    "Strelka2": ("STRELKA",),
    "DeepSomatic": ("DEEPSOMATIC",),
    "VCF normalization": ("BCFTOOLS_NORM", "VT_DECOMPOSE", "VCF_NORMALIZE"),
    "consensus": ("VCF_CONSENSUS",),
    "rescue": ("VCF_RESCUE",),
    "RNA realignment": ("HISAT2", "RNA_REALIGN"),
    "RNA editing": ("RNA_EDITING_ANNOTATION",),
    "COSMIC/gnomAD": ("COSMIC_GNOMAD_ANNOTATION",),
    "VEP": ("ENSEMBLVEP_VEP",),
}

FORBIDDEN_FASTQ_PROCESSES = (
    "CHECK_INPUT_DICTIONARY",
    "CHECK_EXTERNAL_CRAM_DICTIONARY",
    "PICARD_REORDERSAM",
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trace", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def read_trace(path: Path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    names = [" ".join((r.get("process", ""), r.get("name", ""))) for r in rows]
    failed = [r for r in rows if r.get("status", "").upper() == "FAILED"]
    return rows, names, failed


def consensus_vcfs(outdir: Path):
    return sorted(outdir.glob("consensus/**/*.consensus.vcf.gz"))


def validate_vcf_header(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        header = "".join(line for line in handle if line.startswith("#"))
    if "ID=ENS_SUPPORT" not in header:
        raise AssertionError(f"ENS_SUPPORT missing from {path}")
    for removed in ("ENS_CONF_LO", "ENS_CONF_HI"):
        if f"ID={removed}" in header:
            raise AssertionError(f"removed field {removed} present in {path}")


def main() -> int:
    args = parse_args()
    rows, names, failed = read_trace(args.trace)
    if not rows:
        raise AssertionError("execution trace is empty")
    if failed:
        preview = [r.get("process") or r.get("name") for r in failed[:5]]
        raise AssertionError(f"failed tasks in trace: {preview}")

    joined = "\n".join(names).upper()
    missing = [
        label for label, patterns in REQUIRED_PROCESS_GROUPS.items()
        if not any(pattern in joined for pattern in patterns)
    ]
    if missing:
        raise AssertionError(f"required seq2neo process groups missing: {missing}")
    forbidden = [name for name in FORBIDDEN_FASTQ_PROCESSES if name in joined]
    if forbidden:
        raise AssertionError(
            f"external-alignment normalization ran for FASTQ-derived BAMs: {forbidden}"
        )

    vcfs = consensus_vcfs(args.outdir)
    if len(vcfs) < 2:
        raise AssertionError(f"expected DNA and RNA consensus VCFs, found {vcfs}")
    for vcf in vcfs:
        validate_vcf_header(vcf)

    rescue = list(args.outdir.glob("rescue/**/*.filtered.vcf.gz"))
    if not rescue:
        raise AssertionError("final filtered rescue VCF is missing")
    print(
        f"seq2neo regression PASS: {len(rows)} tasks, "
        f"{len(vcfs)} consensus VCFs, {len(rescue)} filtered rescue VCFs"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
