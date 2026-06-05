#!/usr/bin/env python3
"""
Build a unified sample manifest from merged.json and known output directories.

For each eligible sample (partition_set 1-4), discover its rescue VCF path,
verify file existence, and emit a CSV + Parquet manifest that serves as the
single source of truth for downstream statistics tools.

Usage:
    .venv/bin/python examples/seq2neo/build_sample_manifest.py \\
        --merged-json examples/seq2neo/data/processed/merged.json \\
        --output-prefix examples/seq2neo/data/processed/sample_manifest
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path

# ── Known base output directories per set ──────────────────────────────────
SET_TO_BASE_DIR = {
    1: "/t9k/mnt/WorkSpace/data/ngs/zhanlingmin/Rnadnavar/output",
    2: "/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seq2neo/output",
    3: "/t9k/mnt/WorkSpace/data/ngs/liuxin/seq2neo/output",
    4: "/t9k/mnt/WorkSpace/data/ngs/liuxin/seq2neo/output",
}

# Rescue VCF filename pattern (filled with vcf_prefix)
RESCUE_VCF_TEMPLATE = (
    "{vcf_prefix}DT_vs_{vcf_prefix}DN"
    "_rescued_{vcf_prefix}RT_realign_vs_{vcf_prefix}DN"
    ".filtered.vcf.stripped.vcf.gz"
)

# Non-sample directory names to skip
SKIP_DIRS = {"COO8801.shared", ".ipynb_checkpoints"}


def load_merged_json(path: str) -> list[dict]:
    """Load merged.json and return flat list of eligible samples."""
    with open(path) as fh:
        data = json.load(fh)

    samples = []
    for project in data["projects"]:
        for sample in project["samples"]:
            ps = sample.get("partition_set")
            if ps is None:
                continue  # skip incomplete/unassigned
            samples.append(sample)
    return samples


def build_manifest(samples: list[dict]) -> list[dict]:
    """Build manifest rows for every eligible sample."""
    rows = []
    n_complete = 0
    n_incomplete = 0

    for s in samples:
        sample_id = f"{s['_project_id']}_{s['patient_id']}"
        set_number = s["partition_set"]
        patient_id = str(s["patient_id"])
        base_dir = SET_TO_BASE_DIR[set_number]

        # Set 1 has different naming: dir = "{patient_id}_rnadnavar", prefix = "{patient_id}"
        if set_number == 1:
            dir_name = f"{patient_id}_rnadnavar"
            vcf_prefix = patient_id
        else:
            dir_name = sample_id
            vcf_prefix = sample_id

        # Construct rescue VCF path
        rescue_dir = os.path.join(
            base_dir, dir_name, "vcf_realignment", "rescue",
            f"{vcf_prefix}DT_vs_{vcf_prefix}DN_rescued_{vcf_prefix}RT_realign_vs_{vcf_prefix}DN",
        )
        rescue_vcf_filename = RESCUE_VCF_TEMPLATE.format(vcf_prefix=vcf_prefix)
        rescue_vcf_path = os.path.join(rescue_dir, rescue_vcf_filename)
        is_complete = os.path.isfile(rescue_vcf_path)

        if is_complete:
            n_complete += 1
        else:
            n_incomplete += 1

        rows.append({
            "sample_id": sample_id,
            "project_id": s["_project_id"],
            "patient_id": patient_id,
            "set_number": set_number,
            "disease": s["disease"],
            "disease_normalized": s.get("disease_normalized", s["disease"]),
            "status": s["status"],
            "status_reason": s.get("status_reason", ""),
            "base_output_dir": base_dir,
            "dir_name": dir_name,
            "vcf_prefix": vcf_prefix,
            "rescue_vcf_path": rescue_vcf_path,
            "is_complete": is_complete,
        })

    print(f"Manifest: {len(rows)} samples ({n_complete} complete, {n_incomplete} incomplete)")
    return rows


def validate_coverage(samples: list[dict], manifest: list[dict]) -> bool:
    """Cross-check manifest against merged.json samples. Returns True if clean."""
    expected = {
        f"{s['_project_id']}_{s['patient_id']}"
        for s in samples
    }
    actual = {r["sample_id"] for r in manifest}

    missing = expected - actual
    extra = actual - expected

    ok = True
    if missing:
        print(f"WARNING: {len(missing)} samples in sets but missing from manifest:")
        for m in sorted(missing):
            print(f"  - {m}")
        ok = False
    if extra:
        print(f"WARNING: {len(extra)} samples in manifest but not in sets:")
        for e in sorted(extra):
            print(f"  - {e}")
        ok = False

    if ok:
        print(f"Coverage OK: all {len(expected)} expected samples present, no extras.")

    return ok


def write_csv(rows: list[dict], path: str):
    """Write manifest as CSV."""
    if not rows:
        print("No rows to write.")
        return
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"CSV written: {path}")


def write_parquet(rows: list[dict], path: str):
    """Write manifest as Parquet using polars (optional)."""
    try:
        import polars as pl
    except ImportError:
        print("polars not available, skipping Parquet export.")
        return

    df = pl.DataFrame(rows)
    df.write_parquet(path)
    print(f"Parquet written: {path}")


def print_summary(rows: list[dict]):
    """Print per-set completeness summary."""
    print()
    print("=" * 60)
    print("SUMMARY BY SET")
    print("=" * 60)
    by_set = {}
    for r in rows:
        s = r["set_number"]
        if s not in by_set:
            by_set[s] = {"total": 0, "complete": 0, "incomplete": []}
        by_set[s]["total"] += 1
        if r["is_complete"]:
            by_set[s]["complete"] += 1
        else:
            by_set[s]["incomplete"].append(r["sample_id"])

    for set_num in sorted(by_set):
        info = by_set[set_num]
        print(f"Set {set_num}: {info['complete']}/{info['total']} complete")
        if info["incomplete"]:
            for sid in info["incomplete"]:
                print(f"  INCOMPLETE: {sid}")

    # Check non-sample dirs
    print()
    print("=" * 60)
    print("DIRECTORY SCAN (non-sample entries)")
    print("=" * 60)
    all_dirs = set()
    for base_dir in SET_TO_BASE_DIR.values():
        if os.path.isdir(base_dir):
            for entry in os.listdir(base_dir):
                entry_path = os.path.join(base_dir, entry)
                if os.path.isdir(entry_path) and not entry.startswith("."):
                    all_dirs.add(entry)
    manifest_dirs = {r["dir_name"] for r in rows}
    non_sample = all_dirs - manifest_dirs
    if non_sample:
        for d in sorted(non_sample):
            print(f"  Non-sample dir: {d}")
    else:
        print("  No unexpected directories found.")


def main():
    parser = argparse.ArgumentParser(description="Build seq2neo sample manifest")
    parser.add_argument(
        "--merged-json", required=True,
        help="Path to merged.json (e.g., data/processed/merged.json)",
    )
    parser.add_argument(
        "--output-prefix", required=True,
        help="Output prefix for manifest files (without extension)",
    )
    parser.add_argument(
        "--no-parquet", action="store_true",
        help="Skip Parquet output (CSV only)",
    )
    args = parser.parse_args()

    samples = load_merged_json(args.merged_json)
    print(f"Loaded {len(samples)} eligible samples from merged.json")

    rows = build_manifest(samples)
    validate_coverage(samples, rows)
    print_summary(rows)

    csv_path = f"{args.output_prefix}.csv"
    write_csv(rows, csv_path)

    if not args.no_parquet:
        parquet_path = f"{args.output_prefix}.parquet"
        write_parquet(rows, parquet_path)


if __name__ == "__main__":
    main()
