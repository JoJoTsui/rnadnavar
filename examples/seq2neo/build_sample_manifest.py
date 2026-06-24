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
    3: "/t9k/mnt/WorkSpace/data/ngs/liuxin/seq2neo/output_1",
    4: "/t9k/mnt/WorkSpace/data/ngs/liuxin/seq2neo/output_1",
}

# Samples with truncated/incomplete BAM files — temporarily excluded from manifest.
# Data team needs to regenerate/re-transfer these before re-enabling.
EXCLUDED_SAMPLES = {
    "PRJNA298376_4077",
    "PRJNA298376_4110",
    "PRJNA298376_4200",
    "PRJNA298376_4220",
    "PRJNA298376_4264",
}

# Rescue VCF filename pattern (filled with vcf_prefix)
RESCUE_VCF_TEMPLATE = (
    "{vcf_prefix}DT_vs_{vcf_prefix}DN"
    "_rescued_{vcf_prefix}RT_realign_vs_{vcf_prefix}DN"
    ".filtered.vcf.stripped.vcf.gz"
)

# Non-sample directory names to skip
SKIP_DIRS = {"COO8801.shared", ".ipynb_checkpoints"}

# ── Caller normalized VCF path config (mirrors manifest_loader.CALLER_CONFIGS) ─
CALLER_CONFIGS = {
    "DNA_mutect2": {
        "subdir": "normalized/mutect2/{prefix}DT_vs_{prefix}DN",
        "pattern": "*.mutect2.*.dec.norm.vcf.gz",
    },
    "RNA_mutect2": {
        "subdir": "vcf_realignment/normalized/mutect2/{prefix}RT_realign_vs_{prefix}DN",
        "pattern": "*.mutect2.*.dec.norm.vcf.gz",
    },
    "DNA_deepsomatic": {
        "subdir": "normalized/deepsomatic/{prefix}DT_vs_{prefix}DN",
        "pattern": "*.deepsomatic.*.dec.norm.vcf.gz",
    },
    "RNA_deepsomatic": {
        "subdir": "vcf_realignment/normalized/deepsomatic/{prefix}RT_realign_vs_{prefix}DN",
        "pattern": "*.deepsomatic.*.dec.norm.vcf.gz",
    },
    "DNA_strelka": {
        "subdir": "normalized/strelka/{prefix}DT_vs_{prefix}DN",
        "pattern": "*.strelka.*.dec.norm.vcf.gz",
    },
    "RNA_strelka": {
        "subdir": "vcf_realignment/normalized/strelka/{prefix}RT_realign_vs_{prefix}DN",
        "pattern": "*.strelka.*.dec.norm.vcf.gz",
    },
}

# BAM type config (mirrors bam_stats.BAM_TYPES)
BAM_TYPES = {
    "DN": {"suffix": "DN", "label": "DNA Normal"},
    "DT": {"suffix": "DT", "label": "DNA Tumor"},
    "RT": {"suffix": "RT", "label": "RNA Tumor"},
}


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


def _locate_bam_file(base_dir: str, dir_name: str, bam_type: str) -> str | None:
    """Locate BAM file for a given BAM type (DN, DT, or RT)."""
    import glob

    suffix = BAM_TYPES[bam_type]["suffix"]

    if bam_type in ("DN", "DT"):
        search_paths = [
            os.path.join(base_dir, dir_name, "preprocessing", "mapped", f"*{suffix}", "*.sorted.bam"),
            os.path.join(base_dir, dir_name, "preprocessing", "mapped", f"*{suffix}", "*.bam"),
        ]
    else:  # RT
        search_paths = [
            os.path.join(base_dir, dir_name, "preprocessing", "mapped", f"*{suffix}", "*.bam"),
            os.path.join(base_dir, dir_name, "vcf_realignment", "preprocessing", "mapped", f"*{suffix}", "*.bam"),
        ]

    for pattern in search_paths:
        files = glob.glob(pattern)
        if files:
            return files[0]
    return None


def _locate_caller_vcf(base_dir: str, dir_name: str, vcf_prefix: str,
                       caller_name: str, cfg: dict) -> str | None:
    """Locate a single caller normalized VCF file. Returns path or None."""
    import glob

    subdir = cfg["subdir"].format(prefix=vcf_prefix)
    search_path = os.path.join(base_dir, dir_name, subdir, cfg["pattern"])
    files = glob.glob(search_path)
    if files:
        return files[0]
    return None


def build_manifest(samples: list[dict]) -> list[dict]:
    """Build manifest rows for every eligible sample."""
    rows = []
    n_complete = 0
    n_incomplete = 0

    n_excluded = 0
    for s in samples:
        sample_id = f"{s['_project_id']}_{s['patient_id']}"
        if sample_id in EXCLUDED_SAMPLES:
            n_excluded += 1
            continue
        set_number = s.get("partition_set")
        if set_number is None:
            print(f"WARNING: {sample_id} has no partition_set — defaulting to 0")
            set_number = 0
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

        # BAM paths — locate all three BAM types per sample
        bam_paths = {}
        for bt in ["DN", "DT", "RT"]:
            bam_path = _locate_bam_file(base_dir, dir_name, bt)
            bam_paths[f"bam_{bt.lower()}"] = bam_path or ""

        # Caller normalized VCF paths — locate all 6 callers
        caller_paths = {}
        for caller_name, cfg in CALLER_CONFIGS.items():
            vcf_path = _locate_caller_vcf(base_dir, dir_name, vcf_prefix, caller_name, cfg)
            caller_paths[f"caller_{caller_name.lower()}"] = vcf_path or ""

        # Completeness: rescue VCF + at least one BAM + at least one caller VCF
        has_rescue = os.path.isfile(rescue_vcf_path)
        has_bam = any(v for v in bam_paths.values())
        has_caller = any(v for v in caller_paths.values())
        is_complete = has_rescue and has_bam and has_caller

        if is_complete:
            n_complete += 1
        else:
            n_incomplete += 1

        row = {
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
        }
        row.update(bam_paths)
        row.update(caller_paths)
        rows.append(row)

    print(f"Manifest: {len(rows)} samples ({n_complete} complete, {n_incomplete} incomplete, {n_excluded} excluded)")
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


def write_tsv(rows: list[dict], path: str):
    """Write manifest as TSV (tab-separated)."""
    if not rows:
        print("No rows to write.")
        return
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"TSV written: {path}")


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
        help="Skip Parquet output (TSV only)",
    )
    args = parser.parse_args()

    samples = load_merged_json(args.merged_json)
    print(f"Loaded {len(samples)} eligible samples from merged.json")

    rows = build_manifest(samples)
    validate_coverage(samples, rows)
    print_summary(rows)

    tsv_path = f"{args.output_prefix}.tsv"
    write_tsv(rows, tsv_path)

    if not args.no_parquet:
        parquet_path = f"{args.output_prefix}.parquet"
        write_parquet(rows, parquet_path)


if __name__ == "__main__":
    main()
