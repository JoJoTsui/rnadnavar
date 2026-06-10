"""Load and validate the seq2neo sample manifest produced by build_sample_manifest.py.

Provides utilities for loading the manifest CSV/Parquet, determining VCF prefixes
per sample, filtering to complete samples, and constructing caller VCF paths.
"""

import os
from pathlib import Path
from typing import Any

import polars as pl

# ── Caller VCF path configurations ────────────────────────────────────────

CALLER_CONFIGS: dict[str, dict[str, str]] = {
    "DNA_mutect2": {
        "subdir": "normalized/mutect2/{prefix}DT_vs_{prefix}DN",
        "pattern": "*.mutect2.*.dec.norm.vcf.gz",
        "sample_suffix": "DT",
        "format_fields": "GT,AD,AF,DP",
    },
    "RNA_mutect2": {
        "subdir": "vcf_realignment/normalized/mutect2/{prefix}RT_realign_vs_{prefix}DN",
        "pattern": "*.mutect2.*.dec.norm.vcf.gz",
        "sample_suffix": "RT",
        "format_fields": "GT,AD,AF,DP",
    },
    "DNA_deepsomatic": {
        "subdir": "normalized/deepsomatic/{prefix}DT_vs_{prefix}DN",
        "pattern": "*.deepsomatic.*.dec.norm.vcf.gz",
        "sample_suffix": "DT",
        "format_fields": "GT,AD,VAF,DP",
    },
    "RNA_deepsomatic": {
        "subdir": "vcf_realignment/normalized/deepsomatic/{prefix}RT_realign_vs_{prefix}DN",
        "pattern": "*.deepsomatic.*.dec.norm.vcf.gz",
        "sample_suffix": "RT",
        "format_fields": "GT,AD,VAF,DP",
    },
    "DNA_strelka": {
        "subdir": "normalized/strelka/{prefix}DT_vs_{prefix}DN",
        "pattern": "*.strelka.*.dec.norm.vcf.gz",
        "sample_suffix": "TUMOR",
        "format_fields": "DP,TAR,TIR,TOR",
    },
    "RNA_strelka": {
        "subdir": "vcf_realignment/normalized/strelka/{prefix}RT_realign_vs_{prefix}DN",
        "pattern": "*.strelka.*.dec.norm.vcf.gz",
        "sample_suffix": "TUMOR",
        "format_fields": "DP,TAR,TIR,TOR",
    },
}


def load_manifest(path: str | Path) -> pl.DataFrame:
    """Load sample manifest from CSV or Parquet."""
    path = Path(path)
    if path.suffix == ".parquet":
        return pl.read_parquet(path)
    return pl.read_csv(path)


def get_vcf_prefix(set_number: int, patient_id: str, sample_id: str) -> str:
    """Return the VCF prefix for a sample.

    Set 1: patient_id only (e.g., "4060")
    Sets 2-4: full sample_id (e.g., "PRJNA298376_4060")
    """
    if set_number == 1:
        return str(patient_id)
    return str(sample_id)


def get_vcf_prefix_from_row(row: dict[str, Any]) -> str:
    """Return VCF prefix from a manifest row dict."""
    return get_vcf_prefix(
        row["set_number"], str(row["patient_id"]), str(row["sample_id"])
    )


def filter_complete(df: pl.DataFrame) -> pl.DataFrame:
    """Filter manifest to only complete samples."""
    return df.filter(pl.col("is_complete"))


def _find_vcf_file(base_dir: str, subdir: str, pattern: str) -> str | None:
    """Glob for a single VCF file in a subdirectory. Returns first match."""
    import glob

    search_path = os.path.join(base_dir, subdir, pattern)
    files = glob.glob(search_path)
    if files:
        return files[0]
    return None


def get_all_caller_vcf_paths(
    base_output_dir: str,
    dir_name: str,
    vcf_prefix: str,
) -> dict[str, str | None]:
    """Construct all 6 caller VCF paths for a sample.

    Returns a dict mapping caller name -> VCF path (or None if missing).
    """
    base = os.path.join(base_output_dir, dir_name)
    paths = {}
    for caller_name, cfg in CALLER_CONFIGS.items():
        subdir = cfg["subdir"].format(prefix=vcf_prefix)
        paths[caller_name] = _find_vcf_file(base, subdir, cfg["pattern"])
    return paths
