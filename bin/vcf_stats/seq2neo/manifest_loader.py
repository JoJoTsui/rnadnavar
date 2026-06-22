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
        "pre_norm_subdir": "variant_calling/mutect2/{prefix}DT_vs_{prefix}DN",
        "pre_norm_pattern": "*.mutect2.vcf.gz",
    },
    "RNA_mutect2": {
        "subdir": "vcf_realignment/normalized/mutect2/{prefix}RT_realign_vs_{prefix}DN",
        "pattern": "*.mutect2.*.dec.norm.vcf.gz",
        "sample_suffix": "RT",
        "format_fields": "GT,AD,AF,DP",
        "pre_norm_subdir": "variant_calling/mutect2/{prefix}RT_vs_{prefix}DN",
        "pre_norm_pattern": "*.mutect2.vcf.gz",
    },
    "DNA_deepsomatic": {
        "subdir": "normalized/deepsomatic/{prefix}DT_vs_{prefix}DN",
        "pattern": "*.deepsomatic.*.dec.norm.vcf.gz",
        "sample_suffix": "DT",
        "format_fields": "GT,AD,VAF,DP",
        "pre_norm_subdir": "variant_calling/deepsomatic/{prefix}DT_vs_{prefix}DN",
        "pre_norm_pattern": "*.deepsomatic.vcf.gz",
    },
    "RNA_deepsomatic": {
        "subdir": "vcf_realignment/normalized/deepsomatic/{prefix}RT_realign_vs_{prefix}DN",
        "pattern": "*.deepsomatic.*.dec.norm.vcf.gz",
        "sample_suffix": "RT",
        "format_fields": "GT,AD,VAF,DP",
        "pre_norm_subdir": "variant_calling/deepsomatic/{prefix}RT_vs_{prefix}DN",
        "pre_norm_pattern": "*.deepsomatic.vcf.gz",
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
    """Load sample manifest from TSV or Parquet."""
    path = Path(path)
    if path.suffix == ".parquet":
        return pl.read_parquet(path)
    return pl.read_csv(path, separator="\t")


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
    manifest_row: dict[str, Any] | None = None,
) -> dict[str, str | None]:
    """Construct all 6 caller VCF paths for a sample.

    When manifest_row is provided and contains pre-computed caller paths,
    they are used directly (no glob). Falls back to glob discovery if the
    manifest column is missing or empty.

    Returns a dict mapping caller name -> VCF path (or None if missing).
    """
    paths = {}
    for caller_name, cfg in CALLER_CONFIGS.items():
        col = f"caller_{caller_name.lower()}"
        # Pre-computed path from manifest
        if manifest_row and col in manifest_row and manifest_row[col]:
            candidate = manifest_row[col]
            if os.path.isfile(candidate):
                paths[caller_name] = candidate
                continue
        # Fallback: glob discovery
        base = os.path.join(base_output_dir, dir_name)
        subdir = cfg["subdir"].format(prefix=vcf_prefix)
        paths[caller_name] = _find_vcf_file(base, subdir, cfg["pattern"])
    return paths


def get_manifest_bam_paths(
    base_output_dir: str,
    dir_name: str,
    manifest_row: dict[str, Any] | None = None,
) -> dict[str, str | None]:
    """Get BAM paths for a sample from manifest or fallback to discovery.

    Manifest columns: bam_dn, bam_dt, bam_rt
    Fallback: delegates to bam_stats._locate_bam_file

    Returns dict mapping BAM type (DN/DT/RT) -> path (or None if missing).
    """
    from .bam_stats import _locate_bam_file

    paths = {}
    for bt in ["DN", "DT", "RT"]:
        col = f"bam_{bt.lower()}"
        if manifest_row and col in manifest_row and manifest_row[col]:
            candidate = manifest_row[col]
            if os.path.isfile(candidate):
                paths[bt] = candidate
                continue
        paths[bt] = _locate_bam_file(base_output_dir, dir_name, bt)
    return paths


def get_pre_norm_vcf_path(
    base_dir: str,
    dir_name: str,
    vcf_prefix: str,
    caller_name: str,
) -> str | None:
    """Find the pre-normalization VCF path for a caller.

    Returns the un-normalized VCF for Mutect2 and DeepSomatic (useful for
    multi-allelic analysis). Returns None for Strelka or when the VCF is
    not found.

    Args:
        base_dir: Base output directory.
        dir_name: Sample directory name.
        vcf_prefix: VCF prefix for the sample.
        caller_name: Caller key from CALLER_CONFIGS.

    Returns:
        VCF path string, or None if not found or caller is Strelka.
    """
    cfg = CALLER_CONFIGS.get(caller_name)
    if cfg is None:
        return None
    if "strelka" in caller_name.lower():
        return None

    pre_subdir = cfg.get("pre_norm_subdir")
    pre_pattern = cfg.get("pre_norm_pattern")
    if not pre_subdir or not pre_pattern:
        return None

    base = os.path.join(base_dir, dir_name)
    subdir = pre_subdir.format(prefix=vcf_prefix)
    return _find_vcf_file(base, subdir, pre_pattern)
