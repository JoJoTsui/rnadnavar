"""Parse rescue VCF files for metadata, annotations, and classification fields.

Uses cyvcf2 for fast VCF parsing. Extracts INFO fields into a polars DataFrame.
Rescue VCF fields are used for classification/metadata — numeric QC metrics
(DP, AD, VAF) come from caller VCFs (caller_parser.py) as the ground truth.
"""

from typing import Any

import numpy as np
import polars as pl
from cyvcf2 import VCF

# ── INFO fields to extract from the rescue VCF ────────────────────────────

# Numeric (int) fields
RESCUE_INT_FIELDS = [
    "N_SUPPORT_CALLERS",
    "N_CONSENSUS_SUPPORT",
    "N_DNA_CALLERS_SUPPORT",
    "N_RNA_CALLERS_SUPPORT",
    "COSMIC_CNT",
    "DP_MIN",
    "DP_MAX",
]

# Numeric (float) fields
RESCUE_FLOAT_FIELDS = [
    "GNOMAD_AF",
    "QUAL_MEAN",
    "QUAL_MIN",
    "QUAL_MAX",
    "DP_MEAN",
    "VAF_MEAN",
    "VAF_MIN",
    "VAF_MAX",
    "DP_DNA_MEAN",
    "DP_RNA_MEAN",
    "VAF_DNA_MEAN",
    "VAF_RNA_MEAN",
]

# String fields
RESCUE_STRING_FIELDS = [
    "CALLERS_SUPPORT",
    "CALLERS",
    "CONSENSUS_SUPPORT",
    "MODALITIES",
    "CALLERS_BY_MODALITY",
    "CROSS_MODALITY",
    "VC",
    "VC_CALLERS",
    "VC_CONSENSUS",
    "UNIFIED_FILTER",
    "UNIFIED_FILTER_DNA",
    "UNIFIED_FILTER_RNA",
    "PASSES_CONSENSUS",
    "PASSES_CONSENSUS_DNA",
    "PASSES_CONSENSUS_RNA",
    "RESCUED",
    "GNOMAD_RESCUE",
    "COSMIC_RESCUE",
    "FILTERS_ORIGINAL",
    "FILTERS_NORMALIZED",
    "FILTERS_CATEGORY",
    # Individual per-caller per-modality filter fields for tiering
    "FILTER_NORMALIZED_Strelka_DNA_TUMOR",
    "FILTER_NORMALIZED_Strelka_RNA_TUMOR",
    "FILTER_NORMALIZED_Mutect2_DNA_TUMOR",
    "FILTER_NORMALIZED_Mutect2_RNA_TUMOR",
    "FILTER_NORMALIZED_DeepSomatic_DNA_TUMOR",
    "FILTER_NORMALIZED_DeepSomatic_RNA_TUMOR",
    "CONSENSUS_GT",
    "GT_BY_CALLER",
    "DP_BY_CALLER",
    "VAF_BY_CALLER",
    "COSMIC_ID",
    "REDI_ACCESSION",
    "REDI_DB",
    "REDI_TYPE",
    "REDI_REPEAT",
    "REDI_STRAND",
    "REDI_FUNC",
    "REDI_EVIDENCE",
    "REDI_CANONICAL",
    "OLD_MULTIALLELIC",
]

# Boolean flag fields (present = True)
RESCUE_FLAG_FIELDS = [
    "min_alt_reads",
    "gnomad",
    "blacklist",
    "noncoding",
    "ig_pseudo",
    "homopolymer",
    "vc_filter",
    "not_consensus",
    "multiallelic",
    "SOMATIC",  # Strelka-specific flag in combined VCF
]

ALL_RESCUE_FIELDS = (
    RESCUE_INT_FIELDS + RESCUE_FLOAT_FIELDS + RESCUE_STRING_FIELDS + RESCUE_FLAG_FIELDS
)


def rescue_info_fields() -> list[str]:
    """Return the list of all INFO fields extracted from rescue VCFs."""
    return list(ALL_RESCUE_FIELDS)


def _safe_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (ValueError, TypeError):
        return None


def _safe_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        v = float(value)
        if np.isnan(v) or np.isinf(v):
            return None
        return v
    except (ValueError, TypeError):
        return None


def _safe_str(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return "|".join(str(v) for v in value)
    return str(value)


def _add_derived_columns_polars(df: pl.DataFrame) -> pl.DataFrame:
    """Compute variant_type and ti_tv using polars vectorized expressions.

    Replaces the old .to_list() + Python loop approach. 6× faster.
    Used by both the cyvcf2 fallback path and the Rust row-oriented fallback.
    """
    # variant_type: classify as SNV, INS, DEL, or MNV
    df = df.with_columns(
        pl.when(pl.col("ALT").str.len_chars() > pl.col("REF").str.len_chars())
        .then(pl.lit("INS"))
        .when(pl.col("REF").str.len_chars() > pl.col("ALT").str.len_chars())
        .then(pl.lit("DEL"))
        .when((pl.col("REF").str.len_chars() == 1) & (pl.col("ALT").str.len_chars() == 1))
        .then(pl.lit("SNV"))
        .otherwise(pl.lit("MNV"))
        .alias("variant_type")
    )
    # ti_tv: True for transitions (A↔G, C↔T), False for transversions, None for non-SNV
    transitions = ["AG", "GA", "CT", "TC"]
    df = df.with_columns(
        pl.when(
            (pl.col("variant_type") == "SNV")
            & (pl.col("REF").str.to_uppercase() + pl.col("ALT").str.to_uppercase()).is_in(transitions)
        )
        .then(pl.lit(True))
        .when(pl.col("variant_type") == "SNV")
        .then(pl.lit(False))
        .otherwise(pl.lit(None))
        .alias("ti_tv")
    )
    return df


def parse_rescue_vcf(vcf_path: str) -> pl.DataFrame:
    """Parse a rescue VCF file and return a polars DataFrame.

    Extracts core coordinates, FILTER status, all INFO fields listed in
    rescue_info_fields(), plus derived columns variant_type and ti_tv.

    Args:
        vcf_path: Path to the .filtered.vcf.stripped.vcf.gz rescue VCF.

    Returns:
        polars DataFrame with one row per variant.
    """
    records = []

    reader = VCF(vcf_path)
    for record in reader:
        chrom = record.CHROM
        pos = record.POS
        ref = record.REF
        alt = record.ALT[0] if record.ALT else None
        if alt is None:
            continue
        alt_str = str(alt)
        filter_val = (
            record.FILTERS[0]
            if record.FILTERS and record.FILTERS[0] is not None
            else "PASS"
        )
        info = record.INFO

        row: dict[str, Any] = {
            "CHROM": chrom,
            "POS": pos,
            "REF": ref,
            "ALT": alt_str,
            "FILTER": str(filter_val) if filter_val else "PASS",
        }

        # Extract int fields
        for field in RESCUE_INT_FIELDS:
            row[field] = _safe_int(info.get(field))

        # Extract float fields
        for field in RESCUE_FLOAT_FIELDS:
            row[field] = _safe_float(info.get(field))

        # Extract string fields
        for field in RESCUE_STRING_FIELDS:
            row[field] = _safe_str(info.get(field))

        # Extract flag fields (True if present in INFO)
        for field in RESCUE_FLAG_FIELDS:
            row[field] = info.get(field) is True

        records.append(row)

    if not records:
        return pl.DataFrame()

    df = pl.from_dicts(records, infer_schema_length=None)

    # Ensure int columns are int
    for field in RESCUE_INT_FIELDS:
        if field in df.columns:
            df = df.with_columns(pl.col(field).cast(pl.Int64, strict=False))

    # Ensure float columns are float
    for field in RESCUE_FLOAT_FIELDS:
        if field in df.columns:
            df = df.with_columns(pl.col(field).cast(pl.Float64, strict=False))

    # Ensure flag columns are bool
    for field in RESCUE_FLAG_FIELDS:
        if field in df.columns:
            df = df.with_columns(pl.col(field).cast(pl.Boolean, strict=False))

    # Derived columns via polars vectorized expressions (6× faster than per-row Python)
    df = _add_derived_columns_polars(df)

    return df
