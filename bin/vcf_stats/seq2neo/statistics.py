"""Compute variant statistics from rescue + caller data using polars.

Computes per-variant metrics (VAF, means) and aggregate statistics
(per-sample, per-set, per-disease summaries).
"""

from typing import Any

import polars as pl

from .manifest_loader import CALLER_CONFIGS

# Columns needed by cross-sample aggregation functions. The lazy scan over
# per-sample parquet files has ~165 columns; selecting only these ~35 before
# .collect() reduces memory from ~80 GB to ~15 GB (at 58M variants) and
# enables the query optimizer to skip unneeded columns during parquet reads.
_CROSS_SAMPLE_COLS = [
    # Per-sample metadata
    "sample_id", "set_number", "disease", "disease_normalized",
    # Classification / filters
    "FILTER", "VC", "variant_type", "ti_tv",
    # Computed per-variant means (compute_all_per_variant output)
    "DNA_VAF_mean", "RNA_VAF_mean", "DNA_DP_mean", "RNA_DP_mean",
    "DNA_REF_DP_mean", "RNA_REF_DP_mean", "DNA_ALT_DP_mean", "RNA_ALT_DP_mean",
    # Tiering columns
    "final_tier", "caller_tier", "database_tier", "tier_quality",
    # Caller support / cross-modality / rescue flags
    "N_SUPPORT_CALLERS", "CROSS_MODALITY", "RESCUED",
    # Database annotations
    "COSMIC_ID", "GNOMAD_AF", "REDI_EVIDENCE",
    # GT concordance (caller GT columns — only 4 of 6 callers have GT)
    "DNA_mutect2_GT", "RNA_mutect2_GT", "DNA_deepsomatic_GT", "RNA_deepsomatic_GT",
    # Flag filter breakdown fields
    "min_alt_reads", "gnomad", "blacklist", "noncoding", "ig_pseudo",
    "homopolymer", "vc_filter", "not_consensus", "multiallelic",
]


def _ensure_eager(df):
    """Materialize a LazyFrame if needed, loading only cross-sample columns.

    Avoids loading all 165 columns from the per-sample parquet files.
    Accepts both eager (pass-through) and lazy (collect with column pruning).
    """
    if isinstance(df, pl.LazyFrame):
        existing = [c for c in _CROSS_SAMPLE_COLS if c in df.columns]
        return df.select(existing).collect()
    return df


# ── All 6 caller names ────────────────────────────────────────────────────
DNA_CALLERS = ["DNA_deepsomatic", "DNA_mutect2", "DNA_strelka"]
RNA_CALLERS = ["RNA_deepsomatic", "RNA_mutect2", "RNA_strelka"]
ALL_CALLERS = DNA_CALLERS + RNA_CALLERS


def compute_vaf_columns(df: pl.DataFrame) -> pl.DataFrame:
    """Compute per-caller VAF = AD_ALT / DP.

    For Strelka callers, AD_ALT comes from TAR[0] and AD_REF from TOR[0].
    VAF is set to NaN when DP == 0 or DP is null.
    """
    for caller in ALL_CALLERS:
        ad_col = f"{caller}_AD_ALT"
        dp_col = f"{caller}_DP"
        vaf_col = f"{caller}_VAF"

        if ad_col not in df.columns or dp_col not in df.columns:
            continue

        df = df.with_columns(
            pl.when(pl.col(dp_col).is_not_null() & (pl.col(dp_col) > 0))
            .then(pl.col(ad_col).cast(pl.Float64) / pl.col(dp_col).cast(pl.Float64))
            .otherwise(None)
            .alias(vaf_col)
        )
    return df


def compute_mean_columns(df: pl.DataFrame) -> pl.DataFrame:
    """Compute DNA/RNA mean DP, REF, ALT, VAF across callers.

    Uses row-wise mean; ignores null values.
    """
    # DNA means
    dna_dp_cols = [f"{c}_DP" for c in DNA_CALLERS]
    dna_ref_cols = [f"{c}_AD_REF" for c in DNA_CALLERS]
    dna_alt_cols = [f"{c}_AD_ALT" for c in DNA_CALLERS]
    dna_vaf_cols = [f"{c}_VAF" for c in DNA_CALLERS]

    existing_dna_dp = [c for c in dna_dp_cols if c in df.columns]
    existing_dna_ref = [c for c in dna_ref_cols if c in df.columns]
    existing_dna_alt = [c for c in dna_alt_cols if c in df.columns]
    existing_dna_vaf = [c for c in dna_vaf_cols if c in df.columns]

    if existing_dna_dp:
        df = df.with_columns(
            pl.mean_horizontal(existing_dna_dp).alias("DNA_DP_mean")
        )
    if existing_dna_vaf:
        df = df.with_columns(
            pl.mean_horizontal(existing_dna_vaf).alias("DNA_VAF_mean")
        )
    # REF_DP / ALT_DP means
    if existing_dna_ref:
        df = df.with_columns(
            pl.mean_horizontal(existing_dna_ref).alias("DNA_REF_DP_mean")
        )
    if existing_dna_alt:
        df = df.with_columns(
            pl.mean_horizontal(existing_dna_alt).alias("DNA_ALT_DP_mean")
        )

    # RNA means
    rna_dp_cols = [f"{c}_DP" for c in RNA_CALLERS]
    rna_ref_cols = [f"{c}_AD_REF" for c in RNA_CALLERS]
    rna_alt_cols = [f"{c}_AD_ALT" for c in RNA_CALLERS]
    rna_vaf_cols = [f"{c}_VAF" for c in RNA_CALLERS]

    existing_rna_dp = [c for c in rna_dp_cols if c in df.columns]
    existing_rna_ref = [c for c in rna_ref_cols if c in df.columns]
    existing_rna_alt = [c for c in rna_alt_cols if c in df.columns]
    existing_rna_vaf = [c for c in rna_vaf_cols if c in df.columns]

    if existing_rna_dp:
        df = df.with_columns(
            pl.mean_horizontal(existing_rna_dp).alias("RNA_DP_mean")
        )
    if existing_rna_vaf:
        df = df.with_columns(
            pl.mean_horizontal(existing_rna_vaf).alias("RNA_VAF_mean")
        )
    # RNA REF_DP / ALT_DP means
    if existing_rna_ref:
        df = df.with_columns(
            pl.mean_horizontal(existing_rna_ref).alias("RNA_REF_DP_mean")
        )
    if existing_rna_alt:
        df = df.with_columns(
            pl.mean_horizontal(existing_rna_alt).alias("RNA_ALT_DP_mean")
        )

    return df


def compute_all_per_variant(df: pl.DataFrame) -> pl.DataFrame:
    """Run all per-variant computations: VAF + mean columns."""
    df = compute_vaf_columns(df)
    df = compute_mean_columns(df)
    return df


# ── Aggregate statistics ──────────────────────────────────────────────────


def sample_summary(df: pl.DataFrame, sample_id: str) -> dict[str, Any]:
    """Compute per-sample aggregate statistics."""
    n = len(df)
    if n == 0:
        return {"sample_id": sample_id, "total_variants": 0}

    result = {
        "sample_id": sample_id,
        "total_variants": n,
    }

    # Classification counts — use FILTER column (VC is null in rescue VCF)
    if "FILTER" in df.columns:
        for cat in ["Somatic", "Germline", "Reference", "Artifact", "RNAedit", "NoConsensus"]:
            result[f"n_{cat.lower()}"] = df.filter(pl.col("FILTER") == cat).height

    # Variant type counts
    if "variant_type" in df.columns:
        for vt in ["SNV", "INS", "DEL", "MNV"]:
            result[f"n_{vt}"] = df.filter(pl.col("variant_type") == vt).height

    # Ti/Tv ratio (SNV only; MNV/INDEL have ti_tv=None)
    if "ti_tv" in df.columns:
        ti_count = df.filter(pl.col("ti_tv") == True).height
        tv_count = df.filter(pl.col("ti_tv") == False).height
        excluded = df.filter(pl.col("ti_tv").is_null()).height
        result["ti_count"] = ti_count
        result["tv_count"] = tv_count
        result["ti_tv_ratio"] = ti_count / tv_count if tv_count > 0 else None
        result["ti_tv_excluded"] = excluded  # MNV/INS/DEL not included in ratio

    # Mean VAF (DNA and RNA)
    for vaf_col in ["DNA_VAF_mean", "RNA_VAF_mean"]:
        if vaf_col in df.columns:
            result[f"mean_{vaf_col.lower()}"] = df[vaf_col].mean()

    # Mean DP
    for dp_col in ["DNA_DP_mean", "RNA_DP_mean"]:
        if dp_col in df.columns:
            result[f"mean_{dp_col.lower()}"] = df[dp_col].mean()

    # Mean REF_DP and ALT_DP (DNA and RNA)
    for ref_col in ["DNA_REF_DP_mean", "RNA_REF_DP_mean"]:
        if ref_col in df.columns:
            result[f"mean_{ref_col.lower()}"] = df[ref_col].mean()
    for alt_col in ["DNA_ALT_DP_mean", "RNA_ALT_DP_mean"]:
        if alt_col in df.columns:
            result[f"mean_{alt_col.lower()}"] = df[alt_col].mean()

    # Caller support distribution
    if "N_SUPPORT_CALLERS" in df.columns:
        for c in range(1, 7):
            result[f"n_callers_{c}"] = df.filter(pl.col("N_SUPPORT_CALLERS") == c).height

    # Per-tier variant counts
    if "final_tier" in df.columns:
        tier_counts = df.group_by("final_tier").agg(pl.len().alias("n"))
        for row in tier_counts.to_dicts():
            result[f"tier_{row['final_tier']}"] = row["n"]

    # Cross-modality
    if "CROSS_MODALITY" in df.columns:
        result["n_cross_modality"] = df.filter(pl.col("CROSS_MODALITY") == "YES").height
    if "RESCUED" in df.columns:
        result["n_rescued"] = df.filter(pl.col("RESCUED") == "YES").height

    # COSMIC / gnomAD
    if "COSMIC_ID" in df.columns:
        result["n_cosmic"] = df.filter(pl.col("COSMIC_ID").is_not_null()).height
    if "GNOMAD_AF" in df.columns:
        result["n_gnomad"] = df.filter(pl.col("GNOMAD_AF").is_not_null()).height

    # REDIportal
    if "REDI_EVIDENCE" in df.columns:
        for level in ["HIGH", "MEDIUM", "LOW", "NONE"]:
            result[f"n_redi_{level.lower()}"] = df.filter(pl.col("REDI_EVIDENCE") == level).height

    return result


def set_summary(sample_stats: pl.DataFrame) -> pl.DataFrame:
    """Aggregate per-set statistics from per-sample data."""
    if "set_number" not in sample_stats.columns:
        return pl.DataFrame()

    agg_exprs = [
        pl.col("total_variants").sum().alias("total_variants"),
        pl.col("total_variants").mean().alias("mean_variants_per_sample"),
    ]
    # Aggregate all 6 classification categories if present
    for cat in ["somatic", "germline", "reference", "artifact", "rnaedit", "noconsensus"]:
        col = f"n_{cat}"
        if col in sample_stats.columns:
            agg_exprs.append(pl.col(col).sum().alias(col))

    # Add mean VAF/DP columns if present
    for col in ["mean_dna_vaf_mean", "mean_rna_vaf_mean", "mean_dna_dp_mean", "mean_rna_dp_mean"]:
        if col in sample_stats.columns:
            agg_exprs.append(pl.col(col).mean().alias(f"avg_{col}"))

    # Ti/Tv if present
    for col in ["ti_tv_ratio"]:
        if col in sample_stats.columns:
            agg_exprs.append(pl.col(col).mean().alias(f"avg_{col}"))

    return sample_stats.group_by("set_number").agg(agg_exprs).sort("set_number")


def disease_summary(df: pl.DataFrame | pl.LazyFrame) -> pl.DataFrame:
    """Aggregate statistics per disease with VAF/DP/tier metrics."""
    df = _ensure_eager(df)
    if "disease_normalized" not in df.columns:
        return pl.DataFrame()

    agg_exprs = [
        pl.len().alias("n_variants"),
        pl.col("sample_id").n_unique().alias("n_samples"),
    ]
    for vaf_col in ["DNA_VAF_mean", "RNA_VAF_mean"]:
        if vaf_col in df.columns:
            agg_exprs.append(pl.col(vaf_col).mean().alias(f"mean_{vaf_col.lower()}"))
    for dp_col in ["DNA_DP_mean", "RNA_DP_mean"]:
        if dp_col in df.columns:
            agg_exprs.append(pl.col(dp_col).mean().alias(f"mean_{dp_col.lower()}"))

    return (
        df.group_by("disease_normalized")
        .agg(agg_exprs)
        .sort("n_variants", descending=True)
    )


def variant_type_distribution(df: pl.DataFrame | pl.LazyFrame, group_col: str = "set_number") -> pl.DataFrame:
    """Count SNV/INS/DEL/MNV per group."""
    df = _ensure_eager(df)
    if "variant_type" not in df.columns:
        return pl.DataFrame()
    return (
        df.group_by([group_col, "variant_type"])
        .agg(pl.len().alias("count"))
        .sort([group_col, "variant_type"])
    )


def filter_distribution(df: pl.DataFrame | pl.LazyFrame, group_col: str = "set_number") -> pl.DataFrame:
    """Count variants per FILTER value per group."""
    df = _ensure_eager(df)
    return (
        df.group_by([group_col, "FILTER"])
        .agg(pl.len().alias("count"))
        .sort(group_col, "count", descending=[False, True])
    )


def vc_distribution(df: pl.DataFrame | pl.LazyFrame, group_col: str = "set_number") -> pl.DataFrame:
    """Count variants per VC classification per group."""
    df = _ensure_eager(df)
    if "VC" not in df.columns:
        return pl.DataFrame()
    return (
        df.group_by([group_col, "VC"])
        .agg(pl.len().alias("count"))
        .sort([group_col, "VC"])
    )


def caller_overlap_distribution(df: pl.DataFrame | pl.LazyFrame) -> pl.DataFrame:
    """Distribution of final_tier (C1D1..C7D0) — variant tiering support."""
    df = _ensure_eager(df)
    if "final_tier" not in df.columns:
        return pl.DataFrame()
    return (
        df.group_by("final_tier")
        .agg(pl.len().alias("count"))
        .sort("final_tier")
    )


def caller_support_distribution(df: pl.DataFrame | pl.LazyFrame) -> pl.DataFrame:
    """Distribution of N_SUPPORT_CALLERS (raw caller count, 1-6)."""
    df = _ensure_eager(df)
    if "N_SUPPORT_CALLERS" not in df.columns:
        return pl.DataFrame()
    df = df.with_columns(pl.col("N_SUPPORT_CALLERS").cast(pl.Int64, strict=False))
    return (
        df.group_by("N_SUPPORT_CALLERS")
        .agg(pl.len().alias("count"))
        .sort("N_SUPPORT_CALLERS")
    )


def gt_concordance(df: pl.DataFrame | pl.LazyFrame) -> dict[str, int]:
    df = _ensure_eager(df)
    """Compute GT concordance among 4 callers with GT fields (polars-native).

    Returns counts of variants where 2, 3, or 4 callers have valid GT values,
    plus a count of variants with <2 valid GTs.
    """
    gt_cols = [
        "DNA_mutect2_GT", "RNA_mutect2_GT",
        "DNA_deepsomatic_GT", "RNA_deepsomatic_GT",
    ]
    existing = [c for c in gt_cols if c in df.columns]
    if len(existing) < 2:
        return {}

    invalid = ["./.", "./.", "."]
    valid_df = df.select(existing).with_columns(
        pl.sum_horizontal([
            pl.col(c).is_not_null() & ~pl.col(c).is_in(invalid)
            for c in existing
        ]).alias("_n_valid")
    )
    result = {}
    for a in [2, 3, 4]:
        result[str(a)] = valid_df.filter(pl.col("_n_valid") >= a).height
    result["no_agreement"] = valid_df.filter(pl.col("_n_valid") < 2).height
    return result


def flag_filter_breakdown(df: pl.DataFrame | pl.LazyFrame) -> dict[str, int]:
    df = _ensure_eager(df)
    """Count how many variants have each flag filter set."""
    flags = [
        "min_alt_reads", "gnomad", "blacklist", "noncoding",
        "ig_pseudo", "homopolymer", "vc_filter", "not_consensus", "multiallelic",
    ]
    result = {}
    for flag in flags:
        if flag in df.columns:
            result[flag] = df.filter(pl.col(flag) == True).height
    return result


def dataset_summary(df: pl.DataFrame | pl.LazyFrame) -> dict[str, Any]:
    """Compute whole-dataset aggregate statistics across all samples."""
    df = _ensure_eager(df)
    n = len(df)
    if n == 0:
        return {"total_variants": 0}

    result: dict[str, Any] = {
        "total_variants": n,
        "n_samples": df["sample_id"].n_unique() if "sample_id" in df.columns else 0,
    }

    # Classification — use FILTER column (VC is null in rescue VCF)
    if "FILTER" in df.columns:
        for cat in ["Somatic", "Germline", "Reference", "Artifact", "RNAedit", "NoConsensus"]:
            result[f"n_{cat.lower()}"] = df.filter(pl.col("FILTER") == cat).height

    # Variant types
    if "variant_type" in df.columns:
        for vt in ["SNV", "INS", "DEL", "MNV"]:
            result[f"n_{vt}"] = df.filter(pl.col("variant_type") == vt).height

    # Ti/Tv
    if "ti_tv" in df.columns:
        ti = df.filter(pl.col("ti_tv") == True).height
        tv = df.filter(pl.col("ti_tv") == False).height
        result["ti_count"] = ti
        result["tv_count"] = tv
        result["ti_tv_ratio"] = ti / tv if tv > 0 else None

    # Mean VAF/DP
    for col in ["DNA_VAF_mean", "RNA_VAF_mean", "DNA_DP_mean", "RNA_DP_mean",
                "DNA_REF_DP_mean", "RNA_REF_DP_mean", "DNA_ALT_DP_mean", "RNA_ALT_DP_mean"]:
        if col in df.columns:
            result[f"mean_{col.lower()}"] = df[col].mean()

    # Caller support
    if "N_SUPPORT_CALLERS" in df.columns:
        for c in range(1, 7):
            result[f"n_callers_{c}"] = df.filter(pl.col("N_SUPPORT_CALLERS") == c).height

    # Cross-modality
    if "CROSS_MODALITY" in df.columns:
        result["n_cross_modality"] = df.filter(pl.col("CROSS_MODALITY") == "YES").height
    if "RESCUED" in df.columns:
        result["n_rescued"] = df.filter(pl.col("RESCUED") == "YES").height

    # Database annotations
    if "COSMIC_ID" in df.columns:
        result["n_cosmic"] = df.filter(pl.col("COSMIC_ID").is_not_null()).height
    if "GNOMAD_AF" in df.columns:
        result["n_gnomad"] = df.filter(pl.col("GNOMAD_AF").is_not_null()).height

    # Tier distribution
    if "final_tier" in df.columns:
        tier_counts = df.group_by("final_tier").agg(pl.len().alias("count"))
        for row in tier_counts.to_dicts():
            result[f"tier_{row['final_tier']}"] = row["count"]

    return result


def sample_tier_summary(df: pl.DataFrame | pl.LazyFrame) -> pl.DataFrame:
    df = _ensure_eager(df)
    """Level 4: Per-sample × per-tier aggregate statistics.

    For each (sample_id, final_tier) pair, compute: variant count, mean VAF/DP/
    REF_DP/ALT_DP, variant type distribution, Ti/Tv, and N_SUPPORT_CALLERS dist.
    """
    if "sample_id" not in df.columns or "final_tier" not in df.columns:
        return pl.DataFrame()

    agg_exprs = [pl.len().alias("n_variants")]

    for vaf_col in ["DNA_VAF_mean", "RNA_VAF_mean"]:
        if vaf_col in df.columns:
            agg_exprs.append(pl.col(vaf_col).mean().alias(f"mean_{vaf_col.lower()}"))

    for dp_col in ["DNA_DP_mean", "RNA_DP_mean"]:
        if dp_col in df.columns:
            agg_exprs.append(pl.col(dp_col).mean().alias(f"mean_{dp_col.lower()}"))

    for ref_col in ["DNA_REF_DP_mean", "RNA_REF_DP_mean"]:
        if ref_col in df.columns:
            agg_exprs.append(pl.col(ref_col).mean().alias(f"mean_{ref_col.lower()}"))
    for alt_col in ["DNA_ALT_DP_mean", "RNA_ALT_DP_mean"]:
        if alt_col in df.columns:
            agg_exprs.append(pl.col(alt_col).mean().alias(f"mean_{alt_col.lower()}"))

    if "variant_type" in df.columns:
        for vt in ["SNV", "INS", "DEL", "MNV"]:
            agg_exprs.append((pl.col("variant_type") == vt).sum().alias(f"n_{vt}"))

    if "ti_tv" in df.columns:
        agg_exprs.append(pl.col("ti_tv").sum().alias("n_ti"))
        agg_exprs.append((~pl.col("ti_tv")).sum().alias("n_tv"))

    if "N_SUPPORT_CALLERS" in df.columns:
        for c in range(1, 7):
            agg_exprs.append(
                (pl.col("N_SUPPORT_CALLERS") == c).sum().alias(f"n_callers_{c}")
            )

    return (
        df.group_by(["sample_id", "final_tier"])
        .agg(agg_exprs)
        .sort(["sample_id", "final_tier"])
    )
