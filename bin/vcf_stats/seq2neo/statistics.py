"""Compute variant statistics from rescue + caller data using polars.

Computes per-variant metrics (VAF, means) and aggregate statistics
(per-sample, per-set, per-disease summaries).
"""

import re
from pathlib import Path
from typing import Any

import polars as pl

from .manifest_loader import CALLER_CONFIGS


def write_tsv(df: pl.DataFrame, path: str | Path) -> None:
    """Write a polars DataFrame to a TSV (tab-separated) file.

    Uses .write_csv(separator="\t") — polars handles the tab separator
    without issues. TSV avoids column corruption from embedded commas in
    VCF pipe/colon-delimited INFO fields and is natively readable by
    pandas/polars with sep="\t".
    """
    df.write_csv(str(path), separator="\t")


# ═══════════════════════════════════════════════════════════════════════════════
# Wise kernel — unified multi-dimension aggregation
# ═══════════════════════════════════════════════════════════════════════════════

def _metric_columns(expr) -> set[str]:
    """Extract column names referenced in a polars expression string."""
    return set(re.findall(r'col\("([^"]+)"\)', str(expr)))


# Shared metrics applied across all wise dimensions.
# Each entry: (name, polars_expression)
_WISE_METRICS: list[tuple[str, pl.Expr]] = [
    ("n_variants", pl.len()),
    ("n_somatic", (pl.col("FILTER") == "Somatic").cast(pl.Int64).sum()),
    ("n_germline", (pl.col("FILTER") == "Germline").cast(pl.Int64).sum()),
    ("n_reference", (pl.col("FILTER") == "Reference").cast(pl.Int64).sum()),
    ("n_artifact", (pl.col("FILTER") == "Artifact").cast(pl.Int64).sum()),
    ("n_rnaedit", (pl.col("FILTER") == "RNAedit").cast(pl.Int64).sum()),
    ("n_noconsensus", (pl.col("FILTER") == "NoConsensus").cast(pl.Int64).sum()),
    ("n_snv", (pl.col("variant_type") == "SNV").cast(pl.Int64).sum()),
    ("n_ins", (pl.col("variant_type") == "INS").cast(pl.Int64).sum()),
    ("n_del", (pl.col("variant_type") == "DEL").cast(pl.Int64).sum()),
    ("n_mnv", (pl.col("variant_type") == "MNV").cast(pl.Int64).sum()),
    ("mean_dna_vaf", pl.col("DNA_VAF_mean").mean()),
    ("mean_rna_vaf", pl.col("RNA_VAF_mean").mean()),
    ("median_dna_vaf", pl.col("DNA_VAF_mean").median()),
    ("median_rna_vaf", pl.col("RNA_VAF_mean").median()),
    ("mean_dna_dp", pl.col("DNA_DP_mean").mean()),
    ("mean_rna_dp", pl.col("RNA_DP_mean").mean()),
    ("n_cross_modality", (pl.col("modality_evidence_caller") == "cross_modality").cast(pl.Int64).sum()),
    ("n_dna_confident", (pl.col("modality_evidence_caller") == "dna_confident").cast(pl.Int64).sum()),
    ("n_rna_rescued", (pl.col("modality_evidence_caller") == "rna_rescued").cast(pl.Int64).sum()),
    ("n_low_confidence", (pl.col("modality_evidence_caller") == "low_confidence").cast(pl.Int64).sum()),
    ("n_rescued", ((pl.col("modality_evidence_caller") == "cross_modality") | (pl.col("modality_evidence_caller") == "rna_rescued")).cast(pl.Int64).sum()),
    ("n_cosmic", pl.col("COSMIC_ID").is_not_null().cast(pl.Int64).sum()),
    ("n_gnomad", pl.col("GNOMAD_AF").is_not_null().cast(pl.Int64).sum()),
    ("n_ti", pl.col("ti_tv").cast(pl.Int64).sum()),
    ("n_tv", (pl.col("ti_tv") == False).cast(pl.Int64).sum()),
    ("mean_n_support_callers", pl.col("N_SUPPORT_CALLERS").mean()),
]

# Additional per-caller metrics for caller-wise aggregation
_CALLER_WISE_METRICS: list[tuple[str, pl.Expr]] = [
    ("mean_vaf", pl.col("VAF").mean()),
    ("median_vaf", pl.col("VAF").median()),
    ("mean_dp", pl.col("DP").mean()),
    ("median_dp", pl.col("DP").median()),
]


def compute_wise_summary(df, group_cols, extra_metrics=None):
    """Generic kernel for wise-based aggregation.

    Applies _WISE_METRICS (shared across all wises) to a group_by
    aggregation. Only metrics whose columns exist in the dataframe
    are included — missing columns are silently skipped.

    Args:
        df: DataFrame or LazyFrame with variant data.
        group_cols: List of column names to group by.
        extra_metrics: Optional list of (name, expr) tuples for
            wise-specific metrics beyond the shared kernel.

    Returns:
        polars DataFrame grouped by group_cols with metric columns.
    """
    df = _ensure_eager(df)
    df_cols = set(df.columns)

    # Filter to metrics whose columns are available
    metrics = []
    for name, expr in _WISE_METRICS:
        if _metric_columns(expr).issubset(df_cols):
            metrics.append(expr.alias(name))

    if extra_metrics:
        for name, expr in extra_metrics:
            if _metric_columns(expr).issubset(df_cols):
                metrics.append(expr.alias(name))

    if not metrics:
        return pl.DataFrame()

    return df.group_by(group_cols).agg(metrics).sort(group_cols)


def compute_caller_wise_summary(df):
    """Caller-wise aggregation — per-caller VAF/DP metrics.

    Melts per-caller columns (DNA_mutect2_VAF, RNA_strelka_DP, etc.)
    into a long-form caller × metric summary.
    """
    df = _ensure_eager(df)
    callers = ["DNA_mutect2", "DNA_deepsomatic", "DNA_strelka",
               "RNA_mutect2", "RNA_deepsomatic", "RNA_strelka"]

    rows = []
    for caller in callers:
        vaf_col = f"{caller}_VAF"
        dp_col = f"{caller}_DP"
        row = {"caller": caller, "n_with_vaf": 0, "n_with_dp": 0,
               "mean_vaf": None, "median_vaf": None,
               "mean_dp": None, "median_dp": None}

        if vaf_col in df.columns:
            vaf_vals = df[vaf_col].drop_nulls()
            row["n_with_vaf"] = len(vaf_vals)
            if len(vaf_vals) > 0:
                row["mean_vaf"] = vaf_vals.mean()
                row["median_vaf"] = vaf_vals.median()

        if dp_col in df.columns:
            dp_vals = df[dp_col].drop_nulls()
            row["n_with_dp"] = len(dp_vals)
            if len(dp_vals) > 0:
                row["mean_dp"] = dp_vals.mean()
                row["median_dp"] = dp_vals.median()

        rows.append(row)

    return pl.DataFrame(rows)

# Columns needed by cross-sample aggregation functions. The lazy scan over
# per-sample parquet files has ~165 columns; selecting only these ~35 before
# .collect() reduces memory from ~80 GB to ~15 GB (at 58M variants) and
# enables the query optimizer to skip unneeded columns during parquet reads.
_CROSS_SAMPLE_COLS = [
    # Per-sample metadata
    "sample_id", "set_number", "disease", "disease_normalized",
    # Coordinates (needed for chromosome-wise summaries)
    "CHROM",
    # Classification / filters
    "FILTER", "VC", "variant_type", "ti_tv",
    # Computed per-variant means (compute_all_per_variant output)
    "DNA_VAF_mean", "RNA_VAF_mean", "DNA_DP_mean", "RNA_DP_mean",
    "DNA_REF_DP_mean", "RNA_REF_DP_mean", "DNA_ALT_DP_mean", "RNA_ALT_DP_mean",
    # Per-caller VAF columns (6) — needed for caller-wise stats, threshold sweep
    "DNA_mutect2_VAF", "RNA_mutect2_VAF",
    "DNA_deepsomatic_VAF", "RNA_deepsomatic_VAF",
    "DNA_strelka_VAF", "RNA_strelka_VAF",
    # Per-caller DP columns (6) — needed for caller-wise depth stats
    "DNA_mutect2_DP", "RNA_mutect2_DP",
    "DNA_deepsomatic_DP", "RNA_deepsomatic_DP",
    "DNA_strelka_DP", "RNA_strelka_DP",
    # Per-caller AD columns (12) — needed for REF/ALT analysis per caller
    "DNA_mutect2_AD_REF", "DNA_mutect2_AD_ALT",
    "RNA_mutect2_AD_REF", "RNA_mutect2_AD_ALT",
    "DNA_deepsomatic_AD_REF", "DNA_deepsomatic_AD_ALT",
    "RNA_deepsomatic_AD_REF", "RNA_deepsomatic_AD_ALT",
    "DNA_strelka_AD_REF", "DNA_strelka_AD_ALT",
    "RNA_strelka_AD_REF", "RNA_strelka_AD_ALT",
    # Tiering columns
    "final_tier", "caller_tier", "database_tier", "tier_quality",
    # Caller support / cross-modality / rescue flags
    "N_SUPPORT_CALLERS", "N_DNA_CALLERS_SUPPORT", "N_RNA_CALLERS_SUPPORT",
    "modality_evidence_caller", "modality_evidence_dp",
    "n_alleles_at_site", "vaf_sum", "allele_balance_ratio",
    "category_conflict", "multiallelic_class",
    "flag_vaf_overflow", "flag_multi_allelic_heterogeneity",
    "flag_category_conflict", "flag_germline_low_vaf",
    "flag_somatic_high_vaf", "flag_reference_with_signal",
    "flag_rna_rescued",
    # Database annotations
    "COSMIC_ID", "GNOMAD_AF", "REDI_EVIDENCE",
    # GT concordance (caller GT columns — only 4 of 6 callers have GT)
    "DNA_mutect2_GT", "RNA_mutect2_GT", "DNA_deepsomatic_GT", "RNA_deepsomatic_GT",
    # Flag filter breakdown fields
    "min_alt_reads", "gnomad", "blacklist", "noncoding", "ig_pseudo",
    "homopolymer", "vc_filter", "not_consensus",
    # ML partition column (added by cli.py from CHROM → train/val/test)
    "partition",
]


def _ensure_eager(df):
    """Materialize a LazyFrame if needed, loading only cross-sample columns.

    Avoids loading all 165 columns from the per-sample parquet files.
    Accepts both eager (pass-through) and lazy (collect with column pruning).
    Uses collect_schema().names() to avoid the PerformanceWarning triggered
    by accessing .columns on a LazyFrame (which requires full schema resolution).
    """
    if isinstance(df, pl.LazyFrame):
        schema_names = df.collect_schema().names()
        existing = [c for c in _CROSS_SAMPLE_COLS if c in schema_names]
        return df.select(existing).collect()
    return df


# ── All 6 caller names ────────────────────────────────────────────────────
DNA_CALLERS = ["DNA_deepsomatic", "DNA_mutect2", "DNA_strelka"]
RNA_CALLERS = ["RNA_deepsomatic", "RNA_mutect2", "RNA_strelka"]
ALL_CALLERS = DNA_CALLERS + RNA_CALLERS


# Per-caller VAF denominator metadata: documents which depth metric each
# caller uses as the VAF denominator. This matters because Strelka and
# Mutect2/DeepSomatic compute VAF against different depth definitions.
#
#   Strelka:      VAF = TAR[0] / DP_tier1
#     DP = FORMAT/DP which is tier-1 filtered depth (reads passing all
#     internal filters including Q13, minTier=1, no duplicates, no filtered
#     sites). Excludes low-quality and artifact reads. Consequently Strelka
#     VAF tends to be higher than Mutect2/DeepSomatic VAF at the same site
#     because the denominator is smaller.
#
#   Mutect2:      VAF = AD[1] / DP_total
#   DeepSomatic:  VAF = AD[1] / DP_total (or FORMAT/VAF when available)
#     DP = FORMAT/DP which is total unfiltered depth. Includes all reads
#     passing mapping quality filters but not tier-1 filtering.
#
# The difference is CORRECT and EXPECTED — each caller defines VAF against
# its own internal depth metric. Cross-caller VAF comparisons must account
# for this denominator difference.
CALLER_VAF_DENOMINATOR = {
    "DNA_strelka": "tier1_depth",
    "RNA_strelka": "tier1_depth",
    "DNA_mutect2": "total_depth",
    "RNA_mutect2": "total_depth",
    "DNA_deepsomatic": "total_depth",
    "RNA_deepsomatic": "total_depth",
}


def compute_vaf_columns(df: pl.DataFrame) -> pl.DataFrame:
    """Compute per-caller VAF = AD_ALT / DP.

    For Strelka callers, AD_ALT comes from TIR[0] and AD_REF from TAR[0],
    and DP is tier-1 filtered depth (see CALLER_VAF_DENOMINATOR above).
    For Mutect2/DeepSomatic, AD_ALT comes from AD[1] and DP is total depth.
    VAF is set to NaN when DP == 0 or DP is null.
    Values are defensively clipped to [0, 1] to guard against edge cases
    (e.g., rounding, multiallelic sites).
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
            .clip(0.0, 1.0)
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


def compute_modality_evidence(df: pl.DataFrame) -> pl.DataFrame:
    """Compute modality_evidence_caller and modality_evidence_dp.

    modality_evidence_caller — classifies each variant based on caller support:
        cross_modality: both DNA ≥ 1 callers AND RNA ≥ 1 callers support
        dna_confident:  DNA ≥ 2 callers support, RNA < 1
        rna_rescued:    RNA ≥ 2 callers support, DNA < 1
        low_confidence:  everything else

    modality_evidence_dp — classifies based on depth evidence (>30 reads):
        cross_modality: both DNA_DP_mean > 30 AND RNA_DP_mean > 30
        dna_confident:  DNA_DP_mean > 30, RNA_DP_mean <= 30
        rna_rescued:    RNA_DP_mean > 30, DNA_DP_mean <= 30
        low_confidence:  both <= 30
    """
    has_dna_callers = "N_DNA_CALLERS_SUPPORT" in df.columns
    has_rna_callers = "N_RNA_CALLERS_SUPPORT" in df.columns
    has_dna_dp = "DNA_DP_mean" in df.columns
    has_rna_dp = "RNA_DP_mean" in df.columns

    if has_dna_callers and has_rna_callers:
        df = df.with_columns(
            pl.when((pl.col("N_DNA_CALLERS_SUPPORT") >= 1) & (pl.col("N_RNA_CALLERS_SUPPORT") >= 1))
            .then(pl.lit("cross_modality"))
            .when((pl.col("N_DNA_CALLERS_SUPPORT") >= 2) & (pl.col("N_RNA_CALLERS_SUPPORT") < 1))
            .then(pl.lit("dna_confident"))
            .when((pl.col("N_RNA_CALLERS_SUPPORT") >= 2) & (pl.col("N_DNA_CALLERS_SUPPORT") < 1))
            .then(pl.lit("rna_rescued"))
            .otherwise(pl.lit("low_confidence"))
            .alias("modality_evidence_caller")
        )

    if has_dna_dp and has_rna_dp:
        df = df.with_columns(
            pl.when((pl.col("DNA_DP_mean") > 30) & (pl.col("RNA_DP_mean") > 30))
            .then(pl.lit("cross_modality"))
            .when((pl.col("DNA_DP_mean") > 30) & (pl.col("RNA_DP_mean") <= 30))
            .then(pl.lit("dna_confident"))
            .when((pl.col("RNA_DP_mean") > 30) & (pl.col("DNA_DP_mean") <= 30))
            .then(pl.lit("rna_rescued"))
            .otherwise(pl.lit("low_confidence"))
            .alias("modality_evidence_dp")
        )

    return df


def compute_multi_allelic_metrics(df: pl.DataFrame) -> pl.DataFrame:
    """Compute per-site multi-allelic metrics.

    Groups by (CHROM, POS) and computes:
        n_alleles_at_site — number of variant rows at this position
        vaf_sum — sum of DNA_VAF_mean across all alleles at this site
        total_alt_dp — sum of DNA_ALT_DP_mean across all alleles
        allele_balance_ratio — ratio of (2nd highest ALT_DP) / (max ALT_DP)
        category_conflict — True if alleles have different FILTER values
        multiallelic_class — classification:
            single: only 1 allele at site
            normalization_artifact: minor_alt_dp / max_alt_dp < 0.03
            noise: vaf_sum < 0.15
            true_multi_allelic: everything else
    """
    if "CHROM" not in df.columns or "POS" not in df.columns:
        return df

    # Compute per-position group metrics
    group_exprs = [pl.len().alias("n_alleles_at_site")]
    if "DNA_VAF_mean" in df.columns:
        group_exprs.append(pl.col("DNA_VAF_mean").sum().alias("vaf_sum"))
    if "DNA_ALT_DP_mean" in df.columns:
        group_exprs.append(pl.col("DNA_ALT_DP_mean").sum().alias("total_alt_dp"))
    if "FILTER" in df.columns:
        group_exprs.append(
            (pl.col("FILTER").n_unique() > 1).alias("category_conflict")
        )

    # Compute allele_balance_ratio = 2nd_max_alt_dp / max_alt_dp per site.
    # Uses two-pass max aggregation to avoid polars list.get() edge cases.
    if "DNA_ALT_DP_mean" in df.columns:
        non_null = df.filter(pl.col("DNA_ALT_DP_mean").is_not_null())
        # Pass 1: max ALT_DP per position
        max_dp = non_null.group_by(["CHROM", "POS"]).agg(
            pl.col("DNA_ALT_DP_mean").max().alias("alt_dp_max")
        )
        # Pass 2: 2nd max = max of values strictly below the max
        with_max = non_null.join(max_dp, on=["CHROM", "POS"])
        second_max = (
            with_max
            .filter(pl.col("DNA_ALT_DP_mean") < pl.col("alt_dp_max"))
            .group_by(["CHROM", "POS"])
            .agg(pl.col("DNA_ALT_DP_mean").max().alias("alt_dp_2nd_max"))
        )
        ranked = max_dp.join(second_max, on=["CHROM", "POS"], how="left")
        ranked = ranked.with_columns(
            pl.when(
                pl.col("alt_dp_max").is_not_null()
                & (pl.col("alt_dp_max") > 0)
                & pl.col("alt_dp_2nd_max").is_not_null()
            )
            .then(pl.col("alt_dp_2nd_max") / pl.col("alt_dp_max"))
            .otherwise(None)
            .alias("allele_balance_ratio")
        )
        ranked = ranked.select(["CHROM", "POS", "allele_balance_ratio"])
    else:
        ranked = df.select(["CHROM", "POS"]).unique().with_columns(
            pl.lit(None).alias("allele_balance_ratio")
        )

    # Join back group metrics
    site_metrics = df.group_by(["CHROM", "POS"]).agg(group_exprs)
    site_metrics = site_metrics.join(ranked, on=["CHROM", "POS"], how="left")

    # Classify multi-allelic sites
    site_metrics = site_metrics.with_columns(
        pl.when(pl.col("n_alleles_at_site") == 1)
        .then(pl.lit("single"))
        .when(
            (pl.col("n_alleles_at_site") > 1)
            & pl.col("allele_balance_ratio").is_not_null()
            & (pl.col("allele_balance_ratio") < 0.03)
        )
        .then(pl.lit("normalization_artifact"))
        .when(
            (pl.col("n_alleles_at_site") > 1)
            & (pl.col("vaf_sum") < 0.15)
        )
        .then(pl.lit("noise"))
        .when(pl.col("n_alleles_at_site") > 1)
        .then(pl.lit("true_multi_allelic"))
        .otherwise(pl.lit("single"))
        .alias("multiallelic_class")
    )

    df = df.join(site_metrics, on=["CHROM", "POS"], how="left")
    return df


def compute_biological_flags(df: pl.DataFrame) -> pl.DataFrame:
    """Add biological flag columns for quality control and filtering.

    Adds the following boolean columns:
        flag_vaf_overflow — True if vaf_sum > 1.1 across alleles at a site
        flag_multi_allelic_heterogeneity — True if multiallelic_class is
            "true_multi_allelic" or "normalization_artifact"
        flag_category_conflict — True if same position has conflicting FILTER values
        flag_germline_low_vaf — FILTER=Germline AND DNA_VAF_mean < 0.10 AND DNA_DP_mean >= 10
        flag_somatic_high_vaf — FILTER=Somatic AND DNA_VAF_mean > 0.60
        flag_reference_with_signal — FILTER=Reference AND DNA_VAF_mean > 0.05
        flag_rna_rescued — DNA_VAF < 0.05 AND N_DNA_CALLERS_SUPPORT <= 1
                           AND N_RNA_CALLERS_SUPPORT >= 2 AND RNA_DP_mean >= 10
    """
    # flag_vaf_overflow
    if "vaf_sum" in df.columns:
        df = df.with_columns(
            (pl.col("vaf_sum") > 1.1).alias("flag_vaf_overflow")
        )

    # flag_multi_allelic_heterogeneity
    if "multiallelic_class" in df.columns:
        df = df.with_columns(
            pl.col("multiallelic_class").is_in(["true_multi_allelic", "normalization_artifact"])
            .alias("flag_multi_allelic_heterogeneity")
        )

    # flag_category_conflict
    if "category_conflict" in df.columns:
        df = df.with_columns(
            pl.col("category_conflict").alias("flag_category_conflict")
        )

    # flag_germline_low_vaf
    has_filter = "FILTER" in df.columns
    has_vaf = "DNA_VAF_mean" in df.columns
    has_dp = "DNA_DP_mean" in df.columns
    if has_filter and has_vaf and has_dp:
        df = df.with_columns(
            ((pl.col("FILTER") == "Germline")
             & (pl.col("DNA_VAF_mean") < 0.10)
             & (pl.col("DNA_DP_mean") >= 10))
            .alias("flag_germline_low_vaf")
        )
    elif "flag_germline_low_vaf" not in df.columns:
        df = df.with_columns(pl.lit(False).alias("flag_germline_low_vaf"))

    # flag_somatic_high_vaf
    if has_filter and has_vaf:
        df = df.with_columns(
            ((pl.col("FILTER") == "Somatic")
             & (pl.col("DNA_VAF_mean") > 0.60))
            .alias("flag_somatic_high_vaf")
        )
    elif "flag_somatic_high_vaf" not in df.columns:
        df = df.with_columns(pl.lit(False).alias("flag_somatic_high_vaf"))

    # flag_reference_with_signal
    if has_filter and has_vaf:
        df = df.with_columns(
            ((pl.col("FILTER") == "Reference")
             & (pl.col("DNA_VAF_mean") > 0.05))
            .alias("flag_reference_with_signal")
        )
    elif "flag_reference_with_signal" not in df.columns:
        df = df.with_columns(pl.lit(False).alias("flag_reference_with_signal"))

    # flag_rna_rescued
    has_dna_callers = "N_DNA_CALLERS_SUPPORT" in df.columns
    has_rna_callers = "N_RNA_CALLERS_SUPPORT" in df.columns
    has_rna_dp = "RNA_DP_mean" in df.columns
    if has_vaf and has_dna_callers and has_rna_callers and has_rna_dp:
        df = df.with_columns(
            ((pl.col("DNA_VAF_mean") < 0.05)
             & (pl.col("N_DNA_CALLERS_SUPPORT") <= 1)
             & (pl.col("N_RNA_CALLERS_SUPPORT") >= 2)
             & (pl.col("RNA_DP_mean") >= 10))
            .alias("flag_rna_rescued")
        )
    elif "flag_rna_rescued" not in df.columns:
        df = df.with_columns(pl.lit(False).alias("flag_rna_rescued"))

    return df


def compute_all_per_variant(df: pl.DataFrame) -> pl.DataFrame:
    """Run all per-variant computations: VAF + mean columns + modality + flags."""
    df = compute_vaf_columns(df)
    df = compute_mean_columns(df)
    df = compute_modality_evidence(df)
    df = compute_multi_allelic_metrics(df)
    df = compute_biological_flags(df)
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

    # Cross-modality — prefer modality_evidence_caller, fall back to old columns
    if "modality_evidence_caller" in df.columns:
        for cat in ["cross_modality", "dna_confident", "rna_rescued", "low_confidence"]:
            result[f"n_mod_{cat}"] = df.filter(pl.col("modality_evidence_caller") == cat).height
        result["n_rescued"] = df.filter(
            pl.col("modality_evidence_caller").is_in(["cross_modality", "rna_rescued"])
        ).height
    else:
        if "CROSS_MODALITY" in df.columns:
            result["n_cross_modality"] = df.filter(pl.col("CROSS_MODALITY") == "YES").height
        if "RESCUED" in df.columns:
            result["n_rescued"] = df.filter(pl.col("RESCUED") == "YES").height

    # Biological flag counts
    bio_flags = [
        "flag_vaf_overflow", "flag_multi_allelic_heterogeneity",
        "flag_category_conflict", "flag_germline_low_vaf",
        "flag_somatic_high_vaf", "flag_reference_with_signal",
        "flag_rna_rescued",
    ]
    for flag in bio_flags:
        if flag in df.columns:
            result[f"n_{flag}"] = df.filter(pl.col(flag) == True).height

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
    """Compute GT concordance among 4 callers with GT fields (polars-native).

    Returns counts of variants where 2, 3, or 4 callers have valid GT values,
    plus a count of variants with <2 valid GTs.
    """
    df = _ensure_eager(df)
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
    """Count how many variants have each flag filter set."""
    df = _ensure_eager(df)
    flags = [
        "min_alt_reads", "gnomad", "blacklist", "noncoding",
        "ig_pseudo", "homopolymer", "vc_filter", "not_consensus",
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

    # Cross-modality — prefer modality_evidence_caller, fall back to old columns
    if "modality_evidence_caller" in df.columns:
        for cat in ["cross_modality", "dna_confident", "rna_rescued", "low_confidence"]:
            result[f"n_mod_{cat}"] = df.filter(pl.col("modality_evidence_caller") == cat).height
        result["n_rescued"] = df.filter(
            pl.col("modality_evidence_caller").is_in(["cross_modality", "rna_rescued"])
        ).height
    else:
        if "CROSS_MODALITY" in df.columns:
            result["n_cross_modality"] = df.filter(pl.col("CROSS_MODALITY") == "YES").height
        if "RESCUED" in df.columns:
            result["n_rescued"] = df.filter(pl.col("RESCUED") == "YES").height

    # Biological flag counts
    bio_flags = [
        "flag_vaf_overflow", "flag_multi_allelic_heterogeneity",
        "flag_category_conflict", "flag_germline_low_vaf",
        "flag_somatic_high_vaf", "flag_reference_with_signal",
        "flag_rna_rescued",
    ]
    for flag in bio_flags:
        if flag in df.columns:
            result[f"n_{flag}"] = df.filter(pl.col(flag) == True).height

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


def compute_vaf_threshold_sweep(df) -> pl.DataFrame:
    """VAF threshold sweep: retention % vs threshold per caller.

    For each caller and each classification category, counts how many
    variants are retained at VAF thresholds from 0.05 to 0.50 in 0.05
    increments. Returns a long-form DataFrame for charting.

    This answers: "If I require VAF ≥ X, what fraction of variants survive?"
    """
    df = _ensure_eager(df)
    callers = ["DNA_mutect2", "DNA_deepsomatic", "DNA_strelka",
               "RNA_mutect2", "RNA_deepsomatic", "RNA_strelka"]
    thresholds = [0.005, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
    classifications = ["Somatic", "Germline", "Reference", "Artifact", "RNAedit", "NoConsensus"]

    has_filter = "FILTER" in df.columns
    rows = []

    for caller in callers:
        vaf_col = f"{caller}_VAF"
        if vaf_col not in df.columns:
            continue
        total = df[vaf_col].drop_nulls().len()
        if total == 0:
            continue

        for thr in thresholds:
            n_retained = df.filter(pl.col(vaf_col) >= thr).height
            row = {"caller": caller, "threshold": thr,
                   "n_retained": n_retained, "n_total": total,
                   "pct_retained": round(n_retained / total * 100, 2)}
            rows.append(row)

            # Per-classification breakdown (if FILTER column present)
            if has_filter:
                for cat in classifications:
                    cat_total = df.filter(pl.col("FILTER") == cat)[vaf_col].drop_nulls().len()
                    if cat_total == 0:
                        continue
                    n = df.filter((pl.col("FILTER") == cat) & (pl.col(vaf_col) >= thr)).height
                    rows.append({"caller": caller, "threshold": thr,
                                 "classification": cat,
                                 "n_retained": n, "n_total": cat_total,
                                 "pct_retained": round(n / cat_total * 100, 2) if cat_total > 0 else 0})

    if not rows:
        return pl.DataFrame()
    return pl.DataFrame(rows)


def compute_dp_threshold_sweep(df) -> pl.DataFrame:
    """DP threshold sweep: retention % vs depth threshold per caller.

    For each caller and DP metric (total DP, REF_DP, ALT_DP), counts how many
    variants are retained at each depth threshold. Returns a long-form DataFrame
    with columns: caller, metric, threshold, n_retained, n_total, pct_retained.

    Thresholds:
      Total DP: [1, 2, 5, 10, 20, 50, 100, 200]
      REF/ALT DP: [0, 1, 2, 5, 10, 20, 50]

    This answers: "If I require DP >= X, what fraction of variants survive?"
    """
    df = _ensure_eager(df)
    callers = ["DNA_mutect2", "DNA_deepsomatic", "DNA_strelka",
               "RNA_mutect2", "RNA_deepsomatic", "RNA_strelka"]
    dp_thresholds = [1, 2, 5, 10, 20, 50, 100, 200]
    refalt_thresholds = [0, 1, 2, 5, 10, 20, 50]
    classifications = ["Somatic", "Germline", "Reference", "Artifact", "RNAedit", "NoConsensus"]
    has_filter = "FILTER" in df.columns
    rows = []

    # Per-caller total DP sweep
    for caller in callers:
        dp_col = f"{caller}_DP"
        if dp_col not in df.columns:
            continue
        total = df[dp_col].drop_nulls().len()
        if total == 0:
            continue
        for thr in dp_thresholds:
            n = df.filter(pl.col(dp_col) >= thr).height
            rows.append({"caller": caller, "metric": "DP", "threshold": thr,
                         "n_retained": n, "n_total": total,
                         "pct_retained": round(n / total * 100, 2)})

            # Per-classification breakdown (if FILTER column present)
            if has_filter:
                for cat in classifications:
                    cat_total = df.filter(pl.col("FILTER") == cat)[dp_col].drop_nulls().len()
                    if cat_total == 0:
                        continue
                    n_cat = df.filter((pl.col("FILTER") == cat) & (pl.col(dp_col) >= thr)).height
                    rows.append({"caller": caller, "metric": "DP", "threshold": thr,
                                 "classification": cat,
                                 "n_retained": n_cat, "n_total": cat_total,
                                 "pct_retained": round(n_cat / cat_total * 100, 2) if cat_total > 0 else 0})

    # BAM pileup DP sweeps (REF_DP, ALT_DP for DT and RT)
    for bt in ["DT", "RT"]:
        for suffix, metric_name in [("DP", "BAM_DP"), ("REF_DP", "BAM_REF_DP"), ("ALT_DP", "BAM_ALT_DP")]:
            col = f"BAM_{bt}_{suffix}"
            if col not in df.columns:
                continue
            total = df[col].drop_nulls().len()
            if total == 0:
                continue
            thresholds = dp_thresholds if suffix == "DP" else refalt_thresholds
            for thr in thresholds:
                n = df.filter(pl.col(col) >= thr).height
                rows.append({"caller": f"BAM_{bt}", "metric": metric_name, "threshold": thr,
                             "n_retained": n, "n_total": total,
                             "pct_retained": round(n / total * 100, 2)})

    if not rows:
        return pl.DataFrame()
    return pl.DataFrame(rows)


def compute_filter_effectiveness_matrix(df) -> pl.DataFrame:
    """Filter effectiveness matrix: FILTER × Classification.

    For each FILTER value (Somatic, Germline, etc.) counts how many
    variants have each flag filter set (min_alt_reads, gnomad, blacklist,
    etc.). Returns a pivoted matrix for heatmap visualization.

    This answers: "Which filters are most effective for each classification?"
    """
    df = _ensure_eager(df)
    flags = ["min_alt_reads", "gnomad", "blacklist", "noncoding",
             "ig_pseudo", "homopolymer", "vc_filter", "not_consensus"]

    if "FILTER" not in df.columns:
        return pl.DataFrame()

    rows = []
    for cat in ["Somatic", "Germline", "Reference", "Artifact", "RNAedit", "NoConsensus"]:
        cat_df = df.filter(pl.col("FILTER") == cat)
        cat_total = cat_df.height
        if cat_total == 0:
            continue
        for flag in flags:
            if flag not in df.columns:
                continue
            n = cat_df.filter(pl.col(flag) == True).height
            rows.append({
                "classification": cat,
                "filter_flag": flag,
                "n_flagged": n,
                "n_total": cat_total,
                "pct_flagged": round(n / cat_total * 100, 2),
            })

    if not rows:
        return pl.DataFrame()
    return pl.DataFrame(rows)


def sample_tier_summary(df: pl.DataFrame | pl.LazyFrame) -> pl.DataFrame:
    """Level 4: Per-sample x per-tier aggregate statistics.

    For each (sample_id, final_tier) pair, compute: variant count, mean VAF/DP/
    REF_DP/ALT_DP, variant type distribution, Ti/Tv, and N_SUPPORT_CALLERS dist.
    If set_number column exists, it is included in the group-by for per-set faceting.
    """
    df = _ensure_eager(df)
    if "sample_id" not in df.columns or "final_tier" not in df.columns:
        return pl.DataFrame()

    group_cols = ["sample_id", "final_tier"]
    if "set_number" in df.columns:
        group_cols.append("set_number")

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
        df.group_by(group_cols)
        .agg(agg_exprs)
        .sort(group_cols)
    )


# ═══════════════════════════════════════════════════════════════════════════════
# ML Threshold Guidance Statistics (Section 7)
# ═══════════════════════════════════════════════════════════════════════════════


def compute_filter_vaf_dp_cross_tab(df) -> pl.DataFrame:
    """FILTER x VAF_bin x DP_bin cross-tabulation for ML filtering guidance.

    Produces a long-form table of (classification, vaf_bin, dp_bin, count)
    rows, optionally partitioned by the ``partition`` column (if present).
    Useful for understanding the joint distribution of FILTER, VAF, and DP
    which informs ML threshold selection.
    """
    df = _ensure_eager(df)
    vaf_col = "DNA_VAF_mean" if "DNA_VAF_mean" in df.columns else None
    dp_col = "DNA_DP_mean" if "DNA_DP_mean" in df.columns else None
    if not vaf_col or not dp_col or "FILTER" not in df.columns:
        return pl.DataFrame()

    vaf_bins = [0, 0.01, 0.05, 0.10, 0.25, 0.50, 1.0]
    dp_bins = [0, 10, 50, 100, 200, 500, float("inf")]
    vaf_labels = ["<0.01", "0.01-0.05", "0.05-0.10", "0.10-0.25", "0.25-0.50", "0.50-1.0"]
    dp_labels = ["<10", "10-50", "50-100", "100-200", "200-500", "500+"]

    rows = []
    has_partition = "partition" in df.columns
    partitions = df["partition"].unique().to_list() if has_partition else [None]

    for part in partitions:
        sub = df.filter(pl.col("partition") == part) if part else df
        for filt in ["Somatic", "Germline", "Reference", "Artifact", "RNAedit", "NoConsensus"]:
            filt_df = sub.filter(pl.col("FILTER") == filt)
            if filt_df.is_empty():
                continue
            for vi in range(len(vaf_bins) - 1):
                for di in range(len(dp_bins) - 1):
                    n = filt_df.filter(
                        (pl.col(vaf_col) >= vaf_bins[vi]) & (pl.col(vaf_col) < vaf_bins[vi + 1])
                        & (pl.col(dp_col) >= dp_bins[di]) & (pl.col(dp_col) < dp_bins[di + 1])
                    ).height
                    if n > 0:
                        row = {"classification": filt, "vaf_bin": vaf_labels[vi],
                               "dp_bin": dp_labels[di], "count": n}
                        if part:
                            row["partition"] = part
                        rows.append(row)
    return pl.DataFrame(rows) if rows else pl.DataFrame()


def compute_low_vaf_rna_support(df) -> pl.DataFrame:
    """Variants with VAF < 0.05 where N_RNA_CALLERS_SUPPORT >= 2.

    Summarises low-VAF variants that have strong RNA caller support,
    grouped by FILTER classification.  Useful for identifying potential
    true somatic variants rescued by RNA evidence.
    """
    df = _ensure_eager(df)
    if "DNA_VAF_mean" not in df.columns or "N_RNA_CALLERS_SUPPORT" not in df.columns:
        return pl.DataFrame()
    low_vaf = df.filter(
        (pl.col("DNA_VAF_mean") < 0.05) & (pl.col("N_RNA_CALLERS_SUPPORT") >= 2)
    )
    if low_vaf.is_empty():
        return pl.DataFrame()

    agg_exprs = [
        pl.len().alias("n_variants"),
        pl.col("DNA_VAF_mean").mean().alias("mean_vaf"),
    ]
    if "DNA_DP_mean" in df.columns:
        agg_exprs.append(pl.col("DNA_DP_mean").mean().alias("mean_dp"))

    return low_vaf.group_by("FILTER").agg(agg_exprs).sort("n_variants", descending=True)


# ═══════════════════════════════════════════════════════════════════════════════
# FP Cross-Tabulation (Section 8)
# ═══════════════════════════════════════════════════════════════════════════════


def compute_fp_cross_tab(df) -> pl.DataFrame:
    """Cross-tabulation: FILTER x N_SUPPORT_CALLERS x count.

    Shows caller-support profile for all variant classifications including
    Somatic (provides TP contrast against non-somatic FP patterns).
    """
    df = _ensure_eager(df)
    if "FILTER" not in df.columns or "N_SUPPORT_CALLERS" not in df.columns:
        return pl.DataFrame()

    filtered = df.filter(pl.col("FILTER").is_not_null())
    if filtered.is_empty():
        return pl.DataFrame()

    return (
        filtered.group_by(["FILTER", "N_SUPPORT_CALLERS"])
        .agg(pl.len().alias("count"))
        .sort(["FILTER", "N_SUPPORT_CALLERS"])
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Somatic Modality Sub-Classification (Section 9)
# ═══════════════════════════════════════════════════════════════════════════════


def compute_somatic_modality(df) -> pl.DataFrame:
    """Derive somatic modality sub-classification from caller_tier.

    Maps caller tiers to modality categories:
      C1 -> MultiModality   (both DNA and RNA callers agree)
      C2, C5 -> DNA_only    (DNA callers only)
      C3, C6 -> RNA_only    (RNA callers only)
      C4, C7 -> Weak        (low caller support)

    Returns per-modality summary with variant counts and mean VAF/DP.
    """
    df = _ensure_eager(df)
    if "FILTER" not in df.columns or "caller_tier" not in df.columns:
        return pl.DataFrame()

    somatic = df.filter(pl.col("FILTER") == "Somatic")
    if somatic.is_empty():
        return pl.DataFrame()

    somatic = somatic.with_columns(
        pl.when(pl.col("caller_tier") == "C1").then(pl.lit("MultiModality"))
        .when(pl.col("caller_tier").is_in(["C2", "C5"])).then(pl.lit("DNA_only"))
        .when(pl.col("caller_tier").is_in(["C3", "C6"])).then(pl.lit("RNA_only"))
        .when(pl.col("caller_tier").is_in(["C4", "C7"])).then(pl.lit("Weak"))
        .otherwise(pl.lit("Unknown"))
        .alias("somatic_modality")
    )

    agg_exprs = [pl.len().alias("n_variants")]
    if "DNA_VAF_mean" in somatic.columns:
        agg_exprs.append(pl.col("DNA_VAF_mean").mean().alias("mean_dna_vaf"))
    if "RNA_VAF_mean" in somatic.columns:
        agg_exprs.append(pl.col("RNA_VAF_mean").mean().alias("mean_rna_vaf"))
    if "DNA_DP_mean" in somatic.columns:
        agg_exprs.append(pl.col("DNA_DP_mean").mean().alias("mean_dna_dp"))

    return somatic.group_by("somatic_modality").agg(agg_exprs).sort("n_variants", descending=True)


# ═══════════════════════════════════════════════════════════════════════════════
# Rescue Analytics (Section 9)
# ═══════════════════════════════════════════════════════════════════════════════


def _rescue_flag(df: pl.DataFrame) -> pl.DataFrame:
    """Normalize rescue status to YES/NO using modality_evidence_caller.

    Prefers modality_evidence_caller column; falls back to RESCUED column
    for compatibility with older parquet files.
    """
    if "modality_evidence_caller" in df.columns:
        return df.with_columns(
            pl.when(
                pl.col("modality_evidence_caller").is_in(["cross_modality", "rna_rescued"])
            ).then(pl.lit("YES"))
            .otherwise(pl.lit("NO"))
            .alias("RESCUED")
        )
    if "RESCUED" not in df.columns:
        return df.with_columns(pl.lit("NO").alias("RESCUED"))
    return df.with_columns(
        pl.when(pl.col("RESCUED") == "YES").then(pl.lit("YES"))
        .otherwise(pl.lit("NO"))
        .alias("RESCUED")
    )


def _evidence_col(df: pl.DataFrame) -> str:
    """Return the primary evidence grouping column name.

    Prefers modality_evidence_caller for new parquet files;
    falls back to CROSS_MODALITY/RESCUED for backward compat.
    """
    if "modality_evidence_caller" in df.columns:
        return "modality_evidence_caller"
    return "RESCUED"


def compute_rescue_breakdown(df, group_col: str = "set_number") -> pl.DataFrame:
    """Per-group modality evidence category counts and proportions."""
    df = _ensure_eager(df)
    ev_col = _evidence_col(df)
    if ev_col not in df.columns:
        return pl.DataFrame()
    group_cols = [group_col, ev_col] if group_col in df.columns else [ev_col]
    result = df.group_by(group_cols).agg(pl.len().alias("count"))
    if group_col in df.columns:
        totals = result.group_by(group_col).agg(pl.col("count").sum().alias("total"))
        result = result.join(totals, on=group_col, how="left")
    else:
        total = result["count"].sum()
        result = result.with_columns(pl.lit(total).alias("total"))
    result = result.with_columns(
        (pl.col("count") / pl.col("total") * 100).round(2).alias("pct")
    )
    return result.sort(group_cols if group_col in df.columns else [ev_col])


def compute_rescue_by_filter(df) -> pl.DataFrame:
    """Modality evidence × FILTER cross-tabulation with counts and proportions."""
    df = _ensure_eager(df)
    ev_col = _evidence_col(df)
    if ev_col not in df.columns or "FILTER" not in df.columns:
        return pl.DataFrame()
    result = df.group_by(["FILTER", ev_col]).agg(pl.len().alias("count"))
    totals = result.group_by("FILTER").agg(pl.col("count").sum().alias("total_per_filter"))
    result = result.join(totals, on="FILTER", how="left")
    result = result.with_columns(
        (pl.col("count") / pl.col("total_per_filter") * 100).round(2).alias("pct")
    )
    return result.sort(["FILTER", ev_col])


def compute_rescue_cross_tab(df) -> pl.DataFrame:
    """Three-way cross-tabulation: modality evidence × FILTER × set_number."""
    df = _ensure_eager(df)
    ev_col = _evidence_col(df)
    if ev_col not in df.columns or "FILTER" not in df.columns:
        return pl.DataFrame()
    group_cols = [ev_col, "FILTER"]
    if "set_number" in df.columns:
        group_cols.append("set_number")
    return df.group_by(group_cols).agg(pl.len().alias("count")).sort(group_cols)


def compute_rescue_vaf_dp(df) -> pl.DataFrame:
    """VAF/DP distribution statistics by modality evidence (mean, median, Q1, Q3)."""
    df = _ensure_eager(df)
    ev_col = _evidence_col(df)
    if ev_col not in df.columns:
        return pl.DataFrame()
    agg_exprs = [pl.len().alias("n_variants")]
    for col in ["DNA_VAF_mean", "RNA_VAF_mean", "DNA_DP_mean", "RNA_DP_mean"]:
        if col in df.columns:
            agg_exprs.extend([
                pl.col(col).mean().alias(f"{col}_mean"),
                pl.col(col).median().alias(f"{col}_median"),
                pl.col(col).quantile(0.25).alias(f"{col}_q1"),
                pl.col(col).quantile(0.75).alias(f"{col}_q3"),
            ])
    return df.group_by(ev_col).agg(agg_exprs).sort(ev_col)


def sample_rescue_summary(df) -> pl.DataFrame:
    """Per-sample modality evidence counts with FILTER breakdown."""
    df = _ensure_eager(df)
    ev_col = _evidence_col(df)
    if ev_col not in df.columns or "sample_id" not in df.columns:
        return pl.DataFrame()
    group_cols = ["sample_id", ev_col]
    if "set_number" in df.columns:
        group_cols.append("set_number")
    if "FILTER" in df.columns:
        group_cols.append("FILTER")
    return df.group_by(group_cols).agg(pl.len().alias("count")).sort(group_cols)


def compute_rescue_by_tier(df) -> pl.DataFrame:
    """Modality evidence count and rate per final_tier (CxDy)."""
    df = _ensure_eager(df)
    ev_col = _evidence_col(df)
    if ev_col not in df.columns or "final_tier" not in df.columns:
        return pl.DataFrame()
    result = df.group_by(["final_tier", ev_col]).agg(pl.len().alias("count"))
    totals = result.group_by("final_tier").agg(pl.col("count").sum().alias("total"))
    result = result.join(totals, on="final_tier", how="left")
    result = result.with_columns(
        (pl.col("count") / pl.col("total") * 100).round(2).alias("pct")
    )
    return result.sort(["final_tier", ev_col])


def compute_rescue_by_caller_support(df) -> pl.DataFrame:
    """N_SUPPORT_CALLERS distribution by modality evidence."""
    df = _ensure_eager(df)
    ev_col = _evidence_col(df)
    if ev_col not in df.columns or "N_SUPPORT_CALLERS" not in df.columns:
        return pl.DataFrame()
    return df.group_by(["N_SUPPORT_CALLERS", ev_col]).agg(
        pl.len().alias("count")
    ).sort(["N_SUPPORT_CALLERS", ev_col])
