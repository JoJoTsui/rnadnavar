"""Cross-validate rescue VCF pre-computed INFO means against caller ground truth.

Compares values computed from individual caller VCFs (ground truth) against
the pre-computed mean values in the rescue VCF INFO fields. Flags discrepancies
beyond a configurable tolerance.
"""

from typing import Any

import polars as pl

# ── Validation metric pairs: (computed_col, rescue_info_col, label) ───────

VALIDATION_METRICS = [
    ("DNA_DP_mean", "DP_DNA_MEAN", "DNA DP mean"),
    ("RNA_DP_mean", "DP_RNA_MEAN", "RNA DP mean"),
    ("DNA_VAF_mean", "VAF_DNA_MEAN", "DNA VAF mean"),
    ("RNA_VAF_mean", "VAF_RNA_MEAN", "RNA VAF mean"),
    # Overall means (computed from all 6 callers)
    ("DNA_DP_mean", "DP_MEAN", "Overall DP mean (DNA only)"),
    ("DNA_VAF_mean", "VAF_MEAN", "Overall VAF mean (DNA only)"),
]


def validate_single_metric(
    df: pl.DataFrame,
    computed_col: str,
    rescue_col: str,
    tolerance: float = 0.01,
) -> dict[str, Any]:
    """Validate one metric for a single sample.

    Args:
        df: DataFrame with both computed and rescue columns.
        computed_col: Column name of the ground-truth computed value.
        rescue_col: Column name of the rescue VCF INFO value.
        tolerance: Maximum absolute difference to consider a match.

    Returns:
        Dict with validation statistics.
    """
    if computed_col not in df.columns or rescue_col not in df.columns:
        return {
            "metric": rescue_col,
            "n_total": 0,
            "n_match": 0,
            "n_mismatch": 0,
            "n_missing": 0,
            "mismatch_pct": 0.0,
            "max_abs_diff": None,
            "mean_abs_diff": None,
            "status": "missing_columns",
        }

    # Compute absolute difference, filter out rows where either is null
    diff_df = df.select([
        pl.col(computed_col).alias("computed"),
        pl.col(rescue_col).alias("rescue"),
    ]).with_columns(
        (pl.col("computed") - pl.col("rescue")).abs().alias("abs_diff")
    )

    n_total = diff_df.height
    valid = diff_df.filter(
        pl.col("computed").is_not_null() & pl.col("rescue").is_not_null()
    )
    n_valid = valid.height
    n_missing = n_total - n_valid

    if n_valid == 0:
        return {
            "metric": rescue_col,
            "n_total": n_total,
            "n_match": 0,
            "n_mismatch": 0,
            "n_missing": n_missing,
            "mismatch_pct": 0.0,
            "max_abs_diff": None,
            "mean_abs_diff": None,
            "status": "no_valid_pairs",
        }

    mismatches = valid.filter(pl.col("abs_diff") > tolerance)
    n_mismatch = mismatches.height
    n_match = n_valid - n_mismatch

    max_diff = valid["abs_diff"].max()
    mean_diff = valid["abs_diff"].mean()

    return {
        "metric": rescue_col,
        "n_total": n_total,
        "n_match": n_match,
        "n_mismatch": n_mismatch,
        "n_missing": n_missing,
        "mismatch_pct": n_mismatch / n_valid * 100 if n_valid > 0 else 0.0,
        "max_abs_diff": max_diff,
        "mean_abs_diff": mean_diff,
        "status": "ok",
    }


def validate_sample(
    df: pl.DataFrame,
    sample_id: str,
    tolerance: float = 0.01,
) -> list[dict[str, Any]]:
    """Run all validation metrics for a single sample.

    Returns a list of result dicts, one per metric.
    """
    results = []
    for computed_col, rescue_col, label in VALIDATION_METRICS:
        result = validate_single_metric(df, computed_col, rescue_col, tolerance)
        result["sample_id"] = sample_id
        result["label"] = label
        results.append(result)
    return results


def validate_all_samples(
    samples_data: dict[str, pl.DataFrame],
    tolerance: float = 0.01,
) -> pl.DataFrame:
    """Run validation on all samples and return a combined report DataFrame.

    Args:
        samples_data: Dict mapping sample_id -> per-variant DataFrame.
        tolerance: Maximum absolute difference for a match.

    Returns:
        polars DataFrame with one row per metric per sample.
    """
    all_results = []
    for sample_id, df in samples_data.items():
        results = validate_sample(df, sample_id, tolerance)
        all_results.extend(results)

    if not all_results:
        return pl.DataFrame()

    return pl.DataFrame(all_results)


def validation_summary(report: pl.DataFrame) -> pl.DataFrame:
    """Generate a per-metric summary across all samples."""
    if report.is_empty():
        return pl.DataFrame()

    return (
        report.group_by("metric")
        .agg([
            pl.col("n_total").sum().alias("total_variants"),
            pl.col("n_match").sum().alias("total_match"),
            pl.col("n_mismatch").sum().alias("total_mismatch"),
            pl.col("n_missing").sum().alias("total_missing"),
            (pl.col("n_mismatch").sum() * 100.0 / pl.col("n_total").sum()).alias("overall_mismatch_pct"),
            pl.col("max_abs_diff").max().alias("worst_max_abs_diff"),
            pl.col("mean_abs_diff").mean().alias("avg_mean_abs_diff"),
        ])
        .sort("overall_mismatch_pct", descending=True)
    )
