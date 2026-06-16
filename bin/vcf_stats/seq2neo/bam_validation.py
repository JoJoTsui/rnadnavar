"""Cross-validate BAM pileup metrics against caller VCF FORMAT fields.

Computes per-caller per-sample validation statistics: DP comparison,
VAF comparison (Mutect2 + DeepSomatic only), and strand bias comparison.
"""

from typing import Any

import numpy as np
import polars as pl

from .manifest_loader import CALLER_CONFIGS

CALLERS_WITH_GT = ["DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic", "RNA_deepsomatic"]


def validate_bam_vs_caller(
    df: pl.DataFrame,
    sample_id: str,
) -> dict[str, Any]:
    """Validate BAM pileup DP/VAF against caller VCF data for one sample.

    Args:
        df: DataFrame with BAM pileup columns (BAM_DP_DT, etc.) and
            caller columns ({caller}_DP, {caller}_VAF, etc.)
        sample_id: Sample identifier.

    Returns:
        Dict with per-caller validation metrics or empty if no BAM data.
    """
    # Detect available BAM types
    bam_types = []
    for bt in ["DN", "DT", "RT"]:
        if f"BAM_{bt}_DP" in df.columns:
            bam_types.append(bt)

    if not bam_types:
        return {"sample_id": sample_id, "has_bam_data": False}

    results = {"sample_id": sample_id, "has_bam_data": True, "bam_types": ",".join(bam_types)}

    for caller in CALLER_CONFIGS:
        dp_col = f"{caller}_DP"
        if dp_col not in df.columns:
            continue

        # Match BAM type to caller: DNA callers → DT, RNA callers → RT
        if caller.startswith("DNA"):
            bam_dp_col = "BAM_DT_DP"
            bam_alt_col = "BAM_DT_ALT_DP"
            bam_ref_col = "BAM_DT_REF_DP"
        else:
            bam_dp_col = "BAM_RT_DP"
            bam_alt_col = "BAM_RT_ALT_DP"
            bam_ref_col = "BAM_RT_REF_DP"

        if bam_dp_col not in df.columns:
            continue

        # DP comparison: BAM DP vs caller DP
        valid = df.filter(
            pl.col(dp_col).is_not_null() & pl.col(bam_dp_col).is_not_null()
            & (pl.col(dp_col) > 0) & (pl.col(bam_dp_col) > 0)
        )

        if valid.height > 0:
            dp_diff = (valid[dp_col].cast(pl.Float64) - valid[bam_dp_col].cast(pl.Float64)).abs()
            results[f"{caller}_dp_corr"] = _safe_corr(valid[dp_col], valid[bam_dp_col])
            results[f"{caller}_dp_mad"] = dp_diff.mean()
            results[f"{caller}_dp_mismatch_pct"] = (dp_diff > valid[dp_col] * 0.5).sum() / valid.height * 100

        # VAF comparison (Mutect2 + DeepSomatic only — they have VAF in FORMAT)
        if caller in CALLERS_WITH_GT:
            vaf_col = f"{caller}_VAF"
            if vaf_col in df.columns and bam_alt_col in df.columns and bam_dp_col in df.columns:
                # BAM VAF = ALT / DP
                bam_vaf = df.select([
                    pl.when(pl.col(bam_dp_col).is_not_null() & (pl.col(bam_dp_col) > 0))
                    .then(pl.col(bam_alt_col).cast(pl.Float64) / pl.col(bam_dp_col).cast(pl.Float64))
                    .otherwise(None)
                    .alias("_bam_vaf")
                ])["_bam_vaf"]

                vaf_valid = df.select([
                    pl.col(vaf_col).alias("_cv"),
                    bam_vaf.alias("_bv"),
                ]).filter(
                    pl.col("_cv").is_not_null() & pl.col("_bv").is_not_null()
                )

                if vaf_valid.height > 0:
                    results[f"{caller}_vaf_corr"] = _safe_corr(vaf_valid["_cv"], vaf_valid["_bv"])
                    vaf_diff = (vaf_valid["_cv"] - vaf_valid["_bv"]).abs()
                    results[f"{caller}_vaf_mad"] = vaf_diff.mean()

        # Strand bias comparison (Mutect2 only — has SB in FORMAT)
        if "mutect2" in caller.lower():
            sb_col = f"{caller}_SB"
            bam_type_strand = "DT" if caller.startswith("DNA") else "RT"
            if sb_col in df.columns and bam_type_has_strand(df, bam_type_strand):
                # Simple check: BAM has F1R2_alt/F2R1_alt, caller has SB
                results[f"{caller}_has_sb_bam"] = True

    return results


def bam_type_has_strand(df: pl.DataFrame, bam_type: str) -> bool:
    """Check if BAM strand columns exist for the given BAM type (DT/RT)."""
    f1r2_col = f"BAM_{bam_type}_F1R2_alt"
    return f1r2_col in df.columns


def validate_bam_all_samples(
    samples_data: dict[str, pl.DataFrame],
) -> pl.DataFrame:
    """Run BAM validation across all samples.

    Returns a DataFrame with one row per sample containing per-caller validation metrics.
    """
    all_results = []
    for sample_id, df in samples_data.items():
        result = validate_bam_vs_caller(df, sample_id)
        all_results.append(result)

    if not all_results:
        return pl.DataFrame()

    return pl.DataFrame(all_results)


def _safe_corr(a: pl.Series, b: pl.Series) -> float | None:
    """Compute Pearson correlation safely, returning None on failure."""
    try:
        a_arr = a.to_numpy()
        b_arr = b.to_numpy()
        mask = ~(np.isnan(a_arr) | np.isnan(b_arr))
        if mask.sum() < 3:
            return None
        return float(np.corrcoef(a_arr[mask], b_arr[mask])[0, 1])
    except Exception:
        return None
