"""Bridge module: integrate CxDy variant tiering into seq2neo statistics.

Wraps the existing TieringEngine (bin/vcf_stats/tiering_engine.py) to compute
CxDy tiers from already-parsed rescue VCF INFO fields without re-reading VCFs.

The tiering system:
  - C1-C7: Caller support tiers based on concordant DNA/RNA caller counts
  - D0-D1: Database evidence tiers (gnomAD, COSMIC, REDIportal, DARNED)
  - Final tier: CxDy format (e.g., C1D1 = strong both modalities + DB evidence)
"""

import re
from typing import Any

import polars as pl

# Import TieringEngine from sibling tiering_engine.py
from ..tiering_engine import TieringEngine

# ── Column conversion helpers (for Rust FFI) ──────────────────────────────

def _col_to_opt_float(df: pl.DataFrame, col: str) -> list:
    """Convert DataFrame column to list of Optional[float] for Rust."""
    if col not in df.columns:
        return [None] * len(df)
    vals = []
    for v in df[col].to_list():
        if v is None or v == ".":
            vals.append(None)
        else:
            try:
                vals.append(float(v))
            except (ValueError, TypeError):
                vals.append(None)
    return vals

def _col_to_opt_int(df: pl.DataFrame, col: str) -> list:
    """Convert DataFrame column to list of Optional[int] for Rust."""
    if col not in df.columns:
        return [None] * len(df)
    vals = []
    for v in df[col].to_list():
        if v is None or v == ".":
            vals.append(None)
        else:
            try:
                vals.append(int(v))
            except (ValueError, TypeError):
                vals.append(None)
    return vals

def _col_to_opt_str(df: pl.DataFrame, col: str) -> list:
    """Convert DataFrame column to list of Optional[str] for Rust."""
    if col not in df.columns:
        return [None] * len(df)
    vals = []
    for v in df[col].to_list():
        if v is None or v == ".":
            vals.append(None)
        else:
            vals.append(str(v))
    return vals

# ── Tier order constants (from tier_config.py) ───────────────────────────────
CALLER_TIER_ORDER = ["C1", "C2", "C3", "C4", "C5", "C6", "C7"]
FINAL_TIER_ORDER = [f"{c}{d}" for c in CALLER_TIER_ORDER for d in ["D1", "D0"]]

# Single TieringEngine instance (reusable, stateless)
_engine = None


def _get_engine() -> TieringEngine:
    """Get or create the TieringEngine singleton."""
    global _engine
    if _engine is None:
        _engine = TieringEngine()
    return _engine


def _parse_filters_normalized(filters_str: str | None) -> dict[str, str | None]:
    """Parse FILTERS_NORMALIZED string into FILTER_NORMALIZED_<Caller>_<Modality> dict.

    Format: "DNA_strelka:Somatic|RNA_mutect2:Germline|..." (pipe-separated)

    Maps caller names to their normalized form expected by TieringEngine:
      Strelka, DeepSomatic, Mutect2 × DNA_TUMOR, RNA_TUMOR
    """
    result: dict[str, str | None] = {}
    if not filters_str:
        return result

    # Split by pipe or semicolon
    entries = filters_str.split("|") if "|" in filters_str else filters_str.split(";")

    for entry in entries:
        entry = entry.strip()
        if not entry or ":" not in entry:
            continue

        parts = entry.split(":", 1)
        if len(parts) != 2:
            continue

        caller_spec, category = parts
        caller_spec = caller_spec.strip()
        category = category.strip()

        # Parse "DNA_strelka" or "RNA_mutect2" format
        match = re.match(r"^(DNA|RNA)_(\w+)$", caller_spec, re.IGNORECASE)
        if not match:
            continue

        modality_raw = match.group(1).upper()
        caller_raw = match.group(2)

        # Normalize caller name (case-insensitive)
        caller_map = {"strelka": "Strelka", "deepsomatic": "DeepSomatic", "mutect2": "Mutect2"}
        caller_normalized = caller_map.get(caller_raw.lower())
        if not caller_normalized:
            continue

        modality = f"{modality_raw}_TUMOR"
        field_name = f"FILTER_NORMALIZED_{caller_normalized}_{modality}"
        result[field_name] = category

    return result


def compute_tiers_for_dataframe(df: pl.DataFrame) -> pl.DataFrame:
    """Compute CxDy tiers for all variants in a DataFrame.

    Uses Rust stats_core.compute_tiers() when available (fast, GIL-released).
    Falls back to Python TieringEngine when Rust is unavailable.

    Args:
        df: polars DataFrame with rescue VCF columns: FILTER, FILTERS_NORMALIZED,
            GNOMAD_AF, COSMIC_CNT, REDI_EVIDENCE (optional).

    Returns:
        DataFrame with added columns: caller_tier, database_tier, final_tier,
        tier_quality, dna_caller_count, rna_caller_count.
    """
    # Try Rust path first
    try:
        import stats_core
        if hasattr(stats_core, 'compute_tiers'):
            filters = df["FILTER"].to_list()
            filters_norm = df["FILTERS_NORMALIZED"].to_list() if "FILTERS_NORMALIZED" in df.columns else [""] * len(filters)
            gnomad_af = _col_to_opt_float(df, "GNOMAD_AF")
            cosmic_cnt = _col_to_opt_int(df, "COSMIC_CNT")
            redi_ev = _col_to_opt_str(df, "REDI_EVIDENCE")
            dna_sup = _col_to_opt_int(df, "N_DNA_CALLERS_SUPPORT")
            rna_sup = _col_to_opt_int(df, "N_RNA_CALLERS_SUPPORT")

            result = stats_core.compute_tiers(filters, filters_norm, gnomad_af, cosmic_cnt, redi_ev, dna_sup, rna_sup)
            tier_df = pl.DataFrame(result)
            return df.hstack(tier_df)
    except Exception as e:
        print(f"  [TIERING] Rust tiering failed ({e}), falling back to Python")

    # Python fallback
    engine = _get_engine()
    rows = df.to_dicts()
    tier_results: list[dict[str, Any]] = []

    for row in rows:
        final_filter = row.get("FILTER", "PASS")
        if final_filter is None or final_filter == ".":
            final_filter = "PASS"

        # Build FILTER_NORMALIZED_* dict from the legacy string
        filters_str = row.get("FILTERS_NORMALIZED")
        filter_normalized_fields = _parse_filters_normalized(filters_str)

        # Also check for individual FILTER_NORMALIZED_* columns in the row
        for key, value in row.items():
            if key.startswith("FILTER_NORMALIZED_") and value is not None:
                filter_normalized_fields[key] = str(value)

        # Build info_dict for database checking
        info_dict: dict[str, Any] = {}
        gnomad_af = row.get("GNOMAD_AF")
        if gnomad_af is not None:
            try: info_dict["gnomAD_AF"] = float(gnomad_af)
            except (ValueError, TypeError): pass
        cosmic_cnt = row.get("COSMIC_CNT")
        if cosmic_cnt is not None:
            try: info_dict["COSMIC_CNT"] = int(cosmic_cnt)
            except (ValueError, TypeError): pass

        # REDIportal: consider REDI_EVIDENCE != "NONE" as database support
        redi_evidence = row.get("REDI_EVIDENCE")
        if redi_evidence is not None and redi_evidence in ("HIGH", "MEDIUM", "LOW"):
            info_dict["REDIportal"] = True

        try:
            tier_info = engine.compute_tier(
                final_filter=final_filter,
                filter_normalized_fields=filter_normalized_fields,
                info_dict=info_dict if info_dict else None,
            )
        except (ValueError, RuntimeError):
            # Fallback: compute tier from simple caller counts
            dna_count = row.get("N_DNA_CALLERS_SUPPORT")
            rna_count = row.get("N_RNA_CALLERS_SUPPORT")
            try: dna_count = int(dna_count) if dna_count is not None else 0
            except (ValueError, TypeError): dna_count = 0
            try: rna_count = int(rna_count) if rna_count is not None else 0
            except (ValueError, TypeError): rna_count = 0
            has_db = bool(info_dict)
            tier_str = engine.compute_tier_simple(
                dna_caller_count=dna_count,
                rna_caller_count=int(rna_count),
                has_database_support=has_db,
            )
            caller_tier = tier_str[:2]
            database_tier = tier_str[2:]
            tier_info = {
                "final_tier": tier_str,
                "caller_tier": caller_tier,
                "database_tier": database_tier,
                "dna_caller_count": int(dna_count),
                "rna_caller_count": int(rna_count),
                "tier_quality": 0,
            }

        tier_results.append({
            "caller_tier": tier_info["caller_tier"],
            "database_tier": tier_info["database_tier"],
            "final_tier": tier_info["final_tier"],
            "tier_quality": tier_info.get("tier_quality", 0),
            "dna_caller_count": tier_info.get("dna_caller_count", 0),
            "rna_caller_count": tier_info.get("rna_caller_count", 0),
        })

    # Add tier columns to DataFrame
    tier_df = pl.DataFrame(tier_results)
    return df.hstack(tier_df)


def tier_summary(df: pl.DataFrame) -> pl.DataFrame:
    """Compute per-tier aggregate statistics.

    Args:
        df: DataFrame with final_tier column and computed statistics columns.

    Returns:
        DataFrame with one row per CxDy tier containing count, mean VAF, mean DP,
        mean REF_DP, mean ALT_DP, variant type distribution, and Ti/Tv ratio.
    """
    if "final_tier" not in df.columns:
        return pl.DataFrame()

    agg_exprs: list[Any] = [pl.len().alias("n_variants")]

    # Mean VAF
    for vaf_col in ["DNA_VAF_mean", "RNA_VAF_mean"]:
        if vaf_col in df.columns:
            agg_exprs.append(pl.col(vaf_col).mean().alias(f"mean_{vaf_col.lower()}"))

    # Mean DP
    for dp_col in ["DNA_DP_mean", "RNA_DP_mean"]:
        if dp_col in df.columns:
            agg_exprs.append(pl.col(dp_col).mean().alias(f"mean_{dp_col.lower()}"))

    # Mean REF/ALT DP
    for ref_col in ["DNA_REF_DP_mean", "RNA_REF_DP_mean"]:
        if ref_col in df.columns:
            agg_exprs.append(pl.col(ref_col).mean().alias(f"mean_{ref_col.lower()}"))
    for alt_col in ["DNA_ALT_DP_mean", "RNA_ALT_DP_mean"]:
        if alt_col in df.columns:
            agg_exprs.append(pl.col(alt_col).mean().alias(f"mean_{alt_col.lower()}"))

    # Variant type counts
    if "variant_type" in df.columns:
        for vt in ["SNV", "INS", "DEL", "MNV"]:
            agg_exprs.append(
                (pl.col("variant_type") == vt).sum().alias(f"n_{vt}")
            )

    # Ti/Tv
    if "ti_tv" in df.columns:
        agg_exprs.append(pl.col("ti_tv").sum().alias("n_ti"))
        agg_exprs.append((~pl.col("ti_tv")).sum().alias("n_tv"))
        agg_exprs.append(
            (pl.col("ti_tv").sum() / (~pl.col("ti_tv")).sum()).alias("ti_tv_ratio")
        )

    # N_SUPPORT_CALLERS distribution
    if "N_SUPPORT_CALLERS" in df.columns:
        for c in range(1, 7):
            agg_exprs.append(
                (pl.col("N_SUPPORT_CALLERS") == c).sum().alias(f"n_callers_{c}")
            )

    summary = df.group_by("final_tier").agg(agg_exprs)

    # Sort by tier order
    tier_order_expr = pl.col("final_tier").replace_strict(
        {t: i for i, t in enumerate(FINAL_TIER_ORDER)},
        default=len(FINAL_TIER_ORDER),
    )
    summary = summary.with_columns(tier_order_expr.alias("_tier_order")).sort("_tier_order").drop("_tier_order")

    return summary
