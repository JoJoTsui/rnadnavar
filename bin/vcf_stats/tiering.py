#!/usr/bin/env python3
"""
Rescue Variant Tiering

Utilities to compute tiers for variants in the final rescue VCF based on:
- Category-aware caller support (C1-C7)
- Database evidence (D0-D1)
- Combined hybrid tiering (CxDy format)

Uses the TieringEngine for consistent tier assignment with the production pipeline.
"""

from __future__ import annotations

import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

# Add bin/common to path for shared config imports (sibling directory)
_bin_common_path = Path(__file__).parent.parent / "common"
if str(_bin_common_path) not in sys.path:
    sys.path.insert(0, str(_bin_common_path))

# Import from shared config
try:
    from tier_config import (
        CALLER_TIER_ORDER,
        CALLER_TIER_RULES,
        TIER_COLORS,
        TIER_ORDER,
        TIER_QUALITY_SCORES,
        get_tier_color,
        get_tier_quality,
    )
    from vcf_config import CATEGORY_ORDER, FILTER_FIELD_VALUES
except ImportError:
    # Fallback definitions
    TIER_ORDER = [f"C{c}D{d}" for c in range(1, 8) for d in [1, 0]]
    TIER_COLORS = {}
    TIER_QUALITY_SCORES = {}
    CALLER_TIER_ORDER = [f"C{i}" for i in range(1, 8)]
    CATEGORY_ORDER = [
        "Somatic",
        "Germline",
        "Reference",
        "Artifact",
        "RNAedit",
        "NoConsensus",
    ]
    FILTER_FIELD_VALUES = ["PASS"] + CATEGORY_ORDER

    def get_tier_quality(tier):
        return 0

    def get_tier_color(tier):
        return "#999999"


# Import tiering engine
from .tiering_engine import TieringEngine


@dataclass(frozen=True)
class VariantKey:
    chrom: str
    pos: int
    ref: str
    alt: str


def _parse_callers_from_info(info: Any) -> Tuple[List[str], int]:
    """
    Try to derive the list and count of callers from common INFO field patterns.
    Returns (caller_list, count). Falls back to empty list and 0.
    """
    candidates = []
    try:
        # Common keys seen across pipelines
        for key in [
            "CALLERS",
            "CALLER",
            "TOOLS",
            "SOURCES",
            "SRC",
            "SUPPORTED_BY",
            "CALLER_LIST",
        ]:
            val = info.get(key)
            if val:
                if isinstance(val, (list, tuple)):
                    candidates = [str(x) for x in val if x]
                else:
                    candidates = [
                        x.strip()
                        for x in str(val).replace(";", ",").split(",")
                        if x.strip()
                    ]
                break

        if not candidates:
            # Numeric count keys
            for key in ["NCALLERS", "NUM_CALLERS", "CALLER_COUNT", "N_CALLERS"]:
                val = info.get(key)
                if val is not None:
                    try:
                        return [], int(val)
                    except Exception:
                        pass
    except Exception:
        pass

    # Deduplicate and normalize names
    callers = []
    seen = set()
    for c in candidates:
        low = c.lower()
        if low in seen:
            continue
        seen.add(low)
        callers.append(c)

    return callers, len(callers)


def _parse_modalities_from_info(info: Any) -> List[str]:
    """
    Attempt to infer modality support from INFO fields.
    Returns subset of ["DNA", "RNA"].
    """
    modalities: List[str] = []
    try:
        # Direct list
        for key in ["MODALITIES", "MODAL", "SUPPORT_MODALITIES"]:
            val = info.get(key)
            if val:
                items = [
                    x.strip().upper()
                    for x in str(val).replace(";", ",").split(",")
                    if x.strip()
                ]
                if any("DNA" in x for x in items):
                    modalities.append("DNA")
                if any("RNA" in x for x in items):
                    modalities.append("RNA")
                return sorted(set(modalities))

        # Boolean or numeric flags
        dna_flags = ["DNA_SUPPORT", "HAS_DNA", "DNA", "DNA_ALT_READS", "DNA_AF"]
        rna_flags = ["RNA_SUPPORT", "HAS_RNA", "RNA", "RNA_ALT_READS", "RNA_AF"]

        for key in dna_flags:
            if key in info and info.get(key) not in (None, 0, False, "0"):
                modalities.append("DNA")
                break
        for key in rna_flags:
            if key in info and info.get(key) not in (None, 0, False, "0"):
                modalities.append("RNA")
                break
    except Exception:
        pass

    return sorted(set(modalities))


def _per_modality_caller_counts(info: Any) -> Tuple[int, int]:
    """
    Parse counts of DNA-supporting callers and RNA-supporting callers from VCF INFO.
    Returns (dna_callers, rna_callers). Fallback to 0 if not found.
    """
    dna, rna = 0, 0
    try:
        # Direct per-modality caller support count fields
        dna_val = info.get("N_DNA_CALLERS_SUPPORT")
        if dna_val is not None:
            dna = int(dna_val)

        rna_val = info.get("N_RNA_CALLERS_SUPPORT")
        if rna_val is not None:
            rna = int(rna_val)
    except Exception:
        pass
    return dna, rna


def _check_database_support(info: Any) -> bool:
    """
    Check if variant has database support (gnomAD, COSMIC, REDIportal, DARNED).

    Returns:
        True if variant has database support
    """
    try:
        # gnomAD
        gnomad_af = info.get("gnomAD_AF") or info.get("GNOMAD_AF")
        if gnomad_af is not None:
            try:
                if float(gnomad_af) > 0.001:
                    return True
            except (ValueError, TypeError):
                pass

        # COSMIC
        cosmic_cnt = info.get("COSMIC_CNT") or info.get("COSMIC_COUNT")
        if cosmic_cnt is not None:
            try:
                if int(cosmic_cnt) > 0:
                    return True
            except (ValueError, TypeError):
                pass

        # REDIportal / DARNED
        if info.get("REDIportal") or info.get("REDIPORTAL") or info.get("DARNED"):
            return True

    except Exception:
        pass

    return False


def compute_tier(
    dna_callers: int, rna_callers: int, has_database_support: bool = False
) -> str:
    """
    Compute tier using CxDy hybrid tiering system.

    TIERING SYSTEM:
    - Caller tiers (C1-C7): Based on DNA/RNA caller counts
      C1: ≥2 DNA + ≥2 RNA (both strong)
      C2: ≥2 DNA + (0 or 1) RNA (DNA strong)
      C3: ≥2 RNA + (0 or 1) DNA (RNA strong)
      C4: 1 DNA + 1 RNA (both weak)
      C5: 1 DNA + 0 RNA (DNA only weak)
      C6: 0 DNA + 1 RNA (RNA only weak)
      C7: 0 DNA + 0 RNA (no caller support)

    - Database tiers (D0-D1):
      D1: Has database support (gnomAD/COSMIC/REDIportal/DARNED)
      D0: No database support

    - Final tier: CxDy (e.g., C1D1, C2D0, C7D1)

    Args:
        dna_callers: Number of DNA callers supporting the variant
        rna_callers: Number of RNA callers supporting the variant
        has_database_support: Whether variant has database annotation

    Returns:
        Tier string in CxDy format (e.g., "C1D1", "C7D0")
    """
    engine = TieringEngine()
    return engine.compute_tier_simple(
        dna_caller_count=int(dna_callers) if dna_callers is not None else 0,
        rna_caller_count=int(rna_callers) if rna_callers is not None else 0,
        has_database_support=has_database_support,
    )


def load_rescue_variants(rescue_vcf: Path) -> pd.DataFrame:
    """
    Load per-variant details from a rescue VCF into a DataFrame with CxDy tiering.

    Columns: chrom, pos, ref, alt, filter_category, callers, caller_support_count,
             modalities, dna_callers, rna_callers, has_db_support, tier, tier_quality
    """
    try:
        from cyvcf2 import VCF
    except ImportError:
        raise RuntimeError("cyvcf2 is required to parse VCF files for tiering")

    records: List[Dict[str, Any]] = []
    vcf = VCF(str(rescue_vcf))

    for var in vcf:
        callers, n_callers = _parse_callers_from_info(var.INFO)
        modalities = _parse_modalities_from_info(var.INFO)
        dna_callers, rna_callers = _per_modality_caller_counts(var.INFO)
        has_db_support = _check_database_support(var.INFO)

        filt = var.FILTER if var.FILTER and var.FILTER != "." else "PASS"
        if not filt:
            filt = "PASS"

        # Validate filter against known categories
        if filt not in FILTER_FIELD_VALUES:
            # Try case-insensitive match
            filt_lower = filt.lower()
            matched = False
            for valid_cat in FILTER_FIELD_VALUES:
                if valid_cat.lower() == filt_lower:
                    filt = valid_cat
                    matched = True
                    break
            if not matched:
                filt = "Other"

        records.append(
            {
                "chrom": var.CHROM,
                "pos": int(var.POS),
                "ref": var.REF,
                "alt": var.ALT[0]
                if isinstance(var.ALT, (list, tuple)) and var.ALT
                else var.ALT,
                "filter_category": filt,
                "callers": callers,
                "caller_support_count": n_callers,
                "modalities": modalities,
                "dna_callers": dna_callers,
                "rna_callers": rna_callers,
                "has_db_support": has_db_support,
            }
        )

    df = pd.DataFrame.from_records(records)
    if df.empty:
        return df

    # Assign CxDy tiers based on per-modality caller support and database evidence
    df["tier"] = [
        compute_tier(dna, rna, has_db)
        for dna, rna, has_db in zip(
            df["dna_callers"], df["rna_callers"], df["has_db_support"]
        )
    ]

    # Add tier quality scores
    df["tier_quality"] = df["tier"].apply(get_tier_quality)

    # Create variant ID
    df["variant_id"] = (
        df["chrom"].astype(str)
        + ":"
        + df["pos"].astype(str)
        + ":"
        + df["ref"].astype(str)
        + ">"
        + df["alt"].astype(str)
    )
    return df


def tier_rescue_variants(rescue_vcf: Path) -> pd.DataFrame:
    """Public API: load rescue VCF and compute CxDy tiers."""
    return load_rescue_variants(Path(rescue_vcf))


def sample_tier_representatives(
    tiered_variants: pd.DataFrame,
    k: int = 5,
    random_state: Optional[int] = 42,
) -> pd.DataFrame:
    """
    For each filter_category and tier, randomly sample up to k variants.
    Returns a concatenated DataFrame of representatives with columns preserved.
    """
    if tiered_variants.empty:
        return tiered_variants

    rng = random.Random(random_state)
    parts: List[pd.DataFrame] = []
    for (cat, tier), group in tiered_variants.groupby(
        ["filter_category", "tier"], dropna=False
    ):
        if len(group) <= k:
            parts.append(group)
        else:
            idx = list(group.index)
            rng.shuffle(idx)
            take = idx[:k]
            parts.append(group.loc[take])

    if parts:
        return pd.concat(parts, axis=0).reset_index(drop=True)
    return pd.DataFrame(columns=list(tiered_variants.columns))


def get_tier_summary(tiered_df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate summary statistics for tiered variants.

    Args:
        tiered_df: DataFrame from tier_rescue_variants()

    Returns:
        DataFrame with tier distribution summary
    """
    if tiered_df.empty:
        return pd.DataFrame()

    summary = (
        tiered_df.groupby(["filter_category", "tier"])
        .agg(
            {
                "variant_id": "count",
                "dna_callers": "mean",
                "rna_callers": "mean",
                "tier_quality": "first",
            }
        )
        .reset_index()
    )

    summary.columns = [
        "Category",
        "Tier",
        "Count",
        "Mean_DNA_Callers",
        "Mean_RNA_Callers",
        "Tier_Quality",
    ]
    summary = summary.sort_values(["Category", "Tier_Quality"], ascending=[True, False])

    return summary
