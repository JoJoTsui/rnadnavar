#!/usr/bin/env python3
"""
Rescue Variant Tiering

Utilities to compute tiers for variants in the final rescue VCF based on:
- Category-aware caller support (C1-C7)
- Database evidence (D0-D1)
- Combined hybrid tiering (CxDy format)

Uses the new TieringEngine for consistent tier assignment.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import random
import pandas as pd
import sys

# Import new tiering engine
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
            "CALLERS", "CALLER", "TOOLS", "SOURCES", "SRC", "SUPPORTED_BY", "CALLER_LIST"
        ]:
            val = info.get(key)
            if val:
                if isinstance(val, (list, tuple)):
                    candidates = [str(x) for x in val if x]
                else:
                    candidates = [x.strip() for x in str(val).replace(";", ",").split(",") if x.strip()]
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
                items = [x.strip().upper() for x in str(val).replace(";", ",").split(",") if x.strip()]
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


def tier_rule(dna_callers: int, rna_callers: int, has_database_support: bool = False) -> str:
    """
    Assign tier using new hybrid tiering system (CxDy format).
    
    NEW TIERING SYSTEM:
    - Caller tiers (C1-C7): Based on DNA/RNA caller counts
      C1: ≥2 DNA + ≥2 RNA (both strong)
      C2: ≥2 DNA + (0 or 1) RNA (DNA strong)
      C3: ≥2 RNA + (0 or 1) DNA (RNA strong)
      C4: 1 DNA + 1 RNA (both weak)
      C5: 1 DNA + 0 RNA (DNA only weak)
      C6: 0 DNA + 1 RNA (RNA only weak)
      C7: 0 DNA + 0 RNA (no caller support)
    
    - Database tiers (D0-D1):
      D1: Has database support (gnomAD/COSMIC/REDIportal)
      D0: No database support
    
    - Final tier: CxDy (e.g., C1D1, C2D0, C7D1)
    
    Args:
        dna_callers: Number of DNA callers supporting the variant
        rna_callers: Number of RNA callers supporting the variant
        has_database_support: Whether variant has database annotation
    
    Returns:
        Tier string in CxDy format (e.g., "C1D1", "C7D0")
    
    Note:
        For backward compatibility, this uses the legacy counting method
        (N_DNA_CALLERS_SUPPORT, N_RNA_CALLERS_SUPPORT).
        For full category-aware counting, use TieringEngine.compute_tier() directly.
    """
    engine = TieringEngine()
    return engine.compute_tier_simple(
        dna_caller_count=int(dna_callers) if dna_callers is not None else 0,
        rna_caller_count=int(rna_callers) if rna_callers is not None else 0,
        has_database_support=has_database_support
    )


def tier_rule_legacy(dna_callers: int, rna_callers: int) -> str:
    """
    DEPRECATED: Legacy tier rule using T1-T8 system.
    
    This function is kept for backward compatibility only.
    New code should use tier_rule() which returns CxDy format.
    
    Returns:
        Tier string T1-T8 (old system)
    """
    dna_callers = int(dna_callers) if dna_callers is not None else 0
    rna_callers = int(rna_callers) if rna_callers is not None else 0
    
    # Legacy T1-T8 tier assignment
    if dna_callers >= 2 and rna_callers >= 2:
        return "T1"
    if dna_callers >= 2 and rna_callers == 1:
        return "T2"
    if dna_callers >= 2 and rna_callers == 0:
        return "T3"
    if dna_callers == 1 and rna_callers >= 1:
        return "T4"
    if dna_callers == 1 and rna_callers == 0:
        return "T5"
    if dna_callers == 0 and rna_callers >= 2:
        return "T6"
    if dna_callers == 0 and rna_callers == 1:
        return "T7"
    return "T7"


def load_rescue_variants(rescue_vcf: Path) -> pd.DataFrame:
    """
    Load per-variant details from a rescue VCF into a DataFrame.
    Columns: chrom, pos, ref, alt, filter_category, callers, caller_support_count, modalities
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
        filt = var.FILTER if var.FILTER and var.FILTER != "." else "PASS"
        if not filt:
            filt = "PASS"
        if filt not in [
            "PASS", "LowQual", "StrandBias", "Clustered",
            "Somatic", "Germline", "Reference", "Artifact",
            "RNA_Edit", "NoConsensus"
        ]:
            filt = "Other"

        records.append({
            "chrom": var.CHROM,
            "pos": int(var.POS),
            "ref": var.REF,
            "alt": var.ALT[0] if isinstance(var.ALT, (list, tuple)) and var.ALT else var.ALT,
            "filter_category": filt,
            "callers": callers,
            "caller_support_count": n_callers,
            "modalities": modalities,
            "dna_callers": dna_callers,
            "rna_callers": rna_callers,
        })

    df = pd.DataFrame.from_records(records)
    if df.empty:
        return df

    # Assign fine-grained tiers based on per-modality caller support
    df["tier"] = [tier_rule(dna, rna) for dna, rna in zip(df["dna_callers"], df["rna_callers"])]
    df["variant_id"] = (
        df["chrom"].astype(str)
        + ":" + df["pos"].astype(str)
        + ":" + df["ref"].astype(str)
        + ">" + df["alt"].astype(str)
    )
    return df


def tier_rescue_variants(rescue_vcf: Path) -> pd.DataFrame:
    """Public API: load rescue VCF and compute tiers."""
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
    for (cat, tier), group in tiered_variants.groupby(["filter_category", "tier"], dropna=False):
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
