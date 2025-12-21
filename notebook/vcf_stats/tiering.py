#!/usr/bin/env python3
"""
Rescue Variant Tiering

Utilities to compute tiers for variants in the final rescue VCF based on:
- Number of callers detecting/supporting the variant
- Modality support (DNA and/or RNA)

Designed for reuse from notebooks or scripts.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import random
import pandas as pd


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


def tier_rule(dna_callers: int, rna_callers: int) -> str:
    """
    Assign a fine-grained tier based on per-modality caller support counts.
    
    Tiers represent all caller support combinations:
    T1: DNA consensus (≥2) + RNA consensus (≥2) - Strongest: both modalities agree
    T2: DNA consensus (≥2) + RNA support (1)   - DNA strong + single RNA caller
    T3: DNA consensus (≥2) + RNA absent (0)    - DNA strong only
    T4: DNA single (1) + RNA support (≥1)      - Single DNA caller + RNA corroboration
    T5: DNA single (1) + RNA absent (0)        - Single DNA caller only
    T6: DNA absent (0) + RNA consensus (≥2)    - RNA consensus only
    T7: DNA absent (0) + RNA single (1)        - Single RNA caller only
    T8: DNA absent (0) + RNA absent (0)        - No caller support (rare edge case)
    
    Args:
        dna_callers: Number of DNA callers supporting the variant (0-3)
        rna_callers: Number of RNA callers supporting the variant (0-3)
    
    Returns:
        Tier string T1-T8 (T1 strongest, T8 weakest)
    """
    dna_callers = int(dna_callers) if dna_callers is not None else 0
    rna_callers = int(rna_callers) if rna_callers is not None else 0
    
    # T1: Both modalities with consensus support
    if dna_callers >= 2 and rna_callers >= 2:
        return "T1"
    
    # T2: DNA consensus + single RNA support
    if dna_callers >= 2 and rna_callers == 1:
        return "T2"
    
    # T3: DNA consensus only
    if dna_callers >= 2 and rna_callers == 0:
        return "T3"
    
    # T4: Single DNA caller + RNA support
    if dna_callers == 1 and rna_callers >= 1:
        return "T4"
    
    # T5: Single DNA caller only
    if dna_callers == 1 and rna_callers == 0:
        return "T5"
    
    # T6: RNA consensus only (no DNA)
    if dna_callers == 0 and rna_callers >= 2:
        return "T6"
    
    # T7: Single RNA caller only (no DNA)
    if dna_callers == 0 and rna_callers == 1:
        return "T7"
    
    # T8: No caller support (edge case) -> map to weakest defined tier T7
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
