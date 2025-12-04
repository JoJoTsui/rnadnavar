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


def tier_rule(caller_support_count: int, modalities: List[str]) -> str:
    """
    Assign a tier string based on caller support and modality support.
    T1 is strongest; higher numbers mean lower confidence.
    """
    has_dna = "DNA" in modalities
    has_rna = "RNA" in modalities

    if caller_support_count >= 2 and has_dna and has_rna:
        return "T1"
    if caller_support_count >= 2 and (has_dna or has_rna):
        return "T2"
    if caller_support_count == 1 and has_dna and has_rna:
        return "T3"
    if caller_support_count == 1 and (has_dna or has_rna):
        return "T4"
    return "T5"


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
        filt = var.FILTER if var.FILTER and var.FILTER != "." else "PASS"
        if not filt:
            filt = "PASS"
        if filt not in ["PASS", "LowQual", "StrandBias", "Clustered", "Somatic", "Germline", "Reference", "Artifact"]:
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
        })

    df = pd.DataFrame.from_records(records)
    if df.empty:
        return df

    df["tier"] = [tier_rule(c, m) for c, m in zip(df["caller_support_count"], df["modalities"])]
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
