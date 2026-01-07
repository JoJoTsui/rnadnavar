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
import re
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


def _parse_filters_normalized_legacy(filters_normalized_str: str) -> Dict[str, str]:
    """
    Parse legacy FILTERS_NORMALIZED string format to extract per-caller categories.
    
    Format: "DNA_strelka:Somatic|RNA_mutect2:Germline|..."
    Note: Separator is pipe (|) not semicolon (;)
    
    Returns:
        Dict mapping "FILTER_NORMALIZED_<Caller>_<Modality>" to category
    """
    result = {}
    if not filters_normalized_str:
        return result
    
    # Split by pipe (|) for multiple entries - this is the actual VCF format
    # Also try semicolon as fallback for compatibility
    if "|" in filters_normalized_str:
        entries = filters_normalized_str.split("|")
    else:
        entries = filters_normalized_str.split(";")
    
    for entry in entries:
        entry = entry.strip()
        if not entry or ":" not in entry:
            continue
        
        # Parse "DNA_strelka:Somatic" format
        parts = entry.split(":", 1)
        if len(parts) != 2:
            continue
        
        caller_spec, category = parts
        caller_spec = caller_spec.strip()
        category = category.strip()
        
        # Extract modality and caller name
        # Format: DNA_strelka, RNA_mutect2, etc.
        match = re.match(r"^(DNA|RNA)_(\w+)$", caller_spec, re.IGNORECASE)
        if not match:
            continue
        
        modality_raw = match.group(1).upper()
        caller_raw = match.group(2)
        
        # Normalize caller name (case-insensitive match)
        caller_normalized = None
        for known_caller in ["Strelka", "DeepSomatic", "Mutect2"]:
            if known_caller.lower() == caller_raw.lower():
                caller_normalized = known_caller
                break
        
        if not caller_normalized:
            continue
        
        # Map to FILTER_NORMALIZED_* field format
        # Assume DNA means DNA_TUMOR, RNA means RNA_TUMOR for tiering
        modality = f"{modality_raw}_TUMOR"
        field_name = f"FILTER_NORMALIZED_{caller_normalized}_{modality}"
        result[field_name] = category
    
    return result


def tier_variant(vcf_record, tiering_engine=None) -> Dict[str, Any]:
    """
    Compute tier for a single variant using category-concordant counting.
    
    This function implements the unified tiering logic:
    1. Extract final FILTER category
    2. Parse FILTER_NORMALIZED_* fields (or FILTERS_NORMALIZED fallback)
    3. Validate all expected fields are present (fast-fail)
    4. Count category-concordant callers (DNA/RNA)
    5. Check database support
    6. Compute final CxDy tier with quality score
    
    Args:
        vcf_record: cyvcf2.Variant object
        tiering_engine: Optional TieringEngine instance (created if None)
    
    Returns:
        Dictionary with tier information including:
        {
            "final_tier": "C1D1",
            "tier_quality": 160,
            "dna_caller_count": 3,
            "rna_caller_count": 2,
            "concordant_callers": {"DNA": ["Strelka", "Mutect2"], "RNA": ["Strelka"]},
            "supporting_databases": ["gnomAD", "COSMIC"],
            ...
        }
    
    Raises:
        ValueError: If FILTER_NORMALIZED_* fields are missing or invalid
        RuntimeError: If consensus-only variant detected (no individual callers)
    """
    if tiering_engine is None:
        tiering_engine = TieringEngine()
    
    # Step 1: Extract final FILTER category
    final_filter = vcf_record.FILTER
    if final_filter is None or final_filter == "PASS" or final_filter == ".":
        final_filter = "PASS"
    
    # Step 2: Extract INFO fields
    info_dict = dict(vcf_record.INFO) if hasattr(vcf_record, "INFO") else {}
    
    # Step 3: Parse FILTER_NORMALIZED_* fields
    filter_normalized_fields = {}
    
    # Primary: Look for individual FILTER_NORMALIZED_<Caller>_<Modality> fields
    for key in info_dict.keys():
        if key.startswith("FILTER_NORMALIZED_"):
            filter_normalized_fields[key] = info_dict.get(key)
    
    # Fallback: Parse legacy FILTERS_NORMALIZED string
    if not filter_normalized_fields:
        filters_normalized_str = info_dict.get("FILTERS_NORMALIZED", "")
        if filters_normalized_str:
            filter_normalized_fields = _parse_filters_normalized_legacy(filters_normalized_str)
    
    # Step 4: Validate required fields are present (FAST-FAIL)
    expected_callers = ["Strelka", "DeepSomatic", "Mutect2"]
    expected_modalities = ["DNA_TUMOR", "RNA_TUMOR"]
    
    found_fields = set()
    for caller in expected_callers:
        for modality in expected_modalities:
            field_name = f"FILTER_NORMALIZED_{caller}_{modality}"
            if field_name in filter_normalized_fields:
                found_fields.add(field_name)
    
    # Allow partial coverage but require at least some fields
    if not found_fields:
        variant_id = f"{vcf_record.CHROM}:{vcf_record.POS}:{vcf_record.REF}>{vcf_record.ALT[0]}"
        raise ValueError(
            f"Variant {variant_id}: No FILTER_NORMALIZED_* fields found. "
            f"Expected fields like FILTER_NORMALIZED_Strelka_DNA_TUMOR, etc. "
            f"Available INFO keys: {list(info_dict.keys())[:10]}..."
        )
    
    # Step 5: Check for consensus-only variants (FAST-FAIL)
    # These should have DNA_consensus/RNA_consensus in CALLERS_SUPPORT but no individual callers
    callers_support_str = info_dict.get("CALLERS_SUPPORT", "")
    if callers_support_str:
        # Check if only consensus callers are present
        callers_list = [c.strip() for c in callers_support_str.split("|")]
        has_consensus = any("consensus" in c.lower() for c in callers_list)
        has_individual = any("consensus" not in c.lower() for c in callers_list)
        
        if has_consensus and not has_individual:
            variant_id = f"{vcf_record.CHROM}:{vcf_record.POS}:{vcf_record.REF}>{vcf_record.ALT[0]}"
            raise RuntimeError(
                f"Variant {variant_id}: Consensus-only variant detected. "
                f"CALLERS_SUPPORT={callers_support_str}. "
                f"All variants must have at least one individual caller (Strelka/DeepSomatic/Mutect2). "
                f"Consensus-only variants should not exist in rescue VCFs."
            )
    
    # Step 6: Compute tier using TieringEngine
    tier_info = tiering_engine.compute_tier(
        final_filter=final_filter,
        filter_normalized_fields=filter_normalized_fields,
        info_dict=info_dict,
    )
    
    # Step 7: Add variant identification
    tier_info["variant_id"] = f"{vcf_record.CHROM}:{vcf_record.POS}:{vcf_record.REF}>{vcf_record.ALT[0]}"
    tier_info["chrom"] = vcf_record.CHROM
    tier_info["pos"] = vcf_record.POS
    tier_info["ref"] = vcf_record.REF
    tier_info["alt"] = vcf_record.ALT[0] if vcf_record.ALT else "."
    tier_info["filter_category"] = final_filter
    
    return tier_info


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


def tier_rescue_variants(rescue_vcf: Path) -> pd.DataFrame:
    """
    Load rescue VCF and compute CxDy tiers using unified tiering logic.
    
    This function processes rescue VCFs with fast-fail error handling:
    - Validates FILTER_NORMALIZED_* fields are present
    - Rejects consensus-only variants (should not exist)
    - Uses category-concordant caller counting
    - Returns DataFrame with tier metadata including quality scores
    
    Args:
        rescue_vcf: Path to rescue VCF file
    
    Returns:
        DataFrame with columns:
        - variant_id: "chr:pos:ref>alt"
        - chrom, pos, ref, alt: Variant coordinates
        - filter_category: Final FILTER classification
        - tier: CxDy tier (e.g., "C1D1", "C7D0")
        - tier_quality: Quality score (C1D1=160, C7D0=20)
        - dna_caller_count: Number of concordant DNA callers
        - rna_caller_count: Number of concordant RNA callers
        - concordant_callers: Dict {"DNA": [...], "RNA": [...]}
        - supporting_databases: List of database names
        - tier_display: Human-readable tier description
    
    Raises:
        ValueError: If FILTER_NORMALIZED_* fields missing from variants
        RuntimeError: If consensus-only variants encountered
        ImportError: If cyvcf2 not installed
    """
    try:
        from cyvcf2 import VCF
    except ImportError:
        raise ImportError(
            "cyvcf2 is required to parse VCF files for tiering. "
            "Install with: pip install cyvcf2"
        )
    
    # Initialize tiering engine once for all variants
    tiering_engine = TieringEngine()
    
    # Process all variants with fast-fail on first error
    records: List[Dict[str, Any]] = []
    vcf = VCF(str(rescue_vcf))
    
    for i, variant in enumerate(vcf):
        try:
            # Compute tier using new unified function (fast-fail)
            tier_info = tier_variant(variant, tiering_engine=tiering_engine)
            
            # Extract relevant fields for DataFrame
            record = {
                "variant_id": tier_info["variant_id"],
                "chrom": tier_info["chrom"],
                "pos": tier_info["pos"],
                "ref": tier_info["ref"],
                "alt": tier_info["alt"],
                "filter_category": tier_info["filter_category"],
                "tier": tier_info["final_tier"],
                "tier_quality": tier_info["tier_quality"],
                "caller_tier": tier_info["caller_tier"],
                "database_tier": tier_info["database_tier"],
                "dna_caller_count": tier_info["dna_caller_count"],
                "rna_caller_count": tier_info["rna_caller_count"],
                "concordant_callers_dna": ",".join(tier_info["concordant_callers"]["DNA"]),
                "concordant_callers_rna": ",".join(tier_info["concordant_callers"]["RNA"]),
                "supporting_databases": ",".join(tier_info["supporting_databases"]),
                "tier_display": tier_info.get("tier_display_short", tier_info["final_tier"]),
            }
            
            records.append(record)
            
        except (ValueError, RuntimeError) as e:
            # Fast-fail: Stop processing on first validation error
            vcf.close()
            raise type(e)(
                f"Error processing variant #{i+1} in {rescue_vcf.name}: {str(e)}"
            ) from e
    
    vcf.close()
    
    # Convert to DataFrame
    df = pd.DataFrame.from_records(records)
    
    if df.empty:
        return df
    
    # Sort by quality score (highest first) for prioritization
    df = df.sort_values("tier_quality", ascending=False).reset_index(drop=True)
    
    return df


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
