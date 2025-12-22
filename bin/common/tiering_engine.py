#!/usr/bin/env python3
"""
Tiering Engine for Hybrid Variant Tiering System

This module implements the hybrid tiering logic that combines:
1. Caller Support (C1-C7): Based on category-concordant DNA/RNA caller counts
2. Database Evidence (D0-D1): Based on presence in external databases

The final tier uses CxDy notation for clarity (e.g., C1D1, C7D0).

** UNIFIED INFRASTRUCTURE **
This module is the canonical tiering engine used by both Nextflow pipeline post-processing
and Jupyter notebook analysis. It implements the same tier assignment logic in both contexts.

Usage:
    from tiering_engine import TieringEngine
    
    # Create engine instance
    engine = TieringEngine()
    
    # Compute tier from VCF record
    tier = engine.compute_tier_from_vcf_record(variant)
    
    # Or compute from parsed data
    tier_info = engine.compute_tier(
        final_filter="Somatic",
        filter_normalized_fields={...},
        info_dict={...}
    )

Dependencies:
    - tier_config: CALLER_TIER_RULES, tier definitions (centralized)
    - category_matcher: count_concordant_callers, category-aware counting
    - database_checker: compute_database_tier, database evidence detection
"""

import sys
from pathlib import Path
from typing import Dict, Optional, Any, Tuple

# Local imports from bin/common
from tier_config import (
    CALLER_TIER_RULES,
    CALLER_TIER_ORDER,
    DATABASE_TIER_ORDER,
    TIER_ORDER,
    TIER_DISPLAY_NAMES,
    TIER_SHORT_NAMES,
    TIER_QUALITY_SCORES,
    get_tier_quality,
    get_tier_color,
    get_tier_display_name,
    parse_tier,
)

from category_matcher import (
    count_concordant_callers,
    get_concordant_callers_by_modality,
    validate_caller_votes,
)

from database_checker import (
    check_database_support,
    compute_database_tier,
    get_supporting_databases,
    get_database_details,
)


# ============================================================================
# TIERING ENGINE
# ============================================================================

class TieringEngine:
    """
    Main tiering engine that combines caller support and database evidence.
    
    This class provides methods to:
    - Compute caller support tier (C1-C7) using category-aware counting
    - Compute database evidence tier (D0-D1)
    - Combine into final tier (CxDy format)
    - Extract tier information from VCF records
    
    Example:
        >>> engine = TieringEngine()
        >>> tier_info = engine.compute_tier(
        ...     final_filter="Somatic",
        ...     filter_normalized_fields={
        ...         "FILTER_NORMALIZED_Strelka_DNA_TUMOR": "Somatic",
        ...         "FILTER_NORMALIZED_Mutect2_DNA_TUMOR": "Somatic",
        ...     },
        ...     info_dict={"gnomAD_AF": 0.0, "COSMIC_CNT": 5}
        ... )
        >>> print(tier_info["final_tier"])  # "C2D1"
    """
    
    def __init__(self):
        """Initialize the tiering engine."""
        self.caller_tier_rules = CALLER_TIER_RULES
        self.caller_tier_order = CALLER_TIER_ORDER
        self.database_tier_order = DATABASE_TIER_ORDER
    
    def compute_caller_tier(
        self,
        dna_caller_count: int,
        rna_caller_count: int
    ) -> str:
        """
        Compute caller support tier (C1-C7) from caller counts.
        
        Uses mutually exclusive tier assignment - first match wins.
        
        Args:
            dna_caller_count: Number of DNA callers supporting final category
            rna_caller_count: Number of RNA callers supporting final category
        
        Returns:
            Caller tier (C1-C7)
        
        Example:
            >>> engine.compute_caller_tier(3, 2)
            "C1"  # ≥2 DNA + ≥2 RNA
            
            >>> engine.compute_caller_tier(2, 0)
            "C2"  # ≥2 DNA + 0 RNA
        """
        # Check tiers in priority order (C1 → C7)
        for tier in self.caller_tier_order:
            condition = self.caller_tier_rules[tier]["condition"]
            if condition(dna_caller_count, rna_caller_count):
                return tier
        
        # Fallback (should never reach here due to C7 matching 0,0)
        return "C7"
    
    def compute_final_tier(
        self,
        caller_tier: str,
        database_tier: str
    ) -> str:
        """
        Combine caller and database tiers into final CxDy tier.
        
        Args:
            caller_tier: Caller tier (C1-C7)
            database_tier: Database tier (D0-D1)
        
        Returns:
            Final tier in CxDy format (e.g., "C1D1", "C7D0")
        """
        return f"{caller_tier}{database_tier}"
    
    def compute_tier(
        self,
        final_filter: str,
        filter_normalized_fields: Dict[str, Optional[str]],
        info_dict: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """
        Compute complete tier information from variant data.
        
        Args:
            final_filter: Final FILTER category (e.g., "Somatic")
            filter_normalized_fields: Dict of FILTER_NORMALIZED_* fields
            info_dict: Dict of INFO fields for database checking
        
        Returns:
            Dictionary containing:
            {
                "final_tier": "C1D1",
                "caller_tier": "C1",
                "database_tier": "D1",
                "dna_caller_count": 3,
                "rna_caller_count": 2,
                "concordant_callers": {"DNA": ["Strelka", ...], "RNA": [...]},
                "supporting_databases": ["gnomAD", "COSMIC"],
                "tier_quality": 140,
                "tier_display": "C1D1: ≥2 DNA + ≥2 RNA, Database",
            }
        """
        # Step 1: Count category-concordant callers
        dna_count, rna_count = count_concordant_callers(
            final_filter,
            filter_normalized_fields
        )
        
        # Step 2: Get concordant caller names
        concordant_callers = get_concordant_callers_by_modality(
            final_filter,
            filter_normalized_fields
        )
        
        # Step 3: Compute caller tier
        caller_tier = self.compute_caller_tier(dna_count, rna_count)
        
        # Step 4: Compute database tier
        if info_dict is not None:
            database_tier = compute_database_tier(info_dict=info_dict)
            supporting_dbs = get_supporting_databases(info_dict)
            db_details = get_database_details(info_dict)
        else:
            database_tier = "D0"
            supporting_dbs = []
            db_details = {}
        
        # Step 5: Combine into final tier
        final_tier = self.compute_final_tier(caller_tier, database_tier)
        
        # Step 6: Build result dictionary
        result = {
            "final_tier": final_tier,
            "caller_tier": caller_tier,
            "database_tier": database_tier,
            "dna_caller_count": dna_count,
            "rna_caller_count": rna_count,
            "concordant_callers": concordant_callers,
            "supporting_databases": supporting_dbs,
            "database_details": db_details,
            "tier_quality": get_tier_quality(final_tier),
            "tier_display": get_tier_display_name(final_tier),
            "tier_display_short": get_tier_display_name(final_tier, short=True),
            "tier_color": get_tier_color(final_tier),
        }
        
        return result
    
    def compute_tier_from_vcf_record(self, vcf_record) -> Dict[str, Any]:
        """
        Compute tier information directly from a cyvcf2 VCF record.
        
        Args:
            vcf_record: cyvcf2.Variant object
        
        Returns:
            Dictionary with tier information (same as compute_tier)
        
        Example:
            >>> import cyvcf2
            >>> vcf = cyvcf2.VCF("rescue.vcf.gz")
            >>> engine = TieringEngine()
            >>> for variant in vcf:
            ...     tier_info = engine.compute_tier_from_vcf_record(variant)
            ...     print(f"{variant.CHROM}:{variant.POS} - {tier_info['final_tier']}")
        """
        # Extract final FILTER value
        if hasattr(vcf_record, 'FILTER'):
            final_filter = vcf_record.FILTER
            if final_filter is None or final_filter == "PASS" or final_filter == ".":
                # Try to get from INFO if FILTER is PASS
                final_filter = vcf_record.INFO.get("FILTER", "PASS")
        else:
            final_filter = "PASS"
        
        # Extract all FILTER_NORMALIZED_* fields
        filter_normalized_fields = {}
        info_dict = {}
        
        if hasattr(vcf_record, 'INFO'):
            info_dict = dict(vcf_record.INFO)
            for key, value in info_dict.items():
                if key.startswith("FILTER_NORMALIZED_"):
                    filter_normalized_fields[key] = value
        
        # Compute tier
        return self.compute_tier(
            final_filter=final_filter,
            filter_normalized_fields=filter_normalized_fields,
            info_dict=info_dict
        )
    
    def compute_tier_simple(
        self,
        dna_caller_count: int,
        rna_caller_count: int,
        has_database_support: bool = False
    ) -> str:
        """
        Simplified tier computation from caller counts and database flag.
        
        Use this when you already have caller counts and just need the tier string.
        
        Args:
            dna_caller_count: Number of DNA callers
            rna_caller_count: Number of RNA callers
            has_database_support: Whether variant has database support
        
        Returns:
            Final tier string (e.g., "C1D1")
        
        Example:
            >>> engine.compute_tier_simple(3, 2, has_database_support=True)
            "C1D1"
        """
        caller_tier = self.compute_caller_tier(dna_caller_count, rna_caller_count)
        database_tier = "D1" if has_database_support else "D0"
        return self.compute_final_tier(caller_tier, database_tier)


# ============================================================================
# BATCH PROCESSING UTILITIES
# ============================================================================

def compute_tiers_for_vcf(vcf_path: str, limit: Optional[int] = None) -> list:
    """
    Compute tiers for all variants in a VCF file.
    
    Args:
        vcf_path: Path to VCF file
        limit: Optional limit on number of variants to process
    
    Returns:
        List of tier information dictionaries
    
    Example:
        >>> tiers = compute_tiers_for_vcf("rescue.vcf.gz", limit=100)
        >>> for tier_info in tiers:
        ...     print(f"{tier_info['variant_id']}: {tier_info['final_tier']}")
    """
    try:
        import cyvcf2
    except ImportError:
        raise ImportError("cyvcf2 is required for VCF processing. Install with: pip install cyvcf2")
    
    engine = TieringEngine()
    results = []
    
    vcf = cyvcf2.VCF(vcf_path)
    
    for i, variant in enumerate(vcf):
        if limit is not None and i >= limit:
            break
        
        tier_info = engine.compute_tier_from_vcf_record(variant)
        tier_info["variant_id"] = f"{variant.CHROM}:{variant.POS}:{variant.REF}>{variant.ALT[0]}"
        tier_info["chrom"] = variant.CHROM
        tier_info["pos"] = variant.POS
        tier_info["ref"] = variant.REF
        tier_info["alt"] = variant.ALT[0] if variant.ALT else "."
        
        results.append(tier_info)
    
    vcf.close()
    
    return results


def compute_tiers_to_dataframe(vcf_path: str, limit: Optional[int] = None):
    """
    Compute tiers for VCF and return as pandas DataFrame.
    
    Args:
        vcf_path: Path to VCF file
        limit: Optional limit on number of variants
    
    Returns:
        pandas DataFrame with tier information
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required. Install with: pip install pandas")
    
    tiers = compute_tiers_for_vcf(vcf_path, limit=limit)
    
    # Flatten nested structures for DataFrame
    flat_tiers = []
    for tier_info in tiers:
        flat = {
            "variant_id": tier_info["variant_id"],
            "chrom": tier_info["chrom"],
            "pos": tier_info["pos"],
            "ref": tier_info["ref"],
            "alt": tier_info["alt"],
            "final_tier": tier_info["final_tier"],
            "caller_tier": tier_info["caller_tier"],
            "database_tier": tier_info["database_tier"],
            "dna_caller_count": tier_info["dna_caller_count"],
            "rna_caller_count": tier_info["rna_caller_count"],
            "tier_quality": tier_info["tier_quality"],
            "supporting_databases": ",".join(tier_info["supporting_databases"]),
        }
        flat_tiers.append(flat)
    
    return pd.DataFrame(flat_tiers)


# ============================================================================
# VALIDATION AND DEBUGGING
# ============================================================================

def print_tier_summary(tier_info: Dict[str, Any]) -> None:
    """
    Print human-readable summary of tier information.
    
    Args:
        tier_info: Tier information dictionary from compute_tier()
    """
    print(f"Final Tier: {tier_info['final_tier']}")
    print(f"  Caller Tier: {tier_info['caller_tier']}")
    print(f"  Database Tier: {tier_info['database_tier']}")
    print(f"  Quality Score: {tier_info['tier_quality']}")
    print(f"\nCaller Support:")
    print(f"  DNA Callers: {tier_info['dna_caller_count']} {tier_info['concordant_callers']['DNA']}")
    print(f"  RNA Callers: {tier_info['rna_caller_count']} {tier_info['concordant_callers']['RNA']}")
    print(f"\nDatabase Support:")
    dbs = tier_info['supporting_databases']
    print(f"  Supporting Databases: {', '.join(dbs) if dbs else 'None'}")
    
    if tier_info.get('database_details'):
        details = tier_info['database_details']
        if details.get('gnomAD_AF') is not None:
            print(f"    gnomAD AF: {details['gnomAD_AF']}")
        if details.get('COSMIC_CNT') is not None:
            print(f"    COSMIC Count: {details['COSMIC_CNT']}")
        if details.get('REDIportal'):
            print(f"    REDIportal: Yes")
        if details.get('DARNED'):
            print(f"    DARNED: Yes")


if __name__ == "__main__":
    # Example usage and testing
    print("=" * 80)
    print("TIERING ENGINE - EXAMPLES")
    print("=" * 80)
    
    engine = TieringEngine()
    
    # Example 1: High quality somatic variant (C1D1)
    print("\nExample 1: High quality somatic variant")
    tier_info1 = engine.compute_tier(
        final_filter="Somatic",
        filter_normalized_fields={
            "FILTER_NORMALIZED_Strelka_DNA_TUMOR": "Somatic",
            "FILTER_NORMALIZED_Mutect2_DNA_TUMOR": "Somatic",
            "FILTER_NORMALIZED_DeepSomatic_DNA_TUMOR": "Somatic",
            "FILTER_NORMALIZED_Strelka_RNA_TUMOR": "Somatic",
            "FILTER_NORMALIZED_Mutect2_RNA_TUMOR": "Somatic",
        },
        info_dict={"gnomAD_AF": 0.0, "COSMIC_CNT": 15}
    )
    print_tier_summary(tier_info1)
    
    # Example 2: DNA-strong variant without database support (C2D0)
    print("\n" + "=" * 80)
    print("Example 2: DNA-strong novel variant")
    tier_info2 = engine.compute_tier(
        final_filter="Somatic",
        filter_normalized_fields={
            "FILTER_NORMALIZED_Strelka_DNA_TUMOR": "Somatic",
            "FILTER_NORMALIZED_Mutect2_DNA_TUMOR": "Somatic",
            "FILTER_NORMALIZED_DeepSomatic_DNA_TUMOR": "Artifact",
            "FILTER_NORMALIZED_Strelka_RNA_TUMOR": "Artifact",
        },
        info_dict={"gnomAD_AF": 0.0, "COSMIC_CNT": 0}
    )
    print_tier_summary(tier_info2)
    
    # Example 3: RNA editing site (C7D1)
    print("\n" + "=" * 80)
    print("Example 3: RNA editing site")
    tier_info3 = engine.compute_tier(
        final_filter="RNA_Edit",
        filter_normalized_fields={
            "FILTER_NORMALIZED_Strelka_DNA_TUMOR": "Artifact",
            "FILTER_NORMALIZED_Mutect2_DNA_TUMOR": "Artifact",
        },
        info_dict={"REDIportal": True}
    )
    print_tier_summary(tier_info3)
    
    # Example 4: Simplified tier computation
    print("\n" + "=" * 80)
    print("Example 4: Simplified tier computation")
    tier = engine.compute_tier_simple(
        dna_caller_count=2,
        rna_caller_count=1,
        has_database_support=False
    )
    print(f"Simple tier: {tier}")
    print(f"Description: {get_tier_display_name(tier, short=True)}")
    
    print("\n" + "=" * 80)
