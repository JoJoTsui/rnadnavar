#!/usr/bin/env python3
"""
Category-Aware Caller Counting Module

This module implements category-aware caller counting where only callers that
vote for the SAME category as the final FILTER classification are counted
toward tier assignment.

Key principle: A variant's tier should reflect the strength of evidence
supporting its FINAL classification, not conflicting evidence.

Usage:
    from category_matcher import count_concordant_callers
    
    # Parse INFO fields from VCF
    dna_count, rna_count = count_concordant_callers(
        final_filter="Somatic",
        filter_normalized_strelka_dna="Somatic",
        filter_normalized_mutect2_dna="Somatic",
        filter_normalized_deepsomatic_dna="Somatic",
        filter_normalized_strelka_rna="Artifact",
        filter_normalized_mutect2_rna="Somatic",
        filter_normalized_deepsomatic_rna=None
    )
    # Returns: (3, 1) - 3 DNA callers agree, 1 RNA caller agrees
"""

from typing import Dict, List, Optional, Tuple
import re


# ============================================================================
# CALLER DEFINITIONS
# ============================================================================

# List of supported variant callers
VARIANT_CALLERS = ["Strelka", "DeepSomatic", "Mutect2"]

# Modality definitions
DNA_MODALITIES = ["DNA_TUMOR", "DNA_NORMAL"]
RNA_MODALITIES = ["RNA_TUMOR"]

# INFO field patterns for FILTER_NORMALIZED_<CALLER>_<MODALITY>
# Example: FILTER_NORMALIZED_Strelka_DNA_TUMOR, FILTER_NORMALIZED_Mutect2_RNA_TUMOR
FILTER_NORMALIZED_PATTERN = re.compile(
    r"FILTER_NORMALIZED_(?P<caller>\w+)_(?P<modality>DNA_TUMOR|DNA_NORMAL|RNA_TUMOR)"
)


# ============================================================================
# CATEGORY-AWARE CALLER COUNTING
# ============================================================================

def count_concordant_callers(
    final_filter: str,
    filter_normalized_fields: Dict[str, Optional[str]],
) -> Tuple[int, int]:
    """
    Count DNA and RNA callers that agree with final FILTER category.
    
    Only callers whose FILTER_NORMALIZED_<CALLER>_<MODALITY> field matches
    the final FILTER classification are counted as supporting evidence.
    
    Args:
        final_filter: Final FILTER category (e.g., "Somatic", "Germline")
        filter_normalized_fields: Dict mapping field names to values
            Example: {
                "FILTER_NORMALIZED_Strelka_DNA_TUMOR": "Somatic",
                "FILTER_NORMALIZED_Mutect2_DNA_TUMOR": "Somatic",
                "FILTER_NORMALIZED_Strelka_RNA_TUMOR": "Artifact",
                ...
            }
    
    Returns:
        Tuple of (dna_caller_count, rna_caller_count)
        
    Example:
        >>> fields = {
        ...     "FILTER_NORMALIZED_Strelka_DNA_TUMOR": "Somatic",
        ...     "FILTER_NORMALIZED_Mutect2_DNA_TUMOR": "Somatic",
        ...     "FILTER_NORMALIZED_DeepSomatic_DNA_TUMOR": "Somatic",
        ...     "FILTER_NORMALIZED_Strelka_RNA_TUMOR": "Artifact",
        ...     "FILTER_NORMALIZED_Mutect2_RNA_TUMOR": "Somatic",
        ... }
        >>> count_concordant_callers("Somatic", fields)
        (3, 1)  # 3 DNA agree, 1 RNA agrees
    """
    dna_count = 0
    rna_count = 0
    
    # Track which callers we've already counted (avoid double-counting)
    counted_dna_callers = set()
    counted_rna_callers = set()
    
    for field_name, field_value in filter_normalized_fields.items():
        # Skip None/empty values
        if not field_value:
            continue
        
        # Parse field name to extract caller and modality
        match = FILTER_NORMALIZED_PATTERN.match(field_name)
        if not match:
            continue
        
        caller = match.group("caller")
        modality = match.group("modality")
        
        # Check if this caller's vote matches final classification
        if field_value == final_filter:
            # Count DNA modality callers (DNA_TUMOR or DNA_NORMAL)
            if modality in DNA_MODALITIES:
                if caller not in counted_dna_callers:
                    dna_count += 1
                    counted_dna_callers.add(caller)
            
            # Count RNA modality callers (RNA_TUMOR)
            elif modality in RNA_MODALITIES:
                if caller not in counted_rna_callers:
                    rna_count += 1
                    counted_rna_callers.add(caller)
    
    return dna_count, rna_count


def parse_caller_votes_from_info(
    info_dict: Dict[str, str],
    final_filter: str,
) -> Dict[str, Dict[str, Optional[str]]]:
    """
    Parse all FILTER_NORMALIZED fields from INFO dict and organize by caller.
    
    Args:
        info_dict: Dictionary of INFO fields from VCF record
        final_filter: Final FILTER category
    
    Returns:
        Dictionary organized as:
        {
            "Strelka": {
                "DNA_TUMOR": "Somatic",
                "RNA_TUMOR": "Artifact"
            },
            "Mutect2": {
                "DNA_TUMOR": "Somatic",
                "RNA_TUMOR": "Somatic"
            },
            ...
        }
    """
    caller_votes = {caller: {} for caller in VARIANT_CALLERS}
    
    for field_name, field_value in info_dict.items():
        match = FILTER_NORMALIZED_PATTERN.match(field_name)
        if match:
            caller = match.group("caller")
            modality = match.group("modality")
            
            if caller in caller_votes:
                caller_votes[caller][modality] = field_value
    
    return caller_votes


def get_concordant_callers_by_modality(
    final_filter: str,
    filter_normalized_fields: Dict[str, Optional[str]],
) -> Dict[str, List[str]]:
    """
    Get list of concordant callers grouped by modality.
    
    Args:
        final_filter: Final FILTER category
        filter_normalized_fields: Dict mapping field names to values
    
    Returns:
        Dictionary with keys "DNA" and "RNA" containing lists of caller names
        
    Example:
        >>> get_concordant_callers_by_modality("Somatic", fields)
        {
            "DNA": ["Strelka", "Mutect2", "DeepSomatic"],
            "RNA": ["Mutect2"]
        }
    """
    dna_callers = []
    rna_callers = []
    
    counted_dna = set()
    counted_rna = set()
    
    for field_name, field_value in filter_normalized_fields.items():
        if not field_value:
            continue
        
        match = FILTER_NORMALIZED_PATTERN.match(field_name)
        if not match:
            continue
        
        caller = match.group("caller")
        modality = match.group("modality")
        
        if field_value == final_filter:
            if modality in DNA_MODALITIES and caller not in counted_dna:
                dna_callers.append(caller)
                counted_dna.add(caller)
            elif modality in RNA_MODALITIES and caller not in counted_rna:
                rna_callers.append(caller)
                counted_rna.add(caller)
    
    return {"DNA": dna_callers, "RNA": rna_callers}


def count_callers_from_vcf_record(vcf_record) -> Tuple[int, int]:
    """
    Extract concordant caller counts directly from a cyvcf2 VCF record.
    
    This function is designed to work with cyvcf2.Variant objects and
    extracts FILTER_NORMALIZED_* fields from the INFO attribute.
    
    Args:
        vcf_record: cyvcf2.Variant object with INFO fields
    
    Returns:
        Tuple of (dna_caller_count, rna_caller_count)
    
    Example:
        >>> import cyvcf2
        >>> vcf = cyvcf2.VCF("rescue.vcf.gz")
        >>> for variant in vcf:
        ...     dna, rna = count_callers_from_vcf_record(variant)
        ...     print(f"{variant.CHROM}:{variant.POS} - DNA:{dna}, RNA:{rna}")
    """
    # Get final FILTER value
    if hasattr(vcf_record, 'FILTER'):
        final_filter = vcf_record.FILTER
        if final_filter is None or final_filter == "PASS" or final_filter == ".":
            # Try to get from INFO if FILTER is PASS
            final_filter = vcf_record.INFO.get("FILTER", "PASS")
    else:
        final_filter = "PASS"
    
    # Extract all FILTER_NORMALIZED_* fields from INFO
    filter_normalized_fields = {}
    
    if hasattr(vcf_record, 'INFO'):
        info_dict = dict(vcf_record.INFO)
        for key in info_dict.keys():
            if key.startswith("FILTER_NORMALIZED_"):
                filter_normalized_fields[key] = info_dict.get(key)
    
    # Count concordant callers
    return count_concordant_callers(final_filter, filter_normalized_fields)


# ============================================================================
# VALIDATION AND DEBUGGING UTILITIES
# ============================================================================

def validate_caller_votes(
    filter_normalized_fields: Dict[str, Optional[str]],
) -> Dict[str, any]:
    """
    Validate and summarize caller votes for debugging.
    
    Args:
        filter_normalized_fields: Dict mapping field names to values
    
    Returns:
        Dictionary containing:
        - total_fields: Number of FILTER_NORMALIZED fields found
        - dna_votes: Dict of caller -> category for DNA
        - rna_votes: Dict of caller -> category for RNA
        - unique_categories: Set of unique categories voted
    """
    dna_votes = {}
    rna_votes = {}
    unique_categories = set()
    
    for field_name, field_value in filter_normalized_fields.items():
        if not field_value:
            continue
        
        match = FILTER_NORMALIZED_PATTERN.match(field_name)
        if not match:
            continue
        
        caller = match.group("caller")
        modality = match.group("modality")
        unique_categories.add(field_value)
        
        if modality in DNA_MODALITIES:
            dna_votes[caller] = field_value
        elif modality in RNA_MODALITIES:
            rna_votes[caller] = field_value
    
    return {
        "total_fields": len(filter_normalized_fields),
        "dna_votes": dna_votes,
        "rna_votes": rna_votes,
        "unique_categories": sorted(unique_categories),
    }


def print_caller_vote_summary(
    final_filter: str,
    filter_normalized_fields: Dict[str, Optional[str]],
) -> None:
    """
    Print a human-readable summary of caller votes for debugging.
    
    Args:
        final_filter: Final FILTER category
        filter_normalized_fields: Dict mapping field names to values
    """
    validation = validate_caller_votes(filter_normalized_fields)
    dna_count, rna_count = count_concordant_callers(final_filter, filter_normalized_fields)
    concordant = get_concordant_callers_by_modality(final_filter, filter_normalized_fields)
    
    print(f"Final FILTER: {final_filter}")
    print(f"Total FILTER_NORMALIZED fields: {validation['total_fields']}")
    print(f"\nDNA Caller Votes:")
    for caller, category in validation['dna_votes'].items():
        match = "✓" if category == final_filter else "✗"
        print(f"  {match} {caller}: {category}")
    
    print(f"\nRNA Caller Votes:")
    for caller, category in validation['rna_votes'].items():
        match = "✓" if category == final_filter else "✗"
        print(f"  {match} {caller}: {category}")
    
    print(f"\nConcordant Callers:")
    print(f"  DNA: {dna_count} callers {concordant['DNA']}")
    print(f"  RNA: {rna_count} callers {concordant['RNA']}")
    
    print(f"\nUnique Categories: {', '.join(validation['unique_categories'])}")


# ============================================================================
# LEGACY SUPPORT
# ============================================================================

def count_callers_legacy(
    dna_support: int,
    rna_support: int,
) -> Tuple[int, int]:
    """
    Legacy counting method using N_DNA_CALLERS_SUPPORT and N_RNA_CALLERS_SUPPORT.
    
    This method does NOT implement category-aware counting and should be
    used only for backward compatibility or when FILTER_NORMALIZED fields
    are not available.
    
    Args:
        dna_support: Value from N_DNA_CALLERS_SUPPORT INFO field
        rna_support: Value from N_RNA_CALLERS_SUPPORT INFO field
    
    Returns:
        Tuple of (dna_support, rna_support) unchanged
    
    Warning:
        This method counts ALL callers regardless of category concordance.
        Prefer count_concordant_callers() when FILTER_NORMALIZED fields exist.
    """
    return (dna_support or 0, rna_support or 0)


if __name__ == "__main__":
    # Example usage and testing
    print("=" * 80)
    print("CATEGORY-AWARE CALLER COUNTING - EXAMPLE")
    print("=" * 80)
    
    # Example 1: All callers agree on Somatic
    print("\nExample 1: All callers agree on Somatic")
    fields1 = {
        "FILTER_NORMALIZED_Strelka_DNA_TUMOR": "Somatic",
        "FILTER_NORMALIZED_Mutect2_DNA_TUMOR": "Somatic",
        "FILTER_NORMALIZED_DeepSomatic_DNA_TUMOR": "Somatic",
        "FILTER_NORMALIZED_Strelka_RNA_TUMOR": "Somatic",
        "FILTER_NORMALIZED_Mutect2_RNA_TUMOR": "Somatic",
    }
    print_caller_vote_summary("Somatic", fields1)
    
    # Example 2: Mixed votes - some Artifact, some Somatic
    print("\n" + "=" * 80)
    print("Example 2: Mixed votes - 2 DNA agree, 1 RNA agrees")
    fields2 = {
        "FILTER_NORMALIZED_Strelka_DNA_TUMOR": "Somatic",
        "FILTER_NORMALIZED_Mutect2_DNA_TUMOR": "Somatic",
        "FILTER_NORMALIZED_DeepSomatic_DNA_TUMOR": "Artifact",
        "FILTER_NORMALIZED_Strelka_RNA_TUMOR": "Artifact",
        "FILTER_NORMALIZED_Mutect2_RNA_TUMOR": "Somatic",
    }
    print_caller_vote_summary("Somatic", fields2)
    
    # Example 3: RNA_Edit with no caller support
    print("\n" + "=" * 80)
    print("Example 3: RNA_Edit with no caller support (C7 tier expected)")
    fields3 = {
        "FILTER_NORMALIZED_Strelka_DNA_TUMOR": "Artifact",
        "FILTER_NORMALIZED_Mutect2_DNA_TUMOR": "Artifact",
        "FILTER_NORMALIZED_Strelka_RNA_TUMOR": "Artifact",
    }
    print_caller_vote_summary("RNA_Edit", fields3)
    
    print("\n" + "=" * 80)
