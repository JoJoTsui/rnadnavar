#!/usr/bin/env python3
"""
Variant Classifier Module

Functions to classify variants from different callers (Strelka, DeepSomatic, Mutect2)
into biological categories: Somatic, Germline, Reference, Artifact, RNA_Edit, and NoConsensus.

This module supports classification at multiple processing stages:
- Raw caller outputs: Biological classification (Somatic/Germline/Reference/Artifact)
- Consensus/Rescue VCFs: FILTER-based classification
- Annotated VCFs: Detection of special annotations (RNA_Edit, Cosmic_gnomAD)
"""

from typing import Dict, Optional
from cyvcf2 import VCF


def classify_strelka_variant(filter_val, nt_val, normal_dp):
    """
    Classifies a Strelka variant into Somatic, Germline, Reference, or Artifact.

    Args:
        filter_val (str): The value from the VCF FILTER column.
                        Some parsers return None for PASS, Strelka writes "PASS" string.
        nt_val (str): The value of the INFO/NT field (e.g., 'ref', 'het', 'hom').
        normal_dp (int): The read depth of the Normal sample (FORMAT/DP).

    Returns:
        str: One of ['Somatic', 'Germline', 'Reference', 'Artifact']
    """

    # 1. Somatic: Strictly relies on PASS filter
    # Some parsers return None for PASS, Strelka writes "PASS" string.
    if filter_val == "PASS" or filter_val is None:
        # Check NT field for LOF
        if nt_val in ["somatic", "het", "hom"]:
            return "Somatic"
        
        # If NT field indicates something else but filter is PASS, still return Somatic
        return "Somatic"

    # 2. Germline: Filter failed, but NT indicates variant presence in Normal
    if nt_val in ["het", "hom"]:
        # Rescue if Normal Depth is sufficient (>= 2 per Strelka specs)
        if normal_dp >= 2:
            return "Germline"
        return "Artifact"

    # 3. Reference: Filter failed, but NT indicates Normal is Reference
    if nt_val == "ref":
        # Rescue if Normal Depth is sufficient
        if normal_dp >= 2:
            return "Reference"
        return "Artifact"

    # 4. Artifact: Catch-all for everything else (e.g. NT='conflict', LowDepth with dp<2)
    return "Artifact"


def classify_deepsomatic_variant(filter_val):
    """
    Classifies a DeepSomatic variant into Somatic, Germline, Reference, or Artifact.

    DeepSomatic already provides explicit FILTER labels:
    - PASS: Somatic variant
    - germline: Germline variant
    - reference: Reference site
    - Other FILTER values: Artifacts (low quality, etc.)

    Args:
        filter_val (str): The value from the VCF FILTER column.

    Returns:
        str: One of 'Somatic', 'Germline', 'Reference', or 'Artifact'.
    """
    # Normalize filter value
    if filter_val is None or filter_val == "." or filter_val == "PASS":
        return "Somatic"

    # Check explicit DeepSomatic filter labels
    if "GERMLINE" in filter_val.upper():
        return "Germline"

    if "REFCALL" in filter_val.upper() or "REF_CALL" in filter_val.upper():
        return "Reference"

    # Everything else is an artifact (failed filters)
    return "Artifact"


def classify_mutect2_variant(filter_val):
    """
    Classifies a Mutect2 variant into Somatic, Germline, Reference, or Artifact.

    Note: Mutect2 typically only outputs PASS variants in somatic calling mode.
    This function handles edge cases where other filters might appear.

    Args:
        filter_val (str): The value from the VCF FILTER column.

    Returns:
        str: One of 'Somatic', 'Germline', 'Reference', or 'Artifact'.
    """
    # Mutect2 in somatic mode primarily outputs PASS variants
    if filter_val is None or filter_val == "." or filter_val == "PASS":
        return "Somatic"

    # Check for common Mutect2 germline-related filters (if present)
    filter_upper = filter_val.upper()
    if "GERMLINE" in filter_upper or "NORMAL_ARTIFACT" in filter_upper:
        return "Germline"

    # Other failed filters are artifacts
    return "Artifact"


def classify_variant(variant, caller_name, sample_indices=None):
    """
    Universal variant classifier that dispatches to caller-specific functions.

    Args:
        variant: cyvcf2.Variant object
        caller_name (str): Name of the variant caller ('strelka', 'deepsomatic', 'mutect2')
        sample_indices (dict): Optional dict mapping 'tumor' and 'normal' to sample indices

    Returns:
        str: One of 'Somatic', 'Germline', 'Reference', or 'Artifact'.
    """
    filter_val = variant.FILTER if variant.FILTER else None

    if not caller_name:
        # Fallback for unknown caller
        if filter_val == "PASS":
            return "Somatic"
        else:
            return "Artifact"

    caller_lower = caller_name.lower()

    if caller_lower == "strelka":
        try:
            # Extract Strelka-specific fields
            nt_val = None
            normal_dp = 0

            # Try to get NT field (Normal Type)
            # Note: cyvcf2's 'in' operator doesn't work reliably for INFO fields
            # Use .get() directly and check if result is not None
            try:
                nt_val = variant.INFO.get('NT')
            except:
                nt_val = None

            # Get normal sample depth
            if sample_indices and 'normal' in sample_indices:
                normal_idx = sample_indices['normal']
                dp_array = variant.format('DP')
                if dp_array is not None and len(dp_array) > normal_idx:
                    dp_val = dp_array[normal_idx]
                    if dp_val is not None and len(dp_val) > 0:
                        normal_dp = dp_val[0] if dp_val[0] is not None else 0

            return classify_strelka_variant(filter_val, nt_val, normal_dp)
        except Exception as e:
            print(
                f"Warning: Strelka classification failed for variant at {variant.CHROM}:{variant.POS}: {e}"
            )
            return "Artifact"

    elif caller_lower == "deepsomatic":
        return classify_deepsomatic_variant(filter_val)

    elif caller_lower == "mutect2":
        return classify_mutect2_variant(filter_val)

    else:
        # Unknown caller, use generic PASS/FAIL logic
        if filter_val is None or filter_val == "." or filter_val == "PASS":
            return "Somatic"
        else:
            return "Artifact"


def detect_rna_edit_variant(variant) -> bool:
    """
    Detect if a variant is flagged as RNA editing.
    
    Checks for:
    - FILTER field containing 'RNAedit'
    - INFO fields with REDI_* prefix
    - REDI_SITES or RNA_EDITING_SCORE
    
    Args:
        variant: cyvcf2.Variant object
        
    Returns:
        bool: True if variant appears to be RNA-edited
    """
    try:
        # Check FILTER field
        if variant.FILTER and "RNAedit" in str(variant.FILTER):
            return True
        
        # Check INFO fields for RNA editing markers
        info_keys = variant.INFO.keys() if hasattr(variant.INFO, 'keys') else []
        for key in info_keys:
            if 'REDI' in key.upper() or 'RNA' in key.upper() and 'EDIT' in key.upper():
                return True
                
    except Exception:
        pass
    
    return False


def detect_cosmic_gnomad_annotation(variant) -> bool:
    """
    Detect if a variant has Cosmic/gnomAD annotation.
    
    Checks for:
    - COSMIC_ID in INFO
    - gnomAD_AF or gnomAD_* fields
    - COSMICID field
    
    Args:
        variant: cyvcf2.Variant object
        
    Returns:
        bool: True if variant has cosmic/gnomad annotation
    """
    try:
        info = variant.INFO
        
        # Check for common cosmic/gnomad annotation keys
        cosmic_keys = ['COSMIC_ID', 'COSMICID', 'COSMIC', 'gnomAD_AF', 'gnomAD_AFR_AF']
        for key in cosmic_keys:
            if key in info:
                return True
        
        # Check for any gnomAD-prefixed field
        info_keys = info.keys() if hasattr(info, 'keys') else []
        for key in info_keys:
            if 'gnomAD' in str(key) or 'COSMIC' in str(key):
                return True
                
    except Exception:
        pass
    
    return False


def detect_no_consensus_variant(variant) -> bool:
    """
    Detect if a variant is flagged as NoConsensus.
    
    Checks FILTER field for explicit 'NoConsensus' marker.
    
    Args:
        variant: cyvcf2.Variant object
        
    Returns:
        bool: True if variant is marked as NoConsensus
    """
    try:
        if variant.FILTER and str(variant.FILTER) == "NoConsensus":
            return True
    except Exception:
        pass
    
    return False


def classify_annotated_variant(variant) -> str:
    """
    Classify a variant from an annotated VCF (cosmic_gnomad, rna_editing, filtered).
    
    Assigns categories based on FILTER field and special annotations.
    Returns one of: Somatic, Germline, Reference, Artifact, RNA_Edit, NoConsensus
    
    Priority:
    1. NoConsensus (explicit FILTER)
    2. RNA_Edit (FILTER or INFO markers)
    3. PASS/filter-based classification
    
    Args:
        variant: cyvcf2.Variant object
        
    Returns:
        str: Classification category
    """
    # Check for NoConsensus first
    if detect_no_consensus_variant(variant):
        return "NoConsensus"
    
    # Check for RNA editing
    if detect_rna_edit_variant(variant):
        return "RNA_Edit"
    
    # Fall back to FILTER-based classification
    filter_val = variant.FILTER if variant.FILTER else "PASS"
    
    if filter_val == "PASS" or filter_val is None or filter_val == ".":
        return "Somatic"
    elif "GERMLINE" in str(filter_val).upper():
        return "Germline"
    elif "REFERENCE" in str(filter_val).upper() or "REF" in str(filter_val).upper():
        return "Reference"
    else:
        return "Artifact"


def get_sample_indices(vcf_obj, caller_name):
    """
    Determine tumor and normal sample indices from VCF samples.

    Args:
        vcf_obj: cyvcf2.VCF object
        caller_name (str): Name of the variant caller

    Returns:
        dict: {'tumor': int, 'normal': int} or None
    """
    samples = vcf_obj.samples

    if len(samples) < 2:
        return None

    indices = {}

    # Try to identify tumor and normal samples
    for i, sample in enumerate(samples):
        sample_lower = sample.lower()
        if "tumor" in sample_lower or "tumour" in sample_lower:
            indices["tumor"] = i
        elif "normal" in sample_lower:
            indices["normal"] = i

    # Fallback: assume first is tumor, second is normal
    if "tumor" not in indices and len(samples) >= 1:
        indices["tumor"] = 0
    if "normal" not in indices and len(samples) >= 2:
        indices["normal"] = 0
        indices["tumor"] = 1

    return indices


print("✓ Variant classification functions defined")