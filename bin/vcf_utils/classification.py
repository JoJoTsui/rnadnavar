#!/usr/bin/env python3
"""
Variant Classification Module

Provides unified classification of variants from different callers into biological categories:
- Somatic: High-confidence somatic variants
- Germline: Variants present in normal sample
- Reference: Reference calls (no variant)
- Artifact: Low-quality or failed filter variants

This module standardizes variant classification across Strelka, DeepSomatic, and Mutect2.
"""


def classify_strelka_variant(filter_val, nt_val, normal_dp):
    """
    Classifies a Strelka variant into Somatic, Germline, Reference, or Artifact.
    
    Strelka uses:
    - FILTER: PASS indicates high-confidence somatic variant
    - INFO/NT: Normal genotype (ref, het, hom)
    - Normal DP: Read depth in normal sample
    
    Args:
        filter_val (str): The value from the VCF FILTER column (e.g., 'PASS', 'LowDepth').
                          Can be None if the VCF parser returns None for PASS.
        nt_val (str): The value of the INFO/NT field (e.g., 'ref', 'het', 'hom').
        normal_dp (int): The read depth of the Normal sample (FORMAT/DP).
        
    Returns:
        str: One of ['Somatic', 'Germline', 'Reference', 'Artifact']
    """
    
    # 1. Somatic: Strictly relies on PASS filter
    # Some parsers return None for PASS, Strelka writes "PASS" string.
    if filter_val == "PASS" or filter_val is None or filter_val == ".":
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
    - GERMLINE: Germline variant
    - RefCall: Reference (no variant)
    - Other filters: Artifact
    
    Args:
        filter_val (str): The value from the VCF FILTER column.
                          Can be None (treated as PASS by some parsers).
        
    Returns:
        str: One of ['Somatic', 'Germline', 'Reference', 'Artifact']
    """
    
    # Normalize filter value
    if filter_val is None or filter_val == "." or filter_val == "PASS":
        return "Somatic"
    
    # Check explicit DeepSomatic filter labels (case-insensitive)
    filter_upper = str(filter_val).upper()
    
    if "GERMLINE" in filter_upper:
        return "Germline"
    
    if "REFCALL" in filter_upper or "REF_CALL" in filter_upper:
        return "Reference"
    
    # Everything else is an artifact (failed filters)
    return "Artifact"


def classify_mutect2_variant(filter_val, info_dict=None):
    """
    Classifies a Mutect2 variant into Somatic, Germline, Reference, or Artifact.
    
    Note: Mutect2 typically only outputs PASS variants in somatic calling mode.
    This function handles edge cases where other filters might appear.
    
    Args:
        filter_val (str): The value from the VCF FILTER column.
                          Can be None (treated as PASS by some parsers).
        info_dict (dict): Optional INFO field dictionary to check additional flags.
        
    Returns:
        str: One of ['Somatic', 'Germline', 'Reference', 'Artifact']
    """
    
    # Mutect2 in somatic mode primarily outputs PASS variants
    if filter_val is None or filter_val == "." or filter_val == "PASS":
        return "Somatic"
    
    # Check for common Mutect2 germline-related filters (if present)
    filter_upper = str(filter_val).upper()
    
    if "GERMLINE" in filter_upper or "NORMAL_ARTIFACT" in filter_upper:
        return "Germline"
    
    # Other failed filters are artifacts
    return "Artifact"


def classify_variant_from_record(variant, caller_name, sample_indices=None):
    """
    Universal variant classifier that dispatches to caller-specific functions.
    Works with cyvcf2.Variant objects.
    
    Args:
        variant: cyvcf2.Variant object
        caller_name (str): Name of the variant caller ('strelka', 'deepsomatic', 'mutect2')
        sample_indices (dict): Optional dict mapping 'tumor' and 'normal' to sample indices
        
    Returns:
        str: One of ['Somatic', 'Germline', 'Reference', 'Artifact']
    """
    
    filter_val = variant.FILTER
    caller_lower = caller_name.lower()
    
    if caller_lower == "strelka":
        # Strelka requires NT field and normal depth
        try:
            nt_val = variant.INFO.get("NT", "")
            
            # Get normal sample depth
            normal_dp = 0
            if sample_indices and "normal" in sample_indices:
                normal_idx = sample_indices["normal"]
                dp_array = variant.format("DP")
                if dp_array is not None and len(dp_array) > normal_idx:
                    normal_dp = dp_array[normal_idx][0] if dp_array[normal_idx][0] is not None else 0
            
            return classify_strelka_variant(filter_val, nt_val, normal_dp)
        except Exception as e:
            # If classification fails, default to Artifact
            return "Artifact"
    
    elif caller_lower == "deepsomatic":
        return classify_deepsomatic_variant(filter_val)
    
    elif caller_lower == "mutect2":
        return classify_mutect2_variant(filter_val)
    
    else:
        # Unknown caller, use generic PASS/FAIL logic
        if filter_val is None or filter_val == "." or filter_val == "PASS":
            return "Somatic"
        return "Artifact"


def classify_variant_from_dict(variant_dict, caller_name):
    """
    Classify a variant from a dictionary representation (used in aggregation).
    
    Args:
        variant_dict (dict): Dictionary containing variant info with keys:
                            - 'filter': FILTER value
                            - 'info': INFO dict (for Strelka NT field)
                            - 'normal_dp': Normal depth (for Strelka)
        caller_name (str): Name of the variant caller
        
    Returns:
        str: One of ['Somatic', 'Germline', 'Reference', 'Artifact']
    """
    
    filter_val = variant_dict.get('filter')
    caller_lower = caller_name.lower()
    
    if caller_lower == "strelka":
        nt_val = variant_dict.get('info', {}).get('NT', '')
        normal_dp = variant_dict.get('normal_dp', 0)
        return classify_strelka_variant(filter_val, nt_val, normal_dp)
    
    elif caller_lower == "deepsomatic":
        return classify_deepsomatic_variant(filter_val)
    
    elif caller_lower == "mutect2":
        info_dict = variant_dict.get('info', {})
        return classify_mutect2_variant(filter_val, info_dict)
    
    else:
        # Unknown caller
        if filter_val is None or filter_val == "." or filter_val == "PASS":
            return "Somatic"
        return "Artifact"


def get_sample_indices(vcf_obj, caller_name=None):
    """
    Determine tumor and normal sample indices from VCF samples.
    
    Args:
        vcf_obj: cyvcf2.VCF object
        caller_name (str): Optional name of the variant caller
        
    Returns:
        dict: {'tumor': int, 'normal': int} or None if not paired
    """
    samples = vcf_obj.samples
    
    if len(samples) < 2:
        return None
    
    indices = {}
    
    # Try to identify tumor and normal samples by name
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
        indices["normal"] = 1
    
    return indices


def normalize_filter_value(classification):
    """
    Return the biological classification as the normalized filter value.
    
    Args:
        classification (str): One of ['Somatic', 'Germline', 'Reference', 'Artifact']
        
    Returns:
        str: The biological classification itself (Somatic/Germline/Reference/Artifact)
    """
    # Return classification as-is - no conversion needed
    # Valid values: Somatic, Germline, Reference, Artifact
    return classification


def get_unified_filter_headers():
    """
    Get biological category FILTER header lines for unified VCF output.
    
    Returns:
        list: List of FILTER header line dictionaries
    """
    return [
        {
            'ID': 'Somatic',
            'Description': 'High-confidence somatic variant specific to tumor'
        },
        {
            'ID': 'Germline',
            'Description': 'Germline variant detected in normal sample'
        },
        {
            'ID': 'Reference',
            'Description': 'Reference call - no variant detected'
        },
        {
            'ID': 'Artifact',
            'Description': 'Low quality variant or technical artifact'
        }
    ]


def get_classification_info_headers():
    """
    Get INFO header lines for classification metadata.
    
    Returns:
        list: List of INFO header line dictionaries
    """
    return [
        {
            'ID': 'VC',
            'Number': '1',
            'Type': 'String',
            'Description': 'Variant Classification: Somatic, Germline, Reference, or Artifact'
        },
        {
            'ID': 'VC_CALLERS',
            'Number': '.',
            'Type': 'String',
            'Description': 'Callers that classified this variant (with classification)'
        }
    ]
