#!/usr/bin/env python3
"""
Variant Classification Module

Provides unified classification of variants from different callers into biological categories:
- Somatic: High-confidence somatic variants
- Germline: Variants present in normal sample
- Reference: Reference calls (no variant)
- Artifact: Low-quality or failed filter variants, OR variants with inconsistent classifications across callers (consensus mode)
- NoConsensus: Variants that don't meet consensus thresholds or don't fit other categories (rescue mode)

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
        # Check if this is a consensus caller
        if "consensus" in caller_lower:
            # Consensus VCFs already have biological categories in FILTER field
            # Use the FILTER value directly if it's a valid biological category
            if filter_val in ["Somatic", "Germline", "Reference", "Artifact"]:
                return filter_val
            elif filter_val == "NoConsensus":
                # NoConsensus is a threshold filter, not a classification
                # This shouldn't happen in consensus VCFs that are being rescued
                return "Artifact"
            elif filter_val is None or filter_val == "." or filter_val == "PASS":
                # Legacy PASS value
                return "Somatic"
            else:
                # Unknown filter value
                return "Artifact"
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
            'Description': 'Low quality variant, technical artifact, or inconsistent classifications across callers'
        },
        {
            'ID': 'NoConsensus',
            'Description': 'Does not meet consensus threshold or does not fit other classification categories'
        }
    ]


def compute_unified_classification_consensus(variant_data, snv_threshold, indel_threshold, debug=False):
    """
    Compute unified biological classification for consensus mode.
    
    CORRECTED logic for consensus mode:
    1. Check if variant meets consensus threshold based on caller count
    2. If meets threshold but callers disagree -> Artifact  
    3. If meets threshold and callers agree -> Use agreed classification
    4. If doesn't meet threshold -> NoConsensus
    5. Future: Check quality-based Artifact rules (low AD/DP)
    
    Args:
        variant_data (dict): Aggregated variant data with 'callers', 'filters_normalized', 'is_snv'
        snv_threshold (int): SNV consensus threshold
        indel_threshold (int): Indel consensus threshold
        debug (bool): Enable debug output
    
    Returns:
        str: Unified biological classification
    """
    from collections import Counter
    
    # Get individual caller classifications (exclude consensus callers)
    individual_filters = []
    individual_callers = []
    for i, caller in enumerate(variant_data['callers']):
        if not caller.endswith('_consensus'):
            if i < len(variant_data['filters_normalized']):
                individual_filters.append(variant_data['filters_normalized'][i])
                individual_callers.append(caller)
    
    if debug:
        print(f"DEBUG: Individual callers: {individual_callers}")
        print(f"DEBUG: Individual filters: {individual_filters}")
    
    if not individual_filters:
        if debug:
            print("DEBUG: No individual callers -> NoConsensus")
        return 'NoConsensus'  # No individual callers
    
    # Determine required threshold based on variant type
    required_threshold = snv_threshold if variant_data.get('is_snv', True) else indel_threshold
    caller_count = len(individual_filters)
    
    if debug:
        print(f"DEBUG: Caller count: {caller_count}, Required threshold: {required_threshold}")
    
    # Check if variant meets consensus threshold
    if caller_count < required_threshold:
        if debug:
            print(f"DEBUG: Below threshold -> NoConsensus")
        return 'NoConsensus'
    
    # Check for classification consistency among callers
    unique_classifications = set(individual_filters)
    
    if debug:
        print(f"DEBUG: Unique classifications: {unique_classifications}")
    
    if len(unique_classifications) == 1:
        # All callers agree - use the agreed classification
        agreed_classification = list(unique_classifications)[0]
        
        if debug:
            print(f"DEBUG: All callers agree -> {agreed_classification}")
        
        # Future extension point: Add quality-based Artifact rules here
        # if is_low_quality_artifact(variant_data, agreed_classification):
        #     return 'Artifact'
        
        return agreed_classification
    else:
        # Callers disagree - mark as Artifact due to inconsistency
        if debug:
            print(f"DEBUG: Callers disagree -> Artifact")
        return 'Artifact'


def compute_unified_classification_rescue(variant_data, modality_map):
    """
    Compute unified biological classification for rescue mode.
    
    CORRECTED logic for rescue mode:
    1. First check if we have consensus labels from each modality
    2. If both modalities have consensus -> Apply cross-modality rules
    3. If only one modality has consensus -> Use that consensus
    4. If no consensus labels -> Check individual caller patterns
    5. Single callers or insufficient cross-modality support -> NoConsensus
    6. Future: Add quality-based Artifact rules
    
    Args:
        variant_data (dict): Aggregated variant data with 'callers' and 'filters_normalized'
        modality_map (dict): Maps caller names to modalities ('DNA' or 'RNA')
    
    Returns:
        str: Unified biological classification
    """
    # Step 1: Extract consensus labels from DNA and RNA modalities
    dna_consensus_label = None
    rna_consensus_label = None
    
    for i, caller in enumerate(variant_data['callers']):
        if caller.endswith('_consensus'):
            if 'DNA' in caller.upper():
                dna_consensus_label = variant_data['filters_normalized'][i]
            elif 'RNA' in caller.upper():
                rna_consensus_label = variant_data['filters_normalized'][i]
    
    # Step 2: Apply cross-modality consensus rules
    if dna_consensus_label and rna_consensus_label:
        # Both modalities have consensus
        if dna_consensus_label == rna_consensus_label:
            # Agreement across modalities - use agreed classification
            return dna_consensus_label
        elif dna_consensus_label == 'Artifact' and rna_consensus_label == 'Artifact':
            # Both modalities are Artifact
            return 'Artifact'
        elif dna_consensus_label != 'Artifact' and rna_consensus_label != 'Artifact':
            # Cross-modality disagreement on non-Artifact classifications
            return 'Artifact'
        else:
            # One modality is Artifact, other is not - use DNA priority
            return dna_consensus_label
    
    # Step 3: Single modality consensus available
    elif dna_consensus_label and not rna_consensus_label:
        return dna_consensus_label
    elif rna_consensus_label and not dna_consensus_label:
        return rna_consensus_label
    
    # Step 4: No consensus labels - analyze individual caller patterns
    else:
        # Get individual callers by modality (exclude consensus callers)
        dna_callers = []
        rna_callers = []
        
        for i, caller in enumerate(variant_data['callers']):
            if not caller.endswith('_consensus'):
                modality = modality_map.get(caller, 'UNKNOWN')
                if modality == 'DNA' and i < len(variant_data['filters_normalized']):
                    dna_callers.append(variant_data['filters_normalized'][i])
                elif modality == 'RNA' and i < len(variant_data['filters_normalized']):
                    rna_callers.append(variant_data['filters_normalized'][i])
        
        # Check if we have sufficient cross-modality support
        has_dna_support = len(dna_callers) > 0
        has_rna_support = len(rna_callers) > 0
        
        if not (has_dna_support and has_rna_support):
            # Insufficient cross-modality support
            return 'NoConsensus'
        
        # Check for cross-modality consistency
        dna_classifications = set(dna_callers)
        rna_classifications = set(rna_callers)
        
        # If both modalities have consistent internal classifications
        if len(dna_classifications) == 1 and len(rna_classifications) == 1:
            dna_class = list(dna_classifications)[0]
            rna_class = list(rna_classifications)[0]
            
            if dna_class == rna_class:
                # Cross-modality agreement
                # Future extension point: Add quality-based Artifact rules here
                # if is_low_quality_artifact(variant_data, dna_class):
                #     return 'Artifact'
                return dna_class
            else:
                # Cross-modality disagreement
                return 'Artifact'
        else:
            # Internal inconsistency within modalities
            return 'Artifact'


def is_low_quality_artifact(variant_data, proposed_classification, min_dp=10, min_vaf=0.05):
    """
    Future extension: Check if variant should be classified as Artifact based on quality metrics.
    
    This function provides a framework for quality-based Artifact classification that can be
    enabled in the future. Currently returns False (disabled) but can be extended with:
    - Low allelic depth (AD) thresholds
    - Low total depth (DP) thresholds  
    - Low variant allele frequency (VAF) thresholds
    - Strand bias metrics
    - Mapping quality thresholds
    
    Args:
        variant_data (dict): Aggregated variant data with genotype information
        proposed_classification (str): The classification that would be assigned
        min_dp (int): Minimum depth threshold (default: 10)
        min_vaf (float): Minimum VAF threshold (default: 0.05)
    
    Returns:
        bool: True if variant should be classified as Artifact based on quality
    """
    # Currently disabled - return False to use proposed_classification
    # Future implementation could check:
    
    # # Skip quality checks for Reference calls (expected to have low VAF)
    # if proposed_classification == 'Reference':
    #     return False
    # 
    # # Get aggregated genotype statistics
    # gt_agg = variant_data.get('gt_aggregated', {})
    # 
    # # Check depth thresholds
    # mean_dp = gt_agg.get('dp_mean')
    # if mean_dp is not None and mean_dp < min_dp:
    #     return True
    # 
    # # Check VAF thresholds for non-Reference calls
    # mean_vaf = gt_agg.get('vaf_mean')
    # if mean_vaf is not None and mean_vaf < min_vaf:
    #     return True
    # 
    # # Future: Add more quality checks here
    # # - Strand bias detection
    # # - Mapping quality thresholds
    # # - Allelic depth ratios
    # # - Base quality scores
    
    return False


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
            'Description': 'Variant Classification: Somatic, Germline, Reference, Artifact, or NoConsensus'
        },
        {
            'ID': 'VC_CALLERS',
            'Number': '.',
            'Type': 'String',
            'Description': 'Callers that classified this variant (with classification)'
        }
    ]
