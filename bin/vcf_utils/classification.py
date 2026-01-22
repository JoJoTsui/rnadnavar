#!/usr/bin/env python3
"""
Variant Classification Module

Provides unified classification of variants from different callers into biological categories:
- Somatic: High-confidence somatic variants
- Germline: Variants present in normal sample
- Reference: Reference calls (no variant)
- Artifact: Low-quality or failed filter variants, OR variants with inconsistent classifications across callers (consensus mode)
- NoConsensus: Variants that don't meet consensus thresholds or don't fit other categories (rescue mode)
- RNAedit: RNA editing events (added via annotation)

This module standardizes variant classification across Strelka, DeepSomatic, and Mutect2.
Now integrated with the unified classifier for consistency across all pipeline stages.
"""

# =============================================================================
# Mutect2 FILTER Classification Mapping
# =============================================================================
# Maps Mutect2 FILTER labels to biological categories.
# Used for reference - actual classification uses priority-based logic in
# classify_mutect2_variant() function.
#
# Priority Hierarchy:
# 1. Somatic   (FILTER == 'PASS', None, or '.' - unfiltered/passing variants)
# 2. Germline  (Contains 'germline' or 'haplotype')
# 3. Reference (Contains 'panel_of_normals', 'contamination', 'possible_numt')
# 4. Artifact  (Everything else - technical filter failures)
#
# Note: Mutect2 VCFs can be unfiltered (FILTER='.') or filtered (FILTER='PASS'/labels)
# Unfiltered variants are treated as Somatic candidates.

MUTECT2_FILTER_CLASSIFICATION = {
    # Somatic (PASS or unfiltered)
    "PASS": "Somatic",
    ".": "Somatic",  # Unfiltered variant
    
    # Germline (germline-related filters)
    "germline": "Germline",
    "haplotype": "Germline",
    
    # Reference (systematic/external reference issues)
    "panel_of_normals": "Reference",
    "contamination": "Reference",
    "possible_numt": "Reference",
    
    # Artifact (all other technical filters)
    "FAIL": "Artifact",
    "base_qual": "Artifact",
    "clustered_events": "Artifact",
    "duplicate": "Artifact",
    "fragment": "Artifact",
    "low_allele_frac": "Artifact",
    "map_qual": "Artifact",
    "multiallelic": "Artifact",
    "n_ratio": "Artifact",
    "normal_artifact": "Artifact",
    "orientation": "Artifact",
    "position": "Artifact",
    "slippage": "Artifact",
    "strand_bias": "Artifact",
    "strict_strand": "Artifact",
    "weak_evidence": "Artifact",
}

# Priority order for multi-filter classification (highest to lowest)
# When a variant has multiple FILTER values, the classification with the
# highest priority (lowest index) takes precedence.
# Note: For multi-filter classification, the priority is:
#   Germline > Reference > Artifact (Somatic only applies when PASS is alone)
CLASSIFICATION_PRIORITY = ["Somatic", "Germline", "Reference", "Artifact"]

# =============================================================================

# Lazy import for unified classifier to avoid circular imports
_unified_classifier = None


def _get_unified_classifier():
    """Get or create the unified classifier instance (lazy initialization)."""
    global _unified_classifier
    if _unified_classifier is None:
        from .variant_classifier_unified import UnifiedVariantClassifier

        _unified_classifier = UnifiedVariantClassifier()
    return _unified_classifier


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
    Classifies a Mutect2 VCF FILTER string into exactly ONE of 4 categories:
    Somatic, Germline, Reference, Artifact.

    Priority Hierarchy:
    1. Somatic   (FILTER == 'PASS', None, or '.' - unfiltered/passing variants)
    2. Germline  (Contains 'germline' or 'haplotype')
    3. Reference (Contains 'panel_of_normals', 'contamination', 'possible_numt')
    4. Artifact  (Everything else - technical filter failures)

    Note: Mutect2 VCFs can be:
    - Unfiltered (before FilterMutectCalls): FILTER = '.' or None → Somatic candidates
    - Filtered (after FilterMutectCalls): FILTER = 'PASS' or specific filter labels

    Args:
        filter_val (str): The value from the VCF FILTER column.
                          Can be None (unfiltered variant - treated as Somatic).
                          Can be semicolon-separated for multiple filters.
        info_dict (dict): Optional INFO field dictionary (unused, kept for API compatibility).

    Returns:
        str: One of ['Somatic', 'Germline', 'Reference', 'Artifact']
    """
    # 1. Handle None, empty, or '.' as Somatic (unfiltered/passing variants)
    # cyvcf2 returns None for FILTER when VCF has '.' (unfiltered)
    # These are somatic candidates that haven't been filtered out
    if filter_val is None or filter_val == "" or filter_val == ".":
        return "Somatic"
    
    # Convert to string and normalize
    filter_str = str(filter_val).strip()
    
    # Handle empty string after stripping
    if not filter_str or filter_str == ".":
        return "Somatic"
    
    # 2. Explicit PASS → Somatic
    if filter_str == "PASS":
        return "Somatic"
    
    # 3. Parse Filters
    # Split by semicolon to handle multiple filters (e.g., "map_qual;germline")
    filters = [f.strip().lower() for f in filter_str.split(";") if f.strip()]
    
    # If no valid filters after parsing, treat as Somatic
    if not filters:
        return "Somatic"
    
    # 4. Priority Check: Germline
    # If any filter is germline-related, the biological classification is Germline.
    germline_triggers = {"germline", "haplotype"}
    if any(f in germline_triggers for f in filters):
        return "Germline"
    
    # 5. Priority Check: Reference
    # If not germline, check for systematic/external reference issues.
    reference_triggers = {"panel_of_normals", "contamination", "possible_numt"}
    if any(f in reference_triggers for f in filters):
        return "Reference"
    
    # 6. Default: Artifact
    # Falls here if:
    # - It has explicit filter labels (not PASS/None/.)
    # - It has no Germline tags
    # - It has no Reference tags
    # - Contains technical errors (strand_bias, weak_evidence, base_qual, etc.)
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
                    normal_dp = (
                        dp_array[normal_idx][0]
                        if dp_array[normal_idx][0] is not None
                        else 0
                    )

            return classify_strelka_variant(filter_val, nt_val, normal_dp)
        except Exception:
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

    filter_val = variant_dict.get("filter")
    caller_lower = caller_name.lower()

    if caller_lower == "strelka":
        nt_val = variant_dict.get("info", {}).get("NT", "")
        normal_dp = variant_dict.get("normal_dp", 0)
        return classify_strelka_variant(filter_val, nt_val, normal_dp)

    elif caller_lower == "deepsomatic":
        return classify_deepsomatic_variant(filter_val)

    elif caller_lower == "mutect2":
        info_dict = variant_dict.get("info", {})
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

    Includes all classification categories including RNAedit for RNA editing events.

    Returns:
        list: List of FILTER header line dictionaries
    """
    from .classification_config import FILTER_DESCRIPTIONS

    return [
        {
            "ID": "Somatic",
            "Description": FILTER_DESCRIPTIONS.get(
                "Somatic", "High-confidence somatic variant specific to tumor"
            ),
        },
        {
            "ID": "Germline",
            "Description": FILTER_DESCRIPTIONS.get(
                "Germline", "Germline variant detected in normal sample"
            ),
        },
        {
            "ID": "Reference",
            "Description": FILTER_DESCRIPTIONS.get(
                "Reference", "Reference call - no variant detected"
            ),
        },
        {
            "ID": "Artifact",
            "Description": FILTER_DESCRIPTIONS.get(
                "Artifact",
                "Low quality variant, technical artifact, or inconsistent classifications across callers",
            ),
        },
        {
            "ID": "NoConsensus",
            "Description": FILTER_DESCRIPTIONS.get(
                "NoConsensus",
                "Does not meet consensus threshold or does not fit other classification categories",
            ),
        },
        {
            "ID": "RNAedit",
            "Description": FILTER_DESCRIPTIONS.get(
                "RNAedit",
                "RNA editing event confirmed by REDIportal database and RNA evidence",
            ),
        },
    ]


def compute_unified_classification_consensus(
    variant_data, snv_threshold, indel_threshold
):
    """
    Compute unified biological classification for consensus mode.

    This function now delegates to the unified classifier to ensure consistency
    across all classification stages in the pipeline.

    CORRECTED logic for consensus mode:
    1. Check if total caller count meets threshold
    2. If not enough callers -> NoConsensus
    3. If enough callers -> Use majority vote among classifications
    4. If clear majority -> Use majority classification
    5. If tie -> Artifact (disagreement indicates inconsistency)

    Args:
        variant_data (dict): Aggregated variant data with 'callers', 'filters_normalized', 'is_snv'
        snv_threshold (int): SNV consensus threshold
        indel_threshold (int): Indel consensus threshold

    Returns:
        str: Unified biological classification
    """
    # Create classifier with custom thresholds (lazy import to avoid circular refs)
    from .variant_classifier_unified import UnifiedVariantClassifier

    config = {
        "consensus_snv_threshold": snv_threshold,
        "consensus_indel_threshold": indel_threshold,
    }
    classifier = UnifiedVariantClassifier(config)

    # Delegate to unified classifier
    return classifier.classify_consensus_variant(variant_data)


def compute_unified_classification_rescue(variant_data, modality_map):
    """
    Compute unified biological classification for rescue mode.

    This function now delegates to the unified classifier to ensure consistency
    across all classification stages in the pipeline.

    COMPLETE rescue logic handling both consensus and individual callers:
    1. Separate consensus callers from individual callers
    2. If both DNA and RNA consensus available -> Apply cross-modality consensus rules
    3. If only one modality has consensus -> Check for cross-modality support from individual callers
    4. If no consensus labels -> Analyze individual caller cross-modality patterns
    5. Single modality or insufficient cross-modality support -> NoConsensus
    6. Future: Add quality-based Artifact rules

    Args:
        variant_data (dict): Complete variant data with both consensus and individual callers
        modality_map (dict): Maps caller names to modalities ('DNA' or 'RNA')

    Returns:
        str: Unified biological classification
    """
    # Delegate to the global unified classifier instance (lazy initialization)
    return _get_unified_classifier().classify_rescue_variant(variant_data, modality_map)


def is_low_quality_artifact(
    variant_data, proposed_classification, min_dp=10, min_vaf=0.05
):
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
            "ID": "VC",
            "Number": "1",
            "Type": "String",
            "Description": "Variant Classification: Somatic, Germline, Reference, Artifact, or NoConsensus",
        },
        {
            "ID": "VC_CALLERS",
            "Number": ".",
            "Type": "String",
            "Description": "Callers that classified this variant (with classification)",
        },
    ]
