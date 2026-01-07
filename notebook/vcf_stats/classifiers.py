#!/usr/bin/env python3
"""
Variant Classifier Module

Stage-aware classification system based on VCF FILTER field values and caller-specific logic.
Classifies variants into biological categories: Somatic, Germline, Reference,
Artifact, RNAedit, and NoConsensus.

This module provides:
- Unified FILTER-based classification for consensus/rescue/annotation stages
- Caller-specific classification for normalized stage VCFs (Strelka, DeepSomatic, Mutect2)

The classification logic mirrors the production vcf_utils pipeline for consistency.
"""

import sys
from pathlib import Path

# Add bin/common to path for shared config imports
_bin_common_path = Path(__file__).parent.parent.parent / "bin" / "common"
if str(_bin_common_path) not in sys.path:
    sys.path.insert(0, str(_bin_common_path))

try:
    from vcf_config import CATEGORY_ORDER, FILTER_FIELD_VALUES
except ImportError:
    CATEGORY_ORDER = [
        "Somatic",
        "Germline",
        "Reference",
        "Artifact",
        "RNAedit",
        "NoConsensus",
    ]
    FILTER_FIELD_VALUES = [
        "PASS",
        "Somatic",
        "Germline",
        "Reference",
        "Artifact",
        "NoConsensus",
        "RNAedit",
    ]


# ============================================================================
# CALLER-SPECIFIC CLASSIFICATION (for normalized stage)
# ============================================================================


def classify_strelka_variant(filter_val: str, info_dict: dict = None) -> str:
    """
    Classify Strelka variant using NT field and filter values.

    Strelka Classification Rules:
    - PASS filter -> Somatic
    - Failed filter + NT=het/hom + normal_dp>=2 -> Germline
    - Failed filter + NT=het/hom + normal_dp<2 -> Artifact
    - Failed filter + NT=ref + normal_dp>=2 -> Reference
    - Failed filter + NT=ref + normal_dp<2 -> Artifact
    - All other cases -> Artifact

    Args:
        filter_val: FILTER field value
        info_dict: INFO fields dictionary (for NT field)

    Returns:
        Category name
    """
    if filter_val in (None, ".", "PASS"):
        return "Somatic"

    if info_dict is None:
        return "Artifact"

    nt_val = info_dict.get("NT", "").lower()
    normal_dp = info_dict.get("NORMAL_DP", info_dict.get("DP_NORMAL", 0))

    try:
        normal_dp = int(normal_dp) if normal_dp else 0
    except (ValueError, TypeError):
        normal_dp = 0

    if nt_val in ("het", "hom"):
        return "Germline" if normal_dp >= 2 else "Artifact"
    elif nt_val == "ref":
        return "Reference" if normal_dp >= 2 else "Artifact"

    return "Artifact"


def classify_deepsomatic_variant(filter_val: str, info_dict: dict = None) -> str:
    """
    Classify DeepSomatic variant using filter values.

    DeepSomatic Classification Rules:
    - PASS/None/. -> Somatic
    - Contains "GERMLINE" -> Germline
    - Contains "REFCALL" -> Reference
    - Any other filter -> Artifact

    Args:
        filter_val: FILTER field value
        info_dict: INFO fields dictionary (not used, for API consistency)

    Returns:
        Category name
    """
    if filter_val in (None, ".", "PASS"):
        return "Somatic"

    filter_str = str(filter_val).upper()

    if "GERMLINE" in filter_str:
        return "Germline"
    if "REFCALL" in filter_str or "REF_CALL" in filter_str:
        return "Reference"

    return "Artifact"


def classify_mutect2_variant(filter_val: str, info_dict: dict = None) -> str:
    """
    Classify Mutect2 variant using filter values.

    Mutect2 Classification Rules:
    - PASS/None/. -> Somatic
    - Contains "GERMLINE" or "NORMAL_ARTIFACT" -> Germline
    - Any other filter -> Artifact

    Args:
        filter_val: FILTER field value
        info_dict: INFO fields dictionary (not used, for API consistency)

    Returns:
        Category name
    """
    if filter_val in (None, ".", "PASS"):
        return "Somatic"

    filter_str = str(filter_val).upper()

    if "GERMLINE" in filter_str or "NORMAL_ARTIFACT" in filter_str:
        return "Germline"

    return "Artifact"


# Caller classification dispatch table
CALLER_CLASSIFIERS = {
    "strelka": classify_strelka_variant,
    "deepsomatic": classify_deepsomatic_variant,
    "mutect2": classify_mutect2_variant,
}


def classify_caller_variant(
    caller_name: str, filter_val: str, info_dict: dict = None
) -> str:
    """
    Classify variant using caller-specific logic for normalized stage.

    Args:
        caller_name: Name of variant caller (strelka, deepsomatic, mutect2)
        filter_val: FILTER field value
        info_dict: INFO fields dictionary

    Returns:
        Category name
    """
    caller_key = caller_name.lower() if caller_name else ""
    classifier = CALLER_CLASSIFIERS.get(caller_key)

    if classifier:
        return classifier(filter_val, info_dict)

    # Fallback to unified classification
    return classify_by_filter_value(filter_val)


# ============================================================================
# UNIFIED FILTER-BASED CLASSIFICATION (for consensus/rescue/annotation stages)
# ============================================================================


def classify_by_filter_value(filter_val) -> str:
    """
    Classify variant based solely on FILTER field value.

    Used for consensus, rescue, and annotation stage VCFs where the FILTER
    field already contains the biological category.

    Args:
        filter_val: FILTER field value (string or None)

    Returns:
        Category name: 'Somatic', 'Germline', 'Reference', 'Artifact',
                      'RNAedit', or 'NoConsensus'
    """
    # Handle PASS or missing filter (None, ".")
    if filter_val in (None, ".", "PASS"):
        return "Somatic"

    # Convert to lowercase string for case-insensitive matching
    filter_str = str(filter_val).lower()

    # Check for specific biological category values
    if filter_str == "somatic":
        return "Somatic"
    if filter_str == "germline":
        return "Germline"
    if filter_str == "reference":
        return "Reference"
    if filter_str == "rnaedit":
        return "RNAedit"
    if filter_str == "noconsensus":
        return "NoConsensus"
    if filter_str == "artifact":
        return "Artifact"

    # Check for partial matches (legacy/alternative formats)
    if "germline" in filter_str:
        return "Germline"
    if "reference" in filter_str:
        return "Reference"
    if "rnaedit" in filter_str or "rna_edit" in filter_str:
        return "RNAedit"
    if "noconsensus" in filter_str or "no_consensus" in filter_str:
        return "NoConsensus"

    # Default to Artifact for any other filter value
    return "Artifact"


def classify_by_filter(variant, stage: str = None, caller_name: str = None) -> str:
    """
    Classify variant with stage-awareness.

    For normalized stage: Uses caller-specific classification
    For other stages: Uses unified FILTER-based classification

    Args:
        variant: cyvcf2 Variant object
        stage: VCF processing stage (normalized, consensus, rescue, etc.)
        caller_name: Name of variant caller (for normalized stage)

    Returns:
        Category name
    """
    filter_val = variant.FILTER

    # For normalized stage, use caller-specific classification
    if stage == "normalized" and caller_name:
        try:
            info_dict = dict(variant.INFO) if hasattr(variant, "INFO") else {}
        except:
            info_dict = {}
        return classify_caller_variant(caller_name, filter_val, info_dict)

    # For all other stages, use unified FILTER-based classification
    return classify_by_filter_value(filter_val)


def classify_annotated_variant(variant) -> str:
    """
    Classify variants using unified FILTER-based classification.

    Uses the classify_by_filter_value() function for all classification.
    This provides consistent behavior across all VCF types and stages.

    Returns one of: Somatic, Germline, Reference, Artifact, RNAedit, NoConsensus

    Args:
        variant: cyvcf2 Variant object

    Returns:
        Category name based on FILTER value
    """
    return classify_by_filter_value(variant.FILTER)


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


def validate_category(category: str) -> bool:
    """
    Validate that a category is in the defined CATEGORY_ORDER.

    Args:
        category: Category name to validate

    Returns:
        True if valid, False otherwise
    """
    return category in CATEGORY_ORDER


def normalize_category(category: str) -> str:
    """
    Normalize category name to standard form.

    Handles variations like RNA_Edit -> RNAedit.

    Args:
        category: Category name (possibly non-standard)

    Returns:
        Normalized category name
    """
    if not category:
        return "Artifact"

    cat_lower = category.lower().replace("_", "")

    if cat_lower == "somatic":
        return "Somatic"
    elif cat_lower == "germline":
        return "Germline"
    elif cat_lower == "reference":
        return "Reference"
    elif cat_lower == "artifact":
        return "Artifact"
    elif cat_lower in ("rnaedit", "rna_edit"):
        return "RNAedit"
    elif cat_lower in ("noconsensus", "no_consensus"):
        return "NoConsensus"

    return category  # Return as-is if unknown


print("✓ Variant classification functions defined (stage-aware)")
