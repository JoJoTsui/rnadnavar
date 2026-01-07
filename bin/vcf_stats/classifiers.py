#!/usr/bin/env python3
"""
Variant Classifier Module

Stage-aware classification system based on VCF FILTER field values and caller-specific logic.
Classifies variants into biological categories: Somatic, Germline, Reference,
Artifact, RNAedit, and NoConsensus.

This module provides:
- Unified FILTER-based classification for consensus/rescue/annotation stages
- Caller-specific classification for normalized stage VCFs (Strelka, DeepSomatic, Mutect2)

The classification logic imports from the production vcf_utils pipeline for consistency.
"""

import sys
from pathlib import Path

# Add bin/common and bin/vcf_utils to path for shared imports
_bin_common_path = Path(__file__).parent.parent / "common"
if str(_bin_common_path) not in sys.path:
    sys.path.insert(0, str(_bin_common_path))

_bin_vcf_utils_path = Path(__file__).parent.parent / "vcf_utils"
if str(_bin_vcf_utils_path) not in sys.path:
    sys.path.insert(0, str(_bin_vcf_utils_path))

# Import from shared modules
from classification import (
    classify_deepsomatic_variant,
    classify_mutect2_variant,
    classify_strelka_variant,
    classify_variant_from_record,
    get_sample_indices,
)
from vcf_config import validate_category

# Re-export imported functions for backward compatibility
__all__ = [
    "classify_variant_from_record",
    "classify_strelka_variant",
    "classify_deepsomatic_variant",
    "classify_mutect2_variant",
    "get_sample_indices",
    "classify_by_filter",
    "classify_by_filter_value",
    "classify_annotated_variant",
    "normalize_category",
    "validate_category",
]

# ============================================================================
# CALLER-SPECIFIC CLASSIFICATION (for normalized stage)
# ============================================================================
# Note: Classification functions are imported from bin/vcf_utils/classification.py
# to ensure consistency with the production pipeline


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


def classify_by_filter(
    variant, stage: str = None, caller_name: str = None, sample_indices: dict = None
) -> str:
    """
    Classify variant with stage-awareness.

    For normalized stage: Uses caller-specific classification from vcf_utils
    For other stages: Uses unified FILTER-based classification

    Args:
        variant: cyvcf2 Variant object
        stage: VCF processing stage (normalized, consensus, rescue, etc.)
        caller_name: Name of variant caller (for normalized stage)
        sample_indices: Optional dict with 'tumor' and 'normal' sample indices

    Returns:
        Category name
    """
    # For normalized stage, use caller-specific classification from vcf_utils
    if stage == "normalized" and caller_name:
        return classify_variant_from_record(variant, caller_name, sample_indices)

    # For all other stages, use unified FILTER-based classification
    return classify_by_filter_value(variant.FILTER)


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


print("✓ Variant classification functions loaded from shared modules (stage-aware)")
