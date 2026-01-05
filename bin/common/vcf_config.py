#!/usr/bin/env python3
"""
Shared VCF Configuration Module

Central configuration for VCF processing and analysis used by both:
- vcf_utils (Nextflow-based VCF processing pipeline)
- vcf_stats (Jupyter-based VCF statistics analysis)

This module ensures consistent handling of variant categories, colors, and
pipeline stages across the entire workflow.

Usage:
    # In vcf_utils (Nextflow scripts):
    from bin.common.vcf_config import CATEGORY_ORDER, CATEGORY_COLORS

    # In vcf_stats (Jupyter notebooks):
    from vcf_config import CATEGORY_ORDER, CATEGORY_COLORS
"""

# ============================================================================
# VARIANT CLASSIFICATION SYSTEM
# ============================================================================
#
# These categories represent the 6 possible classifications for variants:
# 1. Somatic: Variants detected in tumor/cancer sample only
# 2. Germline: Variants detected in normal/germline sample
# 3. Reference: Variants matching reference genome (homozygous ref)
# 4. Artifact: Low-quality variants, sequencing errors, or technical artifacts
# 5. RNA_Edit: RNA editing sites (A-to-I or C-to-U modifications)
# 6. NoConsensus: Variants where DNA and RNA callers disagree
#

CATEGORY_ORDER = [
    "Somatic",
    "Germline",
    "Reference",
    "Artifact",
    "RNA_Edit",
    "NoConsensus",
]

CATEGORY_COLORS = {
    "Somatic": "#636EFA",  # Blue
    "Germline": "#00CC96",  # Green
    "Reference": "#FFA15A",  # Orange
    "Artifact": "#EF553B",  # Red
    "RNA_Edit": "#AB63FA",  # Purple
    "NoConsensus": "#8A8A8A",  # Gray
}


# ============================================================================
# VCF PROCESSING PIPELINE STAGES
# ============================================================================
#
# These define the complete VCF processing pipeline flow:
# 1. normalized: Individual caller VCFs (Strelka, DeepSomatic, Mutect2)
# 2. consensus: Consensus VCFs (DNA and RNA consensus combined)
# 3. rescue: Combined DNA+RNA variants (rescue analysis)
# 4. cosmic_gnomad: After COSMIC/GnomAD annotation
# 5. rna_editing: After RNA editing detection
# 6. filtered_rescue: Final filtered variants
#

VCF_STAGE_ORDER = [
    "normalized",
    "consensus",
    "rescue",
    "cosmic_gnomad",
    "rna_editing",
    "filtered_rescue",
]

STAGE_DISPLAY_NAMES = {
    "normalized": "Normalized",
    "consensus": "Consensus",
    "rescue": "Rescued",
    "cosmic_gnomad": "COSMIC/GnomAD",
    "rna_editing": "RNA Editing",
    "filtered_rescue": "Filtered",
    # Backward compatibility aliases
    "dna_consensus": "DNA Consensus",
    "rna_consensus": "RNA Consensus",
}


# ============================================================================
# VARIANT CALLER DEFINITIONS
# ============================================================================
#
# Tools used for variant calling in the pipeline
#

TOOLS = ["Strelka", "DeepSomatic", "Mutect2"]
MODALITIES = ["DNA_TUMOR", "DNA_NORMAL", "RNA_TUMOR"]


# ============================================================================
# TIERING SYSTEM FOR RESCUE VARIANTS
# ============================================================================
#
# Rescue variants are tiered based on the number of DNA and RNA callers
# supporting each variant. Tier 1 (T1) is the highest confidence,
# Tier 7 (T7) is the lowest confidence.
#
# Tier Definition Rules:
#   T1: 2+ DNA callers + 2+ RNA callers (strongest consensus)
#   T2: 2+ DNA callers + 1 RNA caller
#   T3: 2+ DNA callers (DNA consensus)
#   T4: 1 DNA caller + 1+ RNA callers
#   T5: 1 DNA caller only
#   T6: 2+ RNA callers (RNA consensus)
#   T7: 1 RNA caller (weakest support)
#

TIER_DEFINITIONS = {
    "T1": {"dna_min": 2, "rna_min": 2, "description": "2+ DNA + 2+ RNA consensus"},
    "T2": {
        "dna_min": 2,
        "rna_min": 1,
        "dna_max": 2,
        "rna_max": 1,
        "description": "2+ DNA + 1 RNA",
    },
    "T3": {"dna_min": 2, "rna_min": 0, "rna_max": 0, "description": "2+ DNA only"},
    "T4": {"dna_min": 1, "dna_max": 1, "rna_min": 1, "description": "1 DNA + 1+ RNA"},
    "T5": {
        "dna_min": 1,
        "dna_max": 1,
        "rna_min": 0,
        "rna_max": 0,
        "description": "1 DNA only",
    },
    "T6": {"dna_min": 0, "dna_max": 0, "rna_min": 2, "description": "2+ RNA only"},
    "T7": {
        "dna_min": 0,
        "dna_max": 0,
        "rna_min": 1,
        "rna_max": 1,
        "description": "1 RNA only",
    },
}

TIER_ORDER = ["T1", "T2", "T3", "T4", "T5", "T6", "T7"]


# ============================================================================
# VCFUTILS / VCF_STATS INTEGRATION SPECIFICATION
# ============================================================================
#
# This section documents how vcf_utils and vcf_stats work together:
#
# vcf_utils (WRITES):
#   - Creates VCF files with FILTER field set to one of CATEGORY_ORDER values
#   - Populates INFO fields: N_DNA_CALLERS_SUPPORT, N_RNA_CALLERS_SUPPORT
#   - Formats: "Somatic", "Germline", "Reference", "Artifact", "NoConsensus"
#   - Special handling for RNA_Edit (may set as FILTER or in custom INFO field)
#
# vcf_stats (READS):
#   - Reads FILTER field and maps to CATEGORY_ORDER values
#   - Parses N_DNA_CALLERS_SUPPORT and N_RNA_CALLERS_SUPPORT for tiering
#   - Performs per-category statistics and visualization
#   - Generates reports with CATEGORY_COLORS for consistent styling
#
# Expected FILTER Field Format in VCF:
#   - One of: "PASS", "Somatic", "Germline", "Reference", "Artifact", "NoConsensus"
#   - RNA_Edit may appear as: FILTER="RNA_Edit" or custom annotation
#
# Expected INFO Fields for Tiering:
#   - N_DNA_CALLERS_SUPPORT=<int>  # Number of DNA callers supporting variant
#   - N_RNA_CALLERS_SUPPORT=<int>  # Number of RNA callers supporting variant
#

FILTER_FIELD_VALUES = [
    "PASS",
    "Somatic",
    "Germline",
    "Reference",
    "Artifact",
    "NoConsensus",
    "RNA_Edit",
]
INFO_FIELD_NAMES = ["N_DNA_CALLERS_SUPPORT", "N_RNA_CALLERS_SUPPORT"]


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================


def get_category_color(category: str) -> str:
    """
    Get the color code for a variant category.

    Args:
        category: One of CATEGORY_ORDER values

    Returns:
        Hex color code or gray (#8A8A8A) if category not found
    """
    return CATEGORY_COLORS.get(category, "#8A8A8A")


def get_stage_display_name(stage: str) -> str:
    """
    Get the human-readable display name for a VCF processing stage.

    Args:
        stage: One of VCF_STAGE_ORDER values

    Returns:
        Display name or stage title if not found
    """
    return STAGE_DISPLAY_NAMES.get(stage, stage.replace("_", " ").title())


def is_annotation_stage(stage: str) -> bool:
    """
    Check if a stage is part of the annotation pipeline (not normalized/consensus).

    Args:
        stage: VCF stage name

    Returns:
        True if stage is in annotation pipeline
    """
    annotation_stages = ["rescue", "cosmic_gnomad", "rna_editing", "filtered_rescue"]
    return stage in annotation_stages


def validate_category(category: str) -> bool:
    """
    Validate that a category is in the defined CATEGORY_ORDER.

    Args:
        category: Category name to validate

    Returns:
        True if valid, False otherwise
    """
    return category in CATEGORY_ORDER


def validate_stage(stage: str) -> bool:
    """
    Validate that a stage is in the defined VCF_STAGE_ORDER.

    Args:
        stage: Stage name to validate

    Returns:
        True if valid, False otherwise
    """
    return stage in VCF_STAGE_ORDER


# ============================================================================
# MIGRATION GUIDE
# ============================================================================
#
# To migrate existing code to use shared config:
#
# Before (local imports):
#     from vcf_stats import CATEGORY_ORDER, CATEGORY_COLORS
#
# After (shared imports):
#     # For notebook code:
#     import sys
#     sys.path.insert(0, "../../bin/common")
#     from vcf_config import CATEGORY_ORDER, CATEGORY_COLORS
#
#     # For Nextflow scripts:
#     from bin.common.vcf_config import CATEGORY_ORDER, CATEGORY_COLORS
#
# Graceful fallback pattern (recommended):
#     try:
#         from vcf_config import CATEGORY_ORDER, CATEGORY_COLORS
#     except ImportError:
#         # Fall back to local definitions if shared config unavailable
#         CATEGORY_ORDER = ["Somatic", "Germline", "Reference", "Artifact", "RNA_Edit", "NoConsensus"]
#         CATEGORY_COLORS = {...}
#

if __name__ == "__main__":
    # Print configuration summary when run directly
    print("=" * 70)
    print("VCF CONFIGURATION SUMMARY")
    print("=" * 70)

    print("\nVariant Categories:")
    for cat in CATEGORY_ORDER:
        color = CATEGORY_COLORS[cat]
        print(f"  {cat:<15} : {color}")

    print("\nPipeline Stages:")
    for i, stage in enumerate(VCF_STAGE_ORDER, 1):
        display = STAGE_DISPLAY_NAMES.get(stage, stage.title())
        print(f"  {i}. {stage:<20} -> {display}")

    print("\nTier System:")
    for tier in TIER_ORDER:
        desc = TIER_DEFINITIONS[tier]["description"]
        print(f"  {tier}: {desc}")

    print("\n" + "=" * 70)
