#!/usr/bin/env python3
"""
Unified Classification and Voting Configuration Module

This module centralizes all classification categories, voting thresholds, and
configuration parameters used throughout the RNADNAvar pipeline. It ensures
consistency across all classification stages: standalone caller classification,
consensus voting, cross-modality rescue, and annotation-based reclassification.

Key Features:
- Centralized classification categories including RNA editing
- Default thresholds for all classification stages
- Canonical chromosome definitions
- Quality filtering parameters (disabled by default)
- Extensible configuration system
"""

# Classification categories (including RNA editing)
CATEGORIES = ["Somatic", "Germline", "Reference", "Artifact", "NoConsensus", "RNAedit"]

# Valid biological categories for FILTER field
BIOLOGICAL_CATEGORIES = ["Somatic", "Germline", "Reference", "Artifact"]

# Default thresholds for all classification stages
DEFAULT_THRESHOLDS = {
    # Consensus voting thresholds
    "consensus_snv_threshold": 2,  # Minimum callers for SNV consensus
    "consensus_indel_threshold": 2,  # Minimum callers for indel consensus
    # Annotation-based reclassification thresholds
    # Note: Maps to Nextflow param cosmic_gnomad_germline_freq_threshold
    "annotation_germline_freq_threshold": 0.001,  # gnomAD AF for germline
    # Note: Maps to Nextflow param cosmic_gnomad_somatic_consensus_threshold
    "annotation_somatic_consensus_threshold": 2,  # Min caller support for somatic
    # Note: Maps to Nextflow param cosmic_gnomad_cosmic_recurrence_threshold
    "annotation_cosmic_recurrence_threshold": 5,  # COSMIC count for somatic rescue
    # Cross-modality thresholds for rescue
    "cross_modality_min_support": 1,  # Min callers per modality
    "cross_modality_min_callers_for_artifact": 2,  # Min callers to declare artifact on disagreement
    # RNA editing classification thresholds
    # Note: Maps to Nextflow param min_rna_support
    "rna_editing_min_rna_support": 2,  # Min RNA callers for RNAedit
    "rna_editing_min_support": 2,  # Alias for compatibility
    # Filtering thresholds (used by filter_vcf.py)
    "filtering_gnomad_threshold": 0.001,  # gnomAD AF threshold for filtering
    "filtering_min_alt_reads": 2,  # Minimum alt reads for filtering
    # Quality thresholds (disabled by default)
    "quality_filtering_enabled": False,
    "min_depth": 10,
    "min_vaf": 0.05,
    # Strelka-specific thresholds
    "strelka_min_normal_depth": 2,  # Min normal depth for germline/reference rescue
}

# Chromosome configuration
CANONICAL_CHROMOSOMES = {
    "autosomes": [str(i) for i in range(1, 23)],  # 1-22
    "sex": ["X", "Y"],
    "mitochondrial": ["M", "MT"],  # Both M and MT are accepted
}

# Build canonical chromosome set for quick lookup
CANONICAL_SET = set()
CANONICAL_SET.update(CANONICAL_CHROMOSOMES["autosomes"])
CANONICAL_SET.update(CANONICAL_CHROMOSOMES["sex"])
CANONICAL_SET.update(CANONICAL_CHROMOSOMES["mitochondrial"])

# Canonical chromosome filtering default
FILTER_NON_CANONICAL_DEFAULT = (
    True  # Filter by default, can disable with --include-non-canonical
)

# Caller name patterns for identification
CALLER_PATTERNS = {
    "strelka": ["strelka"],
    "deepsomatic": ["deepsomatic", "deep_somatic"],
    "mutect2": ["mutect2", "mutect"],
    "consensus": ["consensus"],
}

# Priority order for classification when there are conflicts
CLASSIFICATION_PRIORITY = [
    "Somatic",
    "Germline",
    "Reference",
    "RNAedit",
    "Artifact",
    "NoConsensus",
]

# INFO field names for classification metadata
CLASSIFICATION_INFO_FIELDS = {
    "classification": "VC",  # Variant Classification
    "classification_callers": "VC_CALLERS",  # Callers contributing to classification
    "classification_confidence": "VC_CONF",  # Classification confidence score
    "rna_editing_flag": "RNA_EDIT",  # RNA editing annotation flag
    "cross_modality_flag": "CROSS_MODALITY",  # Cross-modality support flag
    "rescued_flag": "RESCUED",  # Rescue classification flag
}

# Filter flag descriptions for VCF header
FILTER_DESCRIPTIONS = {
    "Somatic": "High-confidence somatic variant specific to tumor",
    "Germline": "Germline variant detected in normal sample",
    "Reference": "Reference call - no variant detected",
    "Artifact": "Low quality variant, technical artifact, or inconsistent classifications across callers",
    "NoConsensus": "Does not meet consensus threshold or does not fit other classification categories",
    "RNAedit": "RNA editing event confirmed by REDIportal database and RNA evidence",
}


def get_config_with_overrides(overrides=None):
    """
    Get configuration with optional overrides.

    Args:
        overrides (dict): Dictionary of threshold overrides

    Returns:
        dict: Configuration with overrides applied
    """
    config = DEFAULT_THRESHOLDS.copy()
    if overrides:
        config.update(overrides)
    return config


def validate_classification(classification):
    """
    Validate that a classification is one of the allowed categories.

    Args:
        classification (str): Classification to validate

    Returns:
        bool: True if valid classification
    """
    return classification in CATEGORIES


def get_biological_classification(classification):
    """
    Convert any classification to a biological category suitable for FILTER field.

    Args:
        classification (str): Input classification

    Returns:
        str: One of the biological categories (Somatic/Germline/Reference/Artifact)
    """
    if classification in BIOLOGICAL_CATEGORIES:
        return classification
    elif classification == "NoConsensus":
        return "Artifact"  # NoConsensus interpreted as Artifact for FILTER field
    elif classification == "RNAedit":
        return "Somatic"  # RNA edits are somatic events
    else:
        return "Artifact"  # Unknown classifications become Artifact


def format_classification_for_filter(classification):
    """
    Format classification for VCF FILTER field according to VCF specifications.

    Args:
        classification (str): Classification category

    Returns:
        str: Formatted string for FILTER field
    """
    # For VCF FILTER field, we only use biological categories
    bio_class = get_biological_classification(classification)
    return bio_class


def parse_caller_name(caller_path_or_name):
    """
    Parse caller name from file path or name string.

    Args:
        caller_path_or_name (str): File path or caller name

    Returns:
        str: Standardized caller name
    """
    name_lower = str(caller_path_or_name).lower()

    for caller, patterns in CALLER_PATTERNS.items():
        for pattern in patterns:
            if pattern in name_lower:
                return caller

    # Default: return cleaned name
    import os

    if "/" in caller_path_or_name:
        return os.path.basename(caller_path_or_name).split(".")[0]
    return caller_path_or_name
