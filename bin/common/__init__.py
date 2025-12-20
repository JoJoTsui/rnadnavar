"""
Common shared configuration and utilities for VCF processing pipeline.

This package contains shared configuration used by both vcf_utils (pipeline)
and vcf_stats (analysis notebooks).
"""

from .vcf_config import (
    # Categories
    CATEGORY_ORDER,
    CATEGORY_COLORS,
    
    # Stages
    VCF_STAGE_ORDER,
    STAGE_DISPLAY_NAMES,
    
    # Tiering
    TIER_DEFINITIONS,
    TIER_ORDER,
    
    # Callers
    TOOLS,
    MODALITIES,
    
    # Utility functions
    get_category_color,
    get_stage_display_name,
    is_annotation_stage,
    validate_category,
    validate_stage,
)

__all__ = [
    "CATEGORY_ORDER",
    "CATEGORY_COLORS",
    "VCF_STAGE_ORDER",
    "STAGE_DISPLAY_NAMES",
    "TIER_DEFINITIONS",
    "TIER_ORDER",
    "TOOLS",
    "MODALITIES",
    "get_category_color",
    "get_stage_display_name",
    "is_annotation_stage",
    "validate_category",
    "validate_stage",
]
