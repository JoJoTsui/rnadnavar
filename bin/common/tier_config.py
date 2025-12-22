#!/usr/bin/env python3
"""
Tier Configuration Module for Hybrid Tiering System

IMPORTANT: This module IMPORTS thresholds from annotation_config.py (centralized source of truth).
Do NOT modify threshold values here - update annotation_config.py instead.

This module defines the two-dimensional tiering system that combines:
1. Caller Support (C1-C7): Based on category-concordant DNA/RNA caller counts
2. Database Evidence (D0-D1): Based on presence in external databases

Final tiers use CxDy notation for clarity and maintainability.

Usage:
    from tier_config import CALLER_TIER_RULES, DATABASE_TIER_RULES, TIER_ORDER
    
    # Get tier display name
    tier_name = TIER_DISPLAY_NAMES["C1D1"]  # "C1D1: ≥2 DNA + ≥2 RNA, Database"
"""

from typing import Dict, List, Tuple, Callable

# Import centralized thresholds from annotation_config
try:
    from annotation_config import GNOMAD_THRESHOLDS, COSMIC_THRESHOLDS, REDIPORTAL_THRESHOLDS
except ImportError:
    # Fallback if annotation_config not available (for backward compatibility)
    GNOMAD_THRESHOLDS = {"germline_frequency": 0.001}
    COSMIC_THRESHOLDS = {"recurrence_minimum": 5}
    REDIPORTAL_THRESHOLDS = {"min_rna_support": 2}

# ============================================================================
# CALLER SUPPORT TIER DEFINITIONS (C1-C7)
# ============================================================================
#
# Caller tiers assess variant quality based on number of DNA and RNA callers
# that vote for the SAME category as the final FILTER classification.
# Only category-concordant callers are counted.
#
# Tier assignment is mutually exclusive - first match wins in order C1→C7
#

CALLER_TIER_RULES = {
    "C1": {
        "description": "≥2 DNA + ≥2 RNA (both modalities strong)",
        "condition": lambda dna, rna: dna >= 2 and rna >= 2,
        "priority": 1,
    },
    "C2": {
        "description": "≥2 DNA + (0 or 1) RNA (DNA-strong, RNA weak/absent)",
        "condition": lambda dna, rna: dna >= 2 and rna <= 1,
        "priority": 2,
    },
    "C3": {
        "description": "≥2 RNA + (0 or 1) DNA (RNA-strong, DNA weak/absent)",
        "condition": lambda dna, rna: rna >= 2 and dna <= 1,
        "priority": 3,
    },
    "C4": {
        "description": "1 DNA + 1 RNA (both modalities weak)",
        "condition": lambda dna, rna: dna == 1 and rna == 1,
        "priority": 4,
    },
    "C5": {
        "description": "1 DNA + 0 RNA (DNA-only weak)",
        "condition": lambda dna, rna: dna == 1 and rna == 0,
        "priority": 5,
    },
    "C6": {
        "description": "0 DNA + 1 RNA (RNA-only weak)",
        "condition": lambda dna, rna: dna == 0 and rna == 1,
        "priority": 6,
    },
    "C7": {
        "description": "0 DNA + 0 RNA (no caller support)",
        "condition": lambda dna, rna: dna == 0 and rna == 0,
        "priority": 7,
    },
}

# Ordered list of caller tiers by priority (for assignment logic)
CALLER_TIER_ORDER = ["C1", "C2", "C3", "C4", "C5", "C6", "C7"]


# ============================================================================
# DATABASE EVIDENCE TIER DEFINITIONS (D0-D1)
# ============================================================================
#
# Database tiers assess whether a variant is present in external databases:
# - gnomAD: Population frequency database (germline variants)
# - COSMIC: Cancer somatic mutation database
# - REDIportal/DARNED: RNA editing databases
#
# D1 = Present in at least one database with significance
# D0 = Not present in any database (novel variant)
#

DATABASE_TIER_RULES = {
    "D1": {
        "description": "Has database support",
        "databases": ["gnomAD", "COSMIC", "REDIportal", "DARNED"],
        "thresholds": {
            "gnomAD_AF": GNOMAD_THRESHOLDS.get("germline_frequency", 0.001),  # From annotation_config
            "COSMIC_CNT": COSMIC_THRESHOLDS.get("recurrence_minimum", 0),      # From annotation_config
            "REDIportal": True,       # Present in REDIportal
            "DARNED": True,           # Present in DARNED
        },
    },
    "D0": {
        "description": "No database support",
    },
}

# Ordered list of database tiers (D1 first for prioritization)
DATABASE_TIER_ORDER = ["D1", "D0"]


# ============================================================================
# COMBINED TIER MATRIX (CxDy)
# ============================================================================
#
# Final tier is the combination of Caller tier × Database tier
# Total: 7 caller tiers × 2 database tiers = 14 final tiers
#
# Format: CxDy where x ∈ {1,2,3,4,5,6,7} and y ∈ {0,1}
#
# Quality ranking (best to worst):
#   C1D1 > C1D0 > C2D1 > C2D0 > C3D1 > C3D0 > ... > C7D1 > C7D0
#

# Generate all 14 possible tier combinations
TIER_COMBINATIONS: List[Tuple[str, str]] = [
    (caller_tier, db_tier)
    for caller_tier in CALLER_TIER_ORDER
    for db_tier in DATABASE_TIER_ORDER
]

# Ordered list of all final tiers (for visualization and sorting)
TIER_ORDER: List[str] = [
    f"{caller}{db}" for caller, db in TIER_COMBINATIONS
]

# Display names for all final tiers
TIER_DISPLAY_NAMES: Dict[str, str] = {}
for caller_tier in CALLER_TIER_ORDER:
    caller_desc = CALLER_TIER_RULES[caller_tier]["description"]
    for db_tier in DATABASE_TIER_ORDER:
        db_desc = DATABASE_TIER_RULES[db_tier]["description"]
        tier_name = f"{caller_tier}{db_tier}"
        TIER_DISPLAY_NAMES[tier_name] = f"{tier_name}: {caller_desc}, {db_desc}"

# Short descriptions for compact display
TIER_SHORT_NAMES: Dict[str, str] = {
    "C1D1": "≥2D+≥2R,DB", "C1D0": "≥2D+≥2R",
    "C2D1": "≥2D,DB", "C2D0": "≥2D",
    "C3D1": "≥2R,DB", "C3D0": "≥2R",
    "C4D1": "1D+1R,DB", "C4D0": "1D+1R",
    "C5D1": "1D,DB", "C5D0": "1D",
    "C6D1": "1R,DB", "C6D0": "1R",
    "C7D1": "0D+0R,DB", "C7D0": "0D+0R",
}


# ============================================================================
# TIER QUALITY SCORES
# ============================================================================
#
# Numeric scores for ranking tiers by quality (higher = better)
# Used for filtering, sorting, and quality assessment
#

TIER_QUALITY_SCORES: Dict[str, int] = {
    # C1 tiers: Highest quality (both modalities strong)
    "C1D1": 140, "C1D0": 130,
    # C2 tiers: DNA-strong
    "C2D1": 120, "C2D0": 110,
    # C3 tiers: RNA-strong
    "C3D1": 100, "C3D0": 90,
    # C4 tiers: Both modalities weak but present
    "C4D1": 80, "C4D0": 70,
    # C5 tiers: DNA-only weak
    "C5D1": 60, "C5D0": 50,
    # C6 tiers: RNA-only weak
    "C6D1": 40, "C6D0": 30,
    # C7 tiers: No caller support (lowest quality)
    "C7D1": 20, "C7D0": 10,
}


# ============================================================================
# CATEGORY-SPECIFIC TIER EXPECTATIONS
# ============================================================================
#
# Some categories typically fall into specific tiers:
# - RNA_Edit: Usually C7D1 (no callers, but in REDIportal/DARNED)
# - NoConsensus: Usually C7D0 or C7D1 (DNA/RNA disagree, no consensus)
# - Somatic: Expected in C1-C6 range (caller-supported)
# - Germline: Often C1D1, C2D1 (high caller + gnomAD support)
#

CATEGORY_TYPICAL_TIERS = {
    "Somatic": ["C1D1", "C1D0", "C2D1", "C2D0", "C3D1", "C3D0"],
    "Germline": ["C1D1", "C2D1", "C3D1", "C4D1"],
    "Reference": ["C1D1", "C2D1"],  # Usually strong caller + gnomAD support
    "Artifact": ["C5D0", "C6D0", "C7D0"],  # Weak caller, no database
    "RNA_Edit": ["C7D1"],  # No callers, REDIportal/DARNED database
    "NoConsensus": ["C7D0", "C7D1"],  # No consensus, variable database
}


# ============================================================================
# TIER COLOR SCHEME
# ============================================================================
#
# Color palette for tier visualization (follows quality gradient)
# Higher quality tiers = cooler colors (blue/green)
# Lower quality tiers = warmer colors (orange/red)
#

TIER_COLORS: Dict[str, str] = {
    # C1 tiers: Dark blue (highest quality)
    "C1D1": "#1f77b4", "C1D0": "#4292c6",
    # C2 tiers: Blue
    "C2D1": "#6baed6", "C2D0": "#9ecae1",
    # C3 tiers: Cyan/Teal
    "C3D1": "#41ab5d", "C3D0": "#74c476",
    # C4 tiers: Green
    "C4D1": "#a1d99b", "C4D0": "#c7e9c0",
    # C5 tiers: Yellow
    "C5D1": "#fdae6b", "C5D0": "#fdd0a2",
    # C6 tiers: Orange
    "C6D1": "#fd8d3c", "C6D0": "#fdae6b",
    # C7 tiers: Red/Gray (lowest quality)
    "C7D1": "#d94801", "C7D0": "#a63603",
}


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def get_tier_quality(tier: str) -> int:
    """
    Get numeric quality score for a tier.
    
    Args:
        tier: Tier name (e.g., "C1D1", "C7D0")
    
    Returns:
        Quality score (higher = better), or 0 if tier not recognized
    """
    return TIER_QUALITY_SCORES.get(tier, 0)


def get_tier_color(tier: str) -> str:
    """
    Get color code for a tier.
    
    Args:
        tier: Tier name (e.g., "C1D1")
    
    Returns:
        Hex color code, or gray if tier not recognized
    """
    return TIER_COLORS.get(tier, "#999999")


def get_tier_display_name(tier: str, short: bool = False) -> str:
    """
    Get human-readable display name for a tier.
    
    Args:
        tier: Tier name (e.g., "C1D1")
        short: If True, return short name (e.g., "≥2D+≥2R,DB")
    
    Returns:
        Display name
    """
    if short:
        return TIER_SHORT_NAMES.get(tier, tier)
    return TIER_DISPLAY_NAMES.get(tier, tier)


def parse_tier(tier: str) -> Tuple[str, str]:
    """
    Parse a tier name into caller and database components.
    
    Args:
        tier: Tier name (e.g., "C1D1")
    
    Returns:
        Tuple of (caller_tier, database_tier), e.g., ("C1", "D1")
    
    Raises:
        ValueError: If tier format is invalid
    """
    if not tier or len(tier) < 4:
        raise ValueError(f"Invalid tier format: {tier}")
    
    # Expected format: CxDy where x is 1-7, y is 0-1
    if tier[0] != 'C' or tier[2] != 'D':
        raise ValueError(f"Invalid tier format: {tier} (expected CxDy)")
    
    caller_tier = tier[:2]  # e.g., "C1"
    db_tier = tier[2:]      # e.g., "D1"
    
    if caller_tier not in CALLER_TIER_ORDER:
        raise ValueError(f"Invalid caller tier: {caller_tier}")
    if db_tier not in DATABASE_TIER_ORDER:
        raise ValueError(f"Invalid database tier: {db_tier}")
    
    return caller_tier, db_tier


def is_high_quality_tier(tier: str, threshold: int = 100) -> bool:
    """
    Check if a tier meets high quality threshold.
    
    Args:
        tier: Tier name (e.g., "C1D1")
        threshold: Minimum quality score (default: 100 = C3D1 and above)
    
    Returns:
        True if tier quality >= threshold
    """
    return get_tier_quality(tier) >= threshold


def filter_tiers_by_quality(tiers: List[str], min_quality: int) -> List[str]:
    """
    Filter tiers by minimum quality score.
    
    Args:
        tiers: List of tier names
        min_quality: Minimum quality score threshold
    
    Returns:
        Filtered list of tiers meeting quality threshold
    """
    return [tier for tier in tiers if get_tier_quality(tier) >= min_quality]


# ============================================================================
# BACKWARD COMPATIBILITY MAPPING (T1-T8 → CxDy)
# ============================================================================
#
# Mapping from old tier system (T1-T8) to new CxDy system
# Used for migrating existing code and data
#

OLD_TO_NEW_TIER_MAPPING = {
    "T1": "C1D0",  # 2+ DNA + 2+ RNA → C1 (assume no DB for old data)
    "T2": "C2D0",  # 2+ DNA + 1 RNA → C2 (DNA strong, RNA weak)
    "T3": "C2D0",  # 2+ DNA only → C2
    "T4": "C4D0",  # 1 DNA + 1+ RNA → C4 (both weak)
    "T5": "C5D0",  # 1 DNA only → C5
    "T6": "C3D0",  # 2+ RNA only → C3
    "T7": "C6D0",  # 1 RNA only → C6
    "T8": "C7D0",  # Special cases → C7
}


def convert_old_tier_to_new(old_tier: str) -> str:
    """
    Convert old tier system (T1-T8) to new CxDy system.
    
    Args:
        old_tier: Old tier name (e.g., "T1")
    
    Returns:
        New tier name (e.g., "C1D0")
    """
    return OLD_TO_NEW_TIER_MAPPING.get(old_tier, "C7D0")


if __name__ == "__main__":
    # Print configuration summary
    print("=" * 80)
    print("TIER CONFIGURATION SUMMARY")
    print("=" * 80)
    
    print("\n📊 CALLER SUPPORT TIERS (C1-C7):")
    for tier in CALLER_TIER_ORDER:
        desc = CALLER_TIER_RULES[tier]["description"]
        print(f"  {tier}: {desc}")
    
    print("\n🗄️  DATABASE EVIDENCE TIERS (D0-D1):")
    for tier in DATABASE_TIER_ORDER:
        desc = DATABASE_TIER_RULES[tier]["description"]
        print(f"  {tier}: {desc}")
    
    print(f"\n🎯 FINAL TIER COMBINATIONS: {len(TIER_ORDER)} total")
    print("  Quality ranking (best to worst):")
    for i, tier in enumerate(TIER_ORDER, 1):
        score = get_tier_quality(tier)
        color = get_tier_color(tier)
        short = get_tier_display_name(tier, short=True)
        print(f"    {i:2d}. {tier} ({score:3d}): {short:<15} {color}")
    
    print("\n" + "=" * 80)
