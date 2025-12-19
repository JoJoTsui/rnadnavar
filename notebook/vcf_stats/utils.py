#!/usr/bin/env python3
"""
Utility Functions for VCF Statistics Module

Reusable functions for stage ordering, legend management, and plot creation
to eliminate code duplication across visualization and analysis modules.
"""

from typing import List, Dict, Any, Optional, Tuple
import pandas as pd

# Import constants (will be available after __init__.py loads this module)
# To avoid circular imports, we'll import these inside functions that need them


def get_stage_order_key(stage: str) -> int:
    """
    Get numeric order for a VCF processing stage.
    
    Args:
        stage: Stage name (e.g., 'rescue', 'dna_consensus', 'cosmic_gnomad')
        
    Returns:
        int: Order index (lower = earlier in pipeline), 999 for unknown stages
    """
    from . import VCF_STAGE_ORDER
    try:
        return VCF_STAGE_ORDER.index(stage)
    except ValueError:
        return 999  # Unknown stages go last


def sort_stages(stages: List[str]) -> List[str]:
    """
    Sort stages by processing order using VCF_STAGE_ORDER.
    
    Args:
        stages: List of stage names to sort
        
    Returns:
        List of stages sorted by pipeline order
        
    Example:
        >>> sort_stages(['filtered_rescue', 'rescue', 'cosmic_gnomad'])
        ['rescue', 'cosmic_gnomad', 'filtered_rescue']
    """
    return sorted(stages, key=get_stage_order_key)


def should_show_legend(subplot_index: int, first_subplot: int = 1) -> bool:
    """
    Determine if legend should be shown for this subplot to avoid duplication.
    
    Args:
        subplot_index: Current subplot index (1-based)
        first_subplot: Which subplot should show the legend (default: 1)
        
    Returns:
        bool: True if this subplot should display the legend
    """
    return subplot_index == first_subplot


def create_legend_config(orientation: str = "v", position: str = "right") -> Dict[str, Any]:
    """
    Create standardized legend configuration for plotly plots.
    
    Args:
        orientation: 'v' for vertical or 'h' for horizontal
        position: 'right', 'top', 'bottom' for legend placement
        
    Returns:
        dict: Plotly legend configuration
    """
    configs = {
        "right": {
            "orientation": "v",
            "yanchor": "top",
            "y": 1.0,
            "xanchor": "left",
            "x": 1.02,
            "title": "Category"
        },
        "top": {
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "right",
            "x": 1,
            "title": None
        },
        "bottom": {
            "orientation": "h",
            "yanchor": "top",
            "y": -0.15,
            "xanchor": "center",
            "x": 0.5,
            "title": None
        }
    }
    
    config = configs.get(position, configs["right"]).copy()
    if orientation and position == "right":
        config["orientation"] = orientation
    
    return config


def filter_no_consensus(df: pd.DataFrame, category_col: str = "Category") -> pd.DataFrame:
    """
    Filter out NoConsensus variants from a DataFrame.
    
    Args:
        df: DataFrame containing variant data
        category_col: Name of the column containing category labels
        
    Returns:
        DataFrame with NoConsensus rows removed
    """
    if category_col not in df.columns:
        return df
    return df[df[category_col] != "NoConsensus"].copy()


def get_stage_display_name(stage: str) -> str:
    """
    Get human-readable display name for a stage.
    
    Args:
        stage: Stage key (e.g., 'dna_consensus', 'cosmic_gnomad')
        
    Returns:
        str: Display name (e.g., 'DNA Consensus', 'COSMIC/GnomAD')
    """
    from . import STAGE_DISPLAY_NAMES
    return STAGE_DISPLAY_NAMES.get(stage, stage.replace("_", " ").title())


def create_dual_view_plots(
    df: pd.DataFrame,
    plot_func: callable,
    category_col: str = "Category",
    exclude_no_consensus: bool = True,
    title_suffix_full: str = "",
    title_suffix_filtered: str = " (excluding NoConsensus)"
) -> Tuple[Any, Optional[Any]]:
    """
    Create dual-view plots: full and NoConsensus-free.
    
    Args:
        df: DataFrame containing plot data
        plot_func: Function that takes df and title_suffix, returns plotly figure
        category_col: Column name containing categories
        exclude_no_consensus: Whether to create second plot
        title_suffix_full: Suffix for full plot title
        title_suffix_filtered: Suffix for filtered plot title
        
    Returns:
        Tuple of (full_fig, filtered_fig) where filtered_fig is None if not requested
    """
    # Create full plot
    full_fig = plot_func(df, title_suffix_full)
    
    # Create filtered plot if requested
    filtered_fig = None
    if exclude_no_consensus:
        df_filtered = filter_no_consensus(df, category_col)
        if not df_filtered.empty:
            filtered_fig = plot_func(df_filtered, title_suffix_filtered)
    
    return full_fig, filtered_fig


def ensure_all_categories_in_legend(
    categories_to_show: List[str],
    categories_in_data: set,
    category_order: Optional[List[str]] = None
) -> List[str]:
    """
    Ensure all expected categories appear in legend even if count is zero.
    
    Args:
        categories_to_show: Categories that should appear in legend
        categories_in_data: Categories actually present in the data
        category_order: Order to use (imports CATEGORY_ORDER if None)
        
    Returns:
        List of categories to include in legend, in correct order
    """
    if category_order is None:
        from . import CATEGORY_ORDER
        category_order = CATEGORY_ORDER
    
    # Include all requested categories that are in the standard order
    result = []
    for cat in category_order:
        if cat in categories_to_show or cat in categories_in_data:
            result.append(cat)
    
    # Add any extra categories not in standard order
    for cat in categories_to_show:
        if cat not in result:
            result.append(cat)
    
    return result


def collect_stage_statistics(
    all_vcf_stats: Dict[str, Any],
    stages: List[str]
) -> Dict[str, Dict[str, int]]:
    """
    Collect classification statistics for specified stages.
    
    Args:
        all_vcf_stats: Complete VCF statistics dictionary
        stages: List of stage names to collect
        
    Returns:
        Dict mapping stage names to classification dictionaries
    """
    stage_stats = {}
    
    for stage in stages:
        classification = {}
        total = 0
        
        if stage in all_vcf_stats:
            for name, data in all_vcf_stats[stage].items():
                if "stats" in data and "basic" in data["stats"]:
                    basic = data["stats"]["basic"]
                    cls = basic.get("classification", {})
                    tot = basic.get("total_variants", 0)
                    
                    # Merge counts
                    total += tot
                    for cat, count in cls.items():
                        classification[cat] = classification.get(cat, 0) + count
        
        stage_stats[stage] = {
            "classification": classification,
            "total_variants": total
        }
    
    return stage_stats


print("✓ VCF statistics utility functions loaded")
