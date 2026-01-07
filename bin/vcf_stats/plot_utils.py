#!/usr/bin/env python3
"""
Shared Plot Utilities

Unified helpers for legend handling, color application, stage ordering,
and heatmap matrix normalization to ensure consistent visuals across plots.
"""

from typing import Dict, List, Tuple, Set, Optional
import numpy as np
import pandas as pd


def build_legend_tracker() -> Set[str]:
    """Return a set to track which categories have been added to the legend."""
    return set()


def should_add_to_legend(categories_seen: Set[str], category: str) -> bool:
    """Return True if this category should be added to legend now (first occurrence)."""
    add = category not in categories_seen
    if add:
        categories_seen.add(category)
    return add


def apply_category_colors(categories: List[str], color_map: Dict[str, str]) -> List[str]:
    """Map categories to colors using provided color_map with gray fallback."""
    return [color_map.get(cat, "#8A8A8A") for cat in categories]


def legend_config(position: str = "right") -> Dict:
    """Return a consistent legend configuration dict for plotly figures."""
    if position == "top":
        return dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
    # default right side vertical
    return dict(orientation="v", yanchor="top", y=1.0, xanchor="left", x=1.02)


def sort_stages_by_order(stages: List[str], order_list: List[str]) -> List[str]:
    """Sort stage names according to their order list, preserving those not found at the end."""
    order_index = {name: idx for idx, name in enumerate(order_list)}
    return sorted(stages, key=lambda s: order_index.get(s, 999))


def heatmap_matrix(df: pd.DataFrame, normalize: Optional[str] = None) -> pd.DataFrame:
    """
    Build a heatmap matrix pivoted by Category x Stage with optional normalization.
    normalize:
      - None: raw counts
      - "stage": column-wise percent normalization (per stage)
    """
    mat = df.pivot_table(index="Category", columns="Stage", values="Count", fill_value=0)
    if normalize == "stage":
        # Avoid divide-by-zero with replace(0,1)
        col_sums = mat.sum(axis=0).replace(0, 1)
        mat = (mat / col_sums) * 100.0
    return mat


def percentile_cap(values: np.ndarray, pct: float = 95.0) -> float:
    """Return the percentile cap for zmax scaling in heatmaps."""
    try:
        return float(np.percentile(values.flatten(), pct))
    except Exception:
        return float(np.max(values) if values.size else 0.0)
