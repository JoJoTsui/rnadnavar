#!/usr/bin/env python3
"""
Indel size distribution analysis utilities.

Provides functions for categorizing indels into fixed bins, computing per-sample
and aggregated statistics, and returning tabular summaries ready for plotting
or export. Designed to be pipeline-agnostic and reusable across notebooks or
workflow steps.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

INDEL_SIZE_BINS: List[str] = ["1", "2", "3", "4-14", "15-29", "30-49", ">50"]

INSERTIONS_KEY = "insertions"
DELETIONS_KEY = "deletions"


def categorize_indel_size(size: int) -> str:
    """Map an indel length to a discrete size bin.

    Args:
        size: Positive for insertions, negative for deletions.

    Returns:
        Bin label drawn from INDEL_SIZE_BINS.
    """
    abs_size = abs(size)

    if abs_size == 0:
        # Treat zero-length as the smallest bin; should be rare in real data.
        return "1"
    if abs_size == 1:
        return "1"
    if abs_size == 2:
        return "2"
    if abs_size == 3:
        return "3"
    if 4 <= abs_size <= 14:
        return "4-14"
    if 15 <= abs_size <= 29:
        return "15-29"
    if 30 <= abs_size <= 49:
        return "30-49"
    return ">50"


def indel_stats_by_bin(indel_sizes: List[int]) -> Dict[str, Dict[str, int]]:
    """Count insertions and deletions by predefined bins."""
    bin_stats = {
        INSERTIONS_KEY: {bin_label: 0 for bin_label in INDEL_SIZE_BINS},
        DELETIONS_KEY: {bin_label: 0 for bin_label in INDEL_SIZE_BINS},
    }

    for size in indel_sizes:
        bin_label = categorize_indel_size(size)
        if size > 0:
            bin_stats[INSERTIONS_KEY][bin_label] += 1
        else:
            bin_stats[DELETIONS_KEY][bin_label] += 1

    return bin_stats


def summarize_indel_sizes(indel_sizes: List[int]) -> Dict[str, float]:
    """Compute summary statistics for a list of indel sizes."""
    if not indel_sizes:
        return {
            "total": 0,
            "insertions": 0,
            "deletions": 0,
            "median_abs": 0.0,
            "max_insertion": 0,
            "max_deletion_abs": 0,
            "mean_insertion": 0.0,
            "mean_deletion_abs": 0.0,
        }

    insertions = [s for s in indel_sizes if s > 0]
    deletions = [s for s in indel_sizes if s < 0]

    return {
        "total": len(indel_sizes),
        "insertions": len(insertions),
        "deletions": len(deletions),
        "median_abs": float(np.median(np.abs(indel_sizes))),
        "max_insertion": max(insertions) if insertions else 0,
        "max_deletion_abs": abs(min(deletions)) if deletions else 0,
        "mean_insertion": float(np.mean(insertions)) if insertions else 0.0,
        "mean_deletion_abs": float(np.mean(np.abs(deletions))) if deletions else 0.0,
    }


def aggregate_indel_statistics(
    all_vcf_stats: Dict[str, Dict],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Aggregate indel statistics across all samples/models.

    Returns:
        summary_df: Per-file summary metrics (counts and size statistics).
        aggregated_binned_df: Global binned counts and percentages across samples.
        per_sample_binned_df: Per-sample binned counts and percentages.
    """
    summary_rows: List[Dict] = []
    aggregated_insertions = {bin_label: 0 for bin_label in INDEL_SIZE_BINS}
    aggregated_deletions = {bin_label: 0 for bin_label in INDEL_SIZE_BINS}
    per_sample_rows: List[Dict] = []

    for _, data in all_vcf_stats.items():
        metadata = data.get("metadata", {})
        sizes: List[int] = data.get("stats", {}).get("indel_sizes", [])
        if not sizes:
            continue

        insertions = [s for s in sizes if s > 0]
        deletions = [s for s in sizes if s < 0]

        bin_stats = indel_stats_by_bin(sizes)

        total_insertions = sum(bin_stats[INSERTIONS_KEY].values())
        total_deletions = sum(bin_stats[DELETIONS_KEY].values())

        for bin_label in INDEL_SIZE_BINS:
            per_sample_rows.append(
                {
                    "Sample": metadata.get("sample", ""),
                    "Model": metadata.get("model", ""),
                    "Size_Bin": bin_label,
                    "Insertions_Count": bin_stats[INSERTIONS_KEY][bin_label],
                    "Deletions_Count": bin_stats[DELETIONS_KEY][bin_label],
                    "Total_Insertions": total_insertions,
                    "Total_Deletions": total_deletions,
                }
            )
            aggregated_insertions[bin_label] += bin_stats[INSERTIONS_KEY][bin_label]
            aggregated_deletions[bin_label] += bin_stats[DELETIONS_KEY][bin_label]

        size_summary = summarize_indel_sizes(sizes)
        summary_rows.append(
            {
                "Sample": metadata.get("sample", ""),
                "Model": metadata.get("model", ""),
                "Total_INDELs": size_summary["total"],
                "Insertions_Count": size_summary["insertions"],
                "Insertions_Pct": (
                    size_summary["insertions"] / size_summary["total"] * 100
                )
                if size_summary["total"]
                else 0,
                "Deletions_Count": size_summary["deletions"],
                "Deletions_Pct": (
                    size_summary["deletions"] / size_summary["total"] * 100
                )
                if size_summary["total"]
                else 0,
                "Median_Size": size_summary["median_abs"],
                "Max_Insertion": size_summary["max_insertion"],
                "Max_Deletion": size_summary["max_deletion_abs"],
                "Mean_Insertion": size_summary["mean_insertion"],
                "Mean_Deletion": size_summary["mean_deletion_abs"],
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    if not summary_df.empty:
        for col in [
            "Insertions_Pct",
            "Deletions_Pct",
            "Median_Size",
            "Mean_Insertion",
            "Mean_Deletion",
        ]:
            if col in summary_df.columns:
                summary_df[col] = summary_df[col].round(2)
        summary_df = summary_df.sort_values(["Sample", "Model"]).reset_index(drop=True)

    total_insertions_all = sum(aggregated_insertions.values())
    total_deletions_all = sum(aggregated_deletions.values())

    aggregated_rows = []
    for bin_label in INDEL_SIZE_BINS:
        aggregated_rows.append(
            {
                "Size_Bin": bin_label,
                "Insertions_Count": aggregated_insertions[bin_label],
                "Insertions_Pct": (
                    aggregated_insertions[bin_label] / total_insertions_all * 100
                )
                if total_insertions_all
                else 0,
                "Deletions_Count": aggregated_deletions[bin_label],
                "Deletions_Pct": (
                    aggregated_deletions[bin_label] / total_deletions_all * 100
                )
                if total_deletions_all
                else 0,
            }
        )

    aggregated_binned_df = pd.DataFrame(aggregated_rows)
    if not aggregated_binned_df.empty:
        aggregated_binned_df["Insertions_Pct"] = aggregated_binned_df[
            "Insertions_Pct"
        ].round(2)
        aggregated_binned_df["Deletions_Pct"] = aggregated_binned_df[
            "Deletions_Pct"
        ].round(2)

    per_sample_binned_df = pd.DataFrame(per_sample_rows)
    if not per_sample_binned_df.empty:
        per_sample_binned_df["Insertions_Pct"] = per_sample_binned_df.apply(
            lambda row: (row["Insertions_Count"] / row["Total_Insertions"] * 100)
            if row["Total_Insertions"]
            else 0,
            axis=1,
        )
        per_sample_binned_df["Deletions_Pct"] = per_sample_binned_df.apply(
            lambda row: (row["Deletions_Count"] / row["Total_Deletions"] * 100)
            if row["Total_Deletions"]
            else 0,
            axis=1,
        )
        per_sample_binned_df["Insertions_Pct"] = per_sample_binned_df[
            "Insertions_Pct"
        ].round(2)
        per_sample_binned_df["Deletions_Pct"] = per_sample_binned_df[
            "Deletions_Pct"
        ].round(2)

    return summary_df, aggregated_binned_df, per_sample_binned_df


__all__ = [
    "INDEL_SIZE_BINS",
    "INSERTIONS_KEY",
    "DELETIONS_KEY",
    "categorize_indel_size",
    "indel_stats_by_bin",
    "summarize_indel_sizes",
    "aggregate_indel_statistics",
]
