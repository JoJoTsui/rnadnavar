#!/usr/bin/env python3
"""
Reusable Plotly visualizations for indel size distributions.
"""

from __future__ import annotations

from typing import Literal

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

INSERTIONS_COLOR = "#00CC96"
DELETIONS_COLOR = "#EF553B"


def _build_text_labels(counts, pcts):
    return [f"{c:,}<br>({p:.1f}%)" for c, p in zip(counts, pcts)]


def plot_indel_binned_distribution(binned_df: pd.DataFrame, title: str) -> go.Figure:
    """Grouped bar chart of aggregated indel bins (insertions vs deletions)."""
    bins = binned_df["Size_Bin"]
    ins_counts = binned_df["Insertions_Count"]
    del_counts = binned_df["Deletions_Count"]
    ins_pct = binned_df["Insertions_Pct"]
    del_pct = binned_df["Deletions_Pct"]

    fig = go.Figure()
    fig.add_bar(
        name="Insertions",
        x=bins,
        y=ins_counts,
        marker_color=INSERTIONS_COLOR,
        text=_build_text_labels(ins_counts, ins_pct),
        textposition="outside",
        offsetgroup=0,
    )
    fig.add_bar(
        name="Deletions",
        x=bins,
        y=del_counts,
        marker_color=DELETIONS_COLOR,
        text=_build_text_labels(del_counts, del_pct),
        textposition="outside",
        offsetgroup=1,
    )

    fig.update_layout(
        title=title,
        xaxis=dict(title="Size Bin (bp)", tickmode="linear"),
        yaxis=dict(title="Count"),
        barmode="group",
        height=600,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        hovermode="x unified",
    )
    return fig


def plot_indel_per_sample_heatmap(
    per_sample_binned_df: pd.DataFrame,
    value_col: str,
    title: str,
    color_scale: str = "Blues",
) -> go.Figure:
    """Heatmap of per-sample indel bin percentages."""
    if "Sample_Model" in per_sample_binned_df.columns:
        index_col = "Sample_Model"
    else:
        index_col = "Sample"

    matrix = per_sample_binned_df.pivot(
        index=index_col, columns="Size_Bin", values=value_col
    ).fillna(0)
    fig = px.imshow(
        matrix,
        color_continuous_scale=color_scale,
        aspect="auto",
        labels=dict(x="Size Bin", y="Sample", color=value_col.replace("_", " ")),
    )
    fig.update_layout(title=title, height=500)
    return fig


def plot_indel_per_sample_stacked(
    per_sample_binned_df: pd.DataFrame,
    variant_type: Literal["Insertions", "Deletions"] = "Insertions",
    value_col: str = "Insertions_Count",
    title: str = "",
    color_map: dict | None = None,
    height: int = 550,
) -> go.Figure:
    """Stacked bar chart by sample for a given variant type (counts or percentages)."""
    df = per_sample_binned_df.copy()
    if "Sample_Model" not in df.columns:
        df["Sample_Model"] = df["Sample"]

    color_map = color_map or {
        "1": "#c7e9c0",
        "2": "#a1d99b",
        "3": "#74c476",
        "4-14": "#41ab5d",
        "15-29": "#238b45",
        "30-49": "#006d2c",
        ">50": "#00441b",
    }

    value_column = value_col
    if variant_type == "Deletions" and value_col.startswith("Insertions"):
        value_column = value_col.replace("Insertions", "Deletions")
    if variant_type == "Deletions" and value_col.startswith("Total_Insertions"):
        value_column = value_col.replace("Insertions", "Deletions")

    filtered = df[["Sample_Model", "Size_Bin", value_column]].rename(
        columns={value_column: "Value"}
    )
    fig = px.bar(
        filtered,
        x="Sample_Model",
        y="Value",
        color="Size_Bin",
        title=title,
        color_discrete_map=color_map,
        text="Value",
    )
    fig.update_traces(
        texttemplate="%{text}",
        textposition="outside",
        hovertemplate="%{x}<br>%{color}: %{y}",
    )
    fig.update_layout(
        barmode="stack", height=height, xaxis_tickangle=45, showlegend=True
    )
    fig.update_yaxes(title_text=value_column.replace("_", " "))
    return fig


__all__ = [
    "INSERTIONS_COLOR",
    "DELETIONS_COLOR",
    "plot_indel_binned_distribution",
    "plot_indel_per_sample_heatmap",
    "plot_indel_per_sample_stacked",
]
