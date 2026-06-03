#!/usr/bin/env python3
"""
Visualize seq2neo FASTQ statistics using seaborn/matplotlib.

Reads sample_stats.tsv and generates 6 figures in the same directory.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

HERE = Path(__file__).resolve().parent
STATS_DIR = HERE / "stats"
SAMPLE_TSV = STATS_DIR / "sample_stats.tsv"

MODALITIES = ["DN", "DT", "RT"]
MOD_COLORS = {"DN": "#4C72B0", "DT": "#DD8452", "RT": "#55A868"}
MOD_LABELS = {"DN": "DNA Normal", "DT": "DNA Tumor", "RT": "RNA Tumor"}
SET_LABELS = {
    "1": "Set 1: Colorectal",
    "2": "Set 2: Colon",
    "3": "Set 3: Bile/Chol/Esoph/Mel",
    "4": "Set 4: Gastric/Lung/Panc/Rect",
}

sns.set_theme(style="whitegrid", context="notebook")
plt.rcParams["figure.dpi"] = 150


def load_and_prepare():
    """Load sample stats, compute primary metrics, build long-form DataFrames."""
    df = pd.read_csv(SAMPLE_TSV, sep="\t")

    # Long-form for primary data (one row per sample per modality)
    records = []
    for _, row in df.iterrows():
        for mod in MODALITIES:
            reads = row.get(f"{mod}_r1_reads")
            if pd.isna(reads) or reads == 0:
                continue
            bases = row[f"{mod}_r1_bases"] + row[f"{mod}_r2_bases"]
            gc = (row[f"{mod}_r1_GC"] + row[f"{mod}_r2_GC"]) / 2
            q20 = (row[f"{mod}_r1_Q20"] + row[f"{mod}_r2_Q20"]) / 2
            q30 = (row[f"{mod}_r1_Q30"] + row[f"{mod}_r2_Q30"]) / 2
            avgqual = (row[f"{mod}_r1_AvgQual"] + row[f"{mod}_r2_AvgQual"]) / 2
            records.append({
                "set": str(int(row["set"])),
                "patient_id": row["patient_id"],
                "disease": row["disease"],
                "status": row["status"],
                "modality": mod,
                "reads": reads,
                "bases": bases,
                "GC": gc,
                "Q20": q20,
                "Q30": q30,
                "AvgQual": avgqual,
            })
    primary = pd.DataFrame(records)

    # Long-form for extra data (DT_extra, RT_extra only)
    extra_records = []
    for _, row in df.iterrows():
        if row["n_extra_pairs"] == 0:
            continue
        for mod in ["DT", "RT"]:
            reads = row.get(f"{mod}_extra_reads")
            if pd.isna(reads) or reads == 0:
                continue
            extra_records.append({
                "set": str(int(row["set"])),
                "patient_id": row["patient_id"],
                "status": row["status"],
                "modality": mod,
                "reads": reads,
                "bases": row[f"{mod}_extra_bases"],
                "GC": row[f"{mod}_extra_GC"],
                "Q20": row[f"{mod}_extra_Q20"],
                "Q30": row[f"{mod}_extra_Q30"],
                "AvgQual": row[f"{mod}_extra_AvgQual"],
            })
    extra = pd.DataFrame(extra_records)

    # Per-sample totals for sorting across all grouped-bar figures
    df["DN_total"] = df["DN_r1_reads"] + df["DN_r2_reads"]
    df["DT_total"] = df["DT_r1_reads"] + df["DT_r2_reads"]
    df["RT_total"] = df["RT_r1_reads"] + df["RT_r2_reads"]
    df["total_reads_sort"] = df["DN_total"] + df["DT_total"] + df["RT_total"]
    df["set_str"] = df["set"].astype(int).astype(str)

    # Pre-compute per-modality metrics for grouped-bar figures
    for mod in MODALITIES:
        df[f"{mod}_reads"] = df[f"{mod}_r1_reads"] + df[f"{mod}_r2_reads"]
        df[f"{mod}_bases"] = df[f"{mod}_r1_bases"] + df[f"{mod}_r2_bases"]
        df[f"{mod}_Q20"] = (df[f"{mod}_r1_Q20"] + df[f"{mod}_r2_Q20"]) / 2
        df[f"{mod}_Q30"] = (df[f"{mod}_r1_Q30"] + df[f"{mod}_r2_Q30"]) / 2

    return df, primary, extra


# ---------------------------------------------------------------------------
# Shared grouped-bar helper
# ---------------------------------------------------------------------------

def _grouped_bar_by_set(df, value_col, title, filename, xlabel):
    """4-panel horizontal grouped bar chart with unified x-axis across panels."""
    fig, axes = plt.subplots(2, 2, figsize=(20, 20))
    max_val = 0

    for ax, (set_id, set_label) in zip(axes.flat, SET_LABELS.items()):
        subset = df[df["set_str"] == set_id].sort_values("total_reads_sort")
        n = len(subset)
        if n == 0:
            ax.set_visible(False)
            continue

        y = np.arange(n)
        bar_height = 0.8 / len(MODALITIES)

        for i, mod in enumerate(MODALITIES):
            vals = subset[f"{mod}_{value_col}"].values
            offset = (i - (len(MODALITIES) - 1) / 2) * bar_height
            ax.barh(y + offset, vals, bar_height,
                    label=MOD_LABELS[mod], color=MOD_COLORS[mod])

        ax.set_yticks(y)
        ax.set_yticklabels(subset["patient_id"].values, fontsize=7)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_title(set_label, fontweight="bold", fontsize=11)
        ax.legend(loc="lower right", fontsize=8)
        max_val = max(max_val, subset[[f"{m}_{value_col}" for m in MODALITIES]].max().max())

    # Unified x-axis across all panels
    for ax in axes.flat:
        if ax.get_visible():
            ax.set_xlim(0, max_val * 1.08)

    fig.suptitle(title, fontweight="bold", fontsize=15, y=1.01)
    fig.tight_layout()
    fig.savefig(STATS_DIR / filename, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {filename}")


# ---------------------------------------------------------------------------
# FIG 1: Primary stats distribution (violin)
# ---------------------------------------------------------------------------

def fig1_primary_overview(primary):
    """2x3 grid of violin plots: DN/DT/RT distributions across all samples."""
    metrics = ["reads", "bases", "GC", "Q20", "Q30", "AvgQual"]
    titles = ["Reads", "Bases", "GC (%)", "Q20 (%)", "Q30 (%)", "AvgQual"]
    log_scale = {"reads", "bases"}

    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    for ax, metric, title in zip(axes.flat, metrics, titles):
        data = primary[["modality", metric]].dropna()
        sns.violinplot(data=data, x="modality", y=metric, hue="modality",
                       palette=MOD_COLORS, ax=ax, linewidth=1,
                       inner="quartile", density_norm="width", legend=False)
        if metric in log_scale:
            ax.set_yscale("log")
        ax.set_xlabel("")
        ax.set_title(title, fontweight="bold", fontsize=12)
        ax.set_xticks(range(len(MODALITIES)))
        ax.set_xticklabels([MOD_LABELS[m] for m in MODALITIES],
                           rotation=15, fontsize=9)
        ax.tick_params(axis="y", labelsize=9)

    fig.suptitle("Primary FASTQ Stats Distribution Across 66 Samples",
                 fontweight="bold", fontsize=16, y=1.01)
    fig.tight_layout()
    fig.savefig(STATS_DIR / "fig1_primary_overview.png",
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  fig1_primary_overview.png")


# ---------------------------------------------------------------------------
# FIG 2–5: Per-sample grouped-bar plots (reads, Q30, Q20, bases)
# ---------------------------------------------------------------------------

def fig2_per_sample_reads(df):
    _grouped_bar_by_set(df, value_col="reads",
        title="Per-Sample Sequencing Reads by Modality and Set",
        filename="fig2_per_sample_reads.png",
        xlabel="Reads (r1 + r2)")


def fig3_per_sample_q30(df):
    _grouped_bar_by_set(df, value_col="Q30",
        title="Per-Sample Q30 (%) by Modality and Set",
        filename="fig3_per_sample_q30.png",
        xlabel="Q30 (%)")


def fig4_per_sample_q20(df):
    _grouped_bar_by_set(df, value_col="Q20",
        title="Per-Sample Q20 (%) by Modality and Set",
        filename="fig4_per_sample_q20.png",
        xlabel="Q20 (%)")


def fig5_per_sample_bases(df):
    _grouped_bar_by_set(df, value_col="bases",
        title="Per-Sample Bases by Modality and Set",
        filename="fig5_per_sample_bases.png",
        xlabel="Bases (r1 + r2)")


# ---------------------------------------------------------------------------
# FIG 6: Extra data overview (violin)
# ---------------------------------------------------------------------------

def fig6_extra_overview(extra):
    """Violin plots for extra data (DT and RT) across samples with extra pairs."""
    if extra.empty:
        print("  fig6_extra_overview.png (skipped: no extra data)")
        return

    metrics = ["reads", "bases", "GC", "Q20", "Q30", "AvgQual"]
    titles = ["Extra Reads", "Extra Bases", "Extra GC (%)",
              "Extra Q20 (%)", "Extra Q30 (%)", "Extra AvgQual"]
    log_scale = {"reads", "bases"}
    n_samples = len(extra["patient_id"].unique())

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    for ax, metric, title in zip(axes.flat, metrics, titles):
        data = extra[["modality", metric]].dropna()
        sns.violinplot(data=data, x="modality", y=metric, hue="modality",
                       palette=MOD_COLORS, ax=ax, linewidth=1,
                       inner="quartile", density_norm="width", legend=False)
        if metric in log_scale:
            ax.set_yscale("log")
        ax.set_xlabel("")
        ax.set_title(title, fontweight="bold", fontsize=12)
        ax.set_xticks(range(2))
        ax.set_xticklabels([f"{MOD_LABELS[m]} (extra)" for m in ["DT", "RT"]],
                           fontsize=9)
        ax.tick_params(axis="y", labelsize=9)

    fig.suptitle(f"Extra FASTQ Stats Distribution ({n_samples} samples)",
                 fontweight="bold", fontsize=16, y=1.01)
    fig.tight_layout()
    fig.savefig(STATS_DIR / "fig6_extra_overview.png",
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  fig6_extra_overview.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading and preparing data...")
    df, primary, extra = load_and_prepare()
    print(f"  {len(primary)} primary records ({len(primary['patient_id'].unique())} samples)")
    print(f"  {len(extra)} extra records ({len(extra['patient_id'].unique())} samples)")

    print("Generating figures...")
    fig1_primary_overview(primary)
    fig2_per_sample_reads(df)
    fig3_per_sample_q30(df)
    fig4_per_sample_q20(df)
    fig5_per_sample_bases(df)
    fig6_extra_overview(extra)

    print(f"\nAll figures saved to {STATS_DIR}/")


if __name__ == "__main__":
    main()
