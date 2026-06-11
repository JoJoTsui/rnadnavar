"""Generate interactive Altair charts and static image exports.

Creates 13 chart types from seq2neo variant statistics. Each chart is
exported as interactive HTML (embeds vega-embed), PNG (via vl-convert),
and SVG (via vl-convert). No Chrome dependency.
"""

from pathlib import Path

import altair as alt
import polars as pl

# Disable altair's 5000-row default limit
alt.data_transformers.disable_max_rows()

def _maybe_collect(df):
    """Collect if lazy, pass-through if eager."""
    if isinstance(df, pl.LazyFrame):
        return df.collect()
    return df


def _maybe_collect(df):
    """Collect if lazy, pass-through if eager."""
    if isinstance(df, pl.LazyFrame):
        return df.collect()
    return df



# ── Color palettes ────────────────────────────────────────────────────────
VC_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
VC_DOMAIN = ["Somatic", "Germline", "Reference", "Artifact"]


def _save_chart(chart: alt.Chart, name: str, output_dir: str):
    """Save a chart as HTML, PNG, and SVG.

    HTML: embeds vega-embed inline (self-contained).
    PNG/SVG: via vl-convert-python (no Chrome needed).
    """
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    html_path = plots_dir / f"{name}.html"
    chart.save(str(html_path))

    png_path = plots_dir / f"{name}.png"
    chart.save(str(png_path), format="png", scale_factor=3)

    svg_path = plots_dir / f"{name}.svg"
    chart.save(str(svg_path), format="svg")


# ── Chart functions ───────────────────────────────────────────────────────


def plot_vc_distribution(df, output_dir: str, group_col: str = "set_number"):
    """Chart 1: Variant counts by VC classification, stacked bar per set."""
    if "VC" not in df.columns:
        return
    counts = (
        df.group_by([group_col, "VC"]).agg(pl.len().alias("count"))
        .sort([group_col, "VC"]).pipe(_maybe_collect).to_pandas()
    )
    counts[group_col] = counts[group_col].astype(str)
    chart = alt.Chart(counts).mark_bar().encode(
        x=alt.X(f"{group_col}:N", title="Set"),
        y=alt.Y("count:Q", title="Number of Variants"),
        color=alt.Color(
            "VC:N",
            title="Variant Classification",
            scale=alt.Scale(domain=VC_DOMAIN, range=VC_COLORS),
            legend=alt.Legend(
                orient="right",
                title="Variant Classification",
                labelFontSize=11,
                titleFontSize=12,
            ),
        ),
    ).properties(title="Variant Classification Distribution by Set")
    _save_chart(chart, "01_vc_distribution", output_dir)
    return chart


def plot_caller_overlap(df, output_dir: str):
    """Chart 2: Variant tier distribution (% by CxDy final tier)."""
    if "final_tier" not in df.columns or "set_number" not in df.columns:
        return
    counts = df.group_by(["set_number", "final_tier"]).agg(pl.len().alias("count"))
    total_per_set = df.group_by("set_number").agg(pl.len().alias("total"))
    counts = counts.join(total_per_set, on="set_number").with_columns(
        (pl.col("count") / pl.col("total") * 100).alias("pct")
    ).pipe(_maybe_collect).to_pandas()
    counts["set_number"] = counts["set_number"].astype(str)
    chart = alt.Chart(counts).mark_bar().encode(
        x=alt.X("final_tier:N", title="Variant Tier (CxDy)"),
        y=alt.Y("pct:Q", title="% of Variants"),
        color=alt.Color("set_number:N", title="Set"),
        column=alt.Column("set_number:N", title="Set"),
    ).properties(title="Variant Tier Distribution by Set")
    _save_chart(chart, "02_caller_overlap", output_dir)
    return chart


def plot_vaf_distribution(df, output_dir: str):
    """Chart 3: Per-caller VAF distribution, box plot."""
    caller_vaf_cols = [
        f"{c}_VAF" for c in [
            "DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic",
            "RNA_deepsomatic", "DNA_strelka", "RNA_strelka",
        ]
    ]
    existing = [c for c in caller_vaf_cols if c in df.columns]
    if not existing:
        return
    melted = df.select(existing).unpivot(
        variable_name="caller", value_name="VAF"
    ).drop_nulls()
    # Sample to avoid vl-convert buffer overflow on large datasets
    if melted.height > 50000:
        melted = melted.sample(50000)
    melted = melted.pipe(_maybe_collect).to_pandas()
    chart = alt.Chart(melted).mark_boxplot().encode(
        x=alt.X("caller:N", title="Caller", axis=alt.Axis(labelAngle=-45)),
        y=alt.Y("VAF:Q", title="Variant Allele Frequency"),
    ).properties(title="VAF Distribution per Caller")
    _save_chart(chart, "03_vaf_distribution", output_dir)
    return chart


def plot_dna_vs_rna_vaf(df, output_dir: str):
    """Chart 4: DNA vs RNA mean VAF scatter."""
    if "DNA_VAF_mean" not in df.columns or "RNA_VAF_mean" not in df.columns:
        return
    cols = ["DNA_VAF_mean", "RNA_VAF_mean"]
    if "VC" in df.columns:
        cols.append("VC")
    pdf = df.select(cols).drop_nulls(subset=["DNA_VAF_mean", "RNA_VAF_mean"])
    if pdf.height > 5000:
        pdf = pdf.sample(5000)
    pdf = pdf.pipe(_maybe_collect).to_pandas()
    color_enc = alt.Color("VC:N", scale=alt.Scale(domain=VC_DOMAIN, range=VC_COLORS)) if "VC" in pdf.columns else alt.value("#1f77b4")
    chart = alt.Chart(pdf).mark_circle(opacity=0.4, size=20).encode(
        x=alt.X("DNA_VAF_mean:Q", title="DNA Mean VAF"),
        y=alt.Y("RNA_VAF_mean:Q", title="RNA Mean VAF"),
        color=color_enc,
    ).properties(title="DNA vs RNA Mean VAF")
    _save_chart(chart, "04_dna_vs_rna_vaf", output_dir)
    return chart


def plot_dna_vs_rna_dp(df, output_dir: str):
    """Chart 5: DNA vs RNA mean DP scatter."""
    if "DNA_DP_mean" not in df.columns or "RNA_DP_mean" not in df.columns:
        return
    pdf = df.select(["DNA_DP_mean", "RNA_DP_mean", "set_number"]).drop_nulls()
    if pdf.height > 5000:
        pdf = pdf.sample(5000)
    pdf = pdf.pipe(_maybe_collect).to_pandas()
    pdf["set_number"] = pdf["set_number"].astype(str)
    chart = alt.Chart(pdf).mark_circle(opacity=0.4, size=20).encode(
        x=alt.X("DNA_DP_mean:Q", title="DNA Mean Depth"),
        y=alt.Y("RNA_DP_mean:Q", title="RNA Mean Depth"),
        color=alt.Color("set_number:N", title="Set"),
    ).properties(title="DNA vs RNA Mean Depth")
    _save_chart(chart, "05_dna_vs_rna_dp", output_dir)
    return chart


def plot_gt_concordance(df, output_dir: str):
    """Chart 6: GT concordance among 4 callers with GT fields."""
    gt_cols = [
        "DNA_mutect2_GT", "RNA_mutect2_GT",
        "DNA_deepsomatic_GT", "RNA_deepsomatic_GT",
    ]
    existing = [c for c in gt_cols if c in df.columns]
    if len(existing) < 2:
        return

    from collections import Counter

    # Compute GT agreement counts using iter_rows() (tuples, not dicts)
    # Avoids df.to_dicts() which creates a Python dict per row — O(n) memory
    agree_counts = {2: 0, 3: 0, 4: 0}
    n_total = 0
    for row in df.select(existing).iter_rows():
        gts = [g for g in row if g is not None and g not in ("./.", "./.", ".")]
        if len(gts) >= 2:
            best = Counter(gts).most_common(1)[0][1]
            if best >= 2:
                agree_counts[best] = agree_counts.get(best, 0) + 1
                n_total += 1
    if n_total == 0:
        return

    counts = pl.DataFrame({
        "agreement_level": ["2", "3", "4"],
        "count": [agree_counts[2], agree_counts[3], agree_counts[4]],
    }).pipe(_maybe_collect).to_pandas()

    chart = alt.Chart(counts).mark_bar().encode(
        x=alt.X("agreement_level:N", title="Number of Callers Agreeing on GT"),
        y=alt.Y("count:Q", title="Number of Variants"),
    ).properties(
        title=f"GT Concordance Among 4 Callers (n={n_total} with ≥2 agreement)"
    )
    _save_chart(chart, "06_gt_concordance", output_dir)
    return chart


def plot_cosmic_gnomad_annotation(df, output_dir: str):
    """Chart 7: COSMIC/gnomAD annotation coverage as side-by-side pies."""
    n = df.height
    has_cosmic = df.filter(pl.col("COSMIC_ID").is_not_null()).height if "COSMIC_ID" in df.columns else -1
    has_gnomad = df.filter(pl.col("GNOMAD_AF").is_not_null()).height if "GNOMAD_AF" in df.columns else -1
    if has_cosmic < 0 and has_gnomad < 0:
        return

    charts = []
    if has_cosmic >= 0:
        cd = pl.DataFrame({
            "category": ["In COSMIC", "Not in COSMIC"],
            "count": [has_cosmic, n - has_cosmic],
        }).pipe(_maybe_collect).to_pandas()
        cosmic_chart = alt.Chart(cd).mark_arc(innerRadius=40).encode(
            theta=alt.Theta("count:Q"),
            color=alt.Color("category:N", scale=alt.Scale(
                domain=["In COSMIC", "Not in COSMIC"], range=["#1f77b4", "#d3d3d3"]
            )),
        ).properties(title="COSMIC Annotation")
        charts.append(cosmic_chart)

    if has_gnomad >= 0:
        gd = pl.DataFrame({
            "category": ["Has gnomAD AF", "No gnomAD AF"],
            "count": [has_gnomad, n - has_gnomad],
        }).pipe(_maybe_collect).to_pandas()
        gnomad_chart = alt.Chart(gd).mark_arc(innerRadius=40).encode(
            theta=alt.Theta("count:Q"),
            color=alt.Color("category:N", scale=alt.Scale(
                domain=["Has gnomAD AF", "No gnomAD AF"], range=["#ff7f0e", "#d3d3d3"]
            )),
        ).properties(title="gnomAD Annotation")
        charts.append(gnomad_chart)

    if len(charts) == 2:
        chart = alt.hconcat(*charts).properties(
            title="COSMIC and gnomAD Annotation Coverage"
        )
    else:
        chart = charts[0]
    _save_chart(chart, "07_cosmic_gnomad", output_dir)
    return chart


def plot_variant_type_distribution(df, output_dir: str):
    """Chart 9: Variant type distribution (SNV/INS/DEL/MNV), stacked bar per set."""
    if "variant_type" not in df.columns:
        return
    pdf = df.group_by(["set_number", "variant_type"]).agg(
        pl.len().alias("count")
    ).pipe(_maybe_collect).to_pandas()
    pdf["set_number"] = pdf["set_number"].astype(str)
    chart = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("set_number:N", title="Set"),
        y=alt.Y("count:Q", title="Number of Variants"),
        color=alt.Color("variant_type:N", title="Variant Type"),
    ).properties(title="Variant Type Distribution by Set")
    _save_chart(chart, "09_variant_type_distribution", output_dir)
    return chart


def plot_ti_tv_ratio(df, output_dir: str):
    """Chart 10: Ti/Tv ratio bar chart per set."""
    if "ti_tv" not in df.columns:
        return
    ti = df.filter(pl.col("ti_tv") == True).group_by("set_number").agg(pl.len().alias("ti"))
    tv = df.filter(pl.col("ti_tv") == False).group_by("set_number").agg(pl.len().alias("tv"))
    ratio = ti.join(tv, on="set_number").with_columns(
        (pl.col("ti") / pl.col("tv")).alias("ratio")).pipe(_maybe_collect).to_pandas()
    ratio["set_number"] = ratio["set_number"].astype(str)
    bars = alt.Chart(ratio).mark_bar().encode(
        x=alt.X("set_number:N", title="Set"),
        y=alt.Y("ratio:Q", title="Ti/Tv Ratio"),
    )
    text = alt.Chart(ratio).mark_text(dy=-8).encode(
        x=alt.X("set_number:N"),
        y=alt.Y("ratio:Q"),
        text=alt.Text("ratio:Q", format=".2f"),
    )
    chart = (bars + text).properties(title="Ti/Tv Ratio by Set")
    _save_chart(chart, "10_ti_tv_ratio", output_dir)
    return chart


def plot_cross_modality(df, output_dir: str):
    """Chart 11: Cross-modality & rescue analysis with percentages."""
    cols_needed = ["CROSS_MODALITY", "RESCUED", "set_number"]
    if not all(c in df.columns for c in cols_needed):
        return

    # Compute total variants per set for percentage calculation
    total_per_set = df.group_by("set_number").agg(pl.len().alias("total"))

    subcharts = []
    for col in ["CROSS_MODALITY", "RESCUED"]:
        pdf = df.group_by(["set_number", col]).agg(
            pl.len().alias("count")
        ).join(total_per_set, on="set_number").with_columns(
            (pl.col("count") / pl.col("total") * 100).round(1).alias("pct")
        ).pipe(_maybe_collect).to_pandas()
        pdf["set_number"] = pdf["set_number"].astype(str)

        bars = alt.Chart(pdf).mark_bar().encode(
            x=alt.X("set_number:N", title="Set"),
            y=alt.Y("count:Q", title="Count"),
            color=alt.Color(f"{col}:N"),
        )
        text = alt.Chart(pdf).mark_text(dy=-8, fontSize=10).encode(
            x=alt.X("set_number:N"),
            y=alt.Y("count:Q"),
            text=alt.Text("pct:N", format=".1f"),
            color=alt.Color(f"{col}:N"),
        ).transform_filter(alt.datum["count"] > 0)
        c = (bars + text).properties(
            title=f"{'Cross-Modality' if col == 'CROSS_MODALITY' else 'Rescued Variants'}"
        )
        subcharts.append(c)

    chart = alt.hconcat(*subcharts).properties(
        title="Cross-Modality and Rescue Analysis (%)"
    )
    _save_chart(chart, "11_cross_modality", output_dir)
    return chart


def plot_per_sample_distribution(sample_stats_df, output_dir: str, top_n: int = 30):
    """Chart 8: Per-sample variant count, horizontal bar chart. Top-N by count."""
    if sample_stats_df is None or (hasattr(sample_stats_df, 'is_empty') and sample_stats_df.is_empty()):
        return
    if "total_variants" not in sample_stats_df.columns:
        return
    if "sample_id" not in sample_stats_df.columns:
        return

    pdf = sample_stats_df.select(["sample_id", "total_variants"]).sort(
        "total_variants", descending=True
    ).head(top_n).pipe(_maybe_collect).to_pandas()
    chart = alt.Chart(pdf).mark_bar().encode(
        y=alt.Y("sample_id:N", title="Sample ID", sort=None),
        x=alt.X("total_variants:Q", title="Total Variants per Sample"),
        tooltip=["sample_id", "total_variants"],
    ).properties(title=f"Per-Sample Variant Counts (Top {top_n})")
    _save_chart(chart, "08_per_sample_distribution", output_dir)
    return chart


# Alias for backward compatibility
plot_per_sample_violin = plot_per_sample_distribution


def plot_validation_heatmap(report, output_dir: str):
    """Chart 13: Rescue VCF validation heatmap (mismatch % per metric × sample)."""
    if report is None or (hasattr(report, 'is_empty') and report.is_empty()):
        return
    pivot = report.pipe(_maybe_collect).to_pandas().pivot(
        index="sample_id", columns="metric", values="mismatch_pct"
    ).reset_index().melt(id_vars="sample_id", var_name="metric", value_name="mismatch_pct")
    chart = alt.Chart(pivot).mark_rect().encode(
        x=alt.X("metric:N", title="Validation Metric"),
        y=alt.Y("sample_id:N", title="Sample"),
        color=alt.Color("mismatch_pct:Q", title="Mismatch %",
                        scale=alt.Scale(scheme="redyellowgreen", reverse=True)),
    ).properties(title="Rescue VCF Validation: Mismatch % by Metric × Sample")
    _save_chart(chart, "13_validation_heatmap", output_dir)
    return chart


# ── Tier-aware charts ────────────────────────────────────────────────────────


def plot_vaf_boxplot_per_tier(df, output_dir: str):
    """VAF distribution per caller, boxplot faceted by caller tier."""
    caller_vaf_cols = [
        f"{c}_VAF" for c in [
            "DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic",
            "RNA_deepsomatic", "DNA_strelka", "RNA_strelka",
        ]
    ]
    existing = [c for c in caller_vaf_cols if c in df.columns]
    if not existing or "caller_tier" not in df.columns:
        return

    melted = df.select(existing + ["caller_tier"]).unpivot(
        index=["caller_tier"], variable_name="caller", value_name="VAF"
    ).drop_nulls()
    if melted.height > 50000:
        melted = melted.sample(50000)
    pdf = melted.pipe(_maybe_collect).to_pandas()
    # Clean caller names for display
    pdf["caller"] = pdf["caller"].str.replace("_VAF", "")

    chart = alt.Chart(pdf).mark_boxplot(extent="min-max").encode(
        x=alt.X("caller:N", title="Caller", axis=alt.Axis(labelAngle=-45)),
        y=alt.Y("VAF:Q", title="Variant Allele Frequency"),
        color=alt.Color("caller:N"),
        column=alt.Column("caller_tier:N", title="Caller Tier"),
    ).properties(title="VAF Distribution per Caller × Caller Tier")
    _save_chart(chart, "14_vaf_violin_per_tier", output_dir)
    return chart


def plot_dp_boxplot_per_tier(df, output_dir: str):
    """DP distribution per caller, boxplot faceted by caller tier."""
    caller_dp_cols = [
        f"{c}_DP" for c in [
            "DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic",
            "RNA_deepsomatic", "DNA_strelka", "RNA_strelka",
        ]
    ]
    existing = [c for c in caller_dp_cols if c in df.columns]
    if not existing or "caller_tier" not in df.columns:
        return

    melted = df.select(existing + ["caller_tier"]).unpivot(
        index=["caller_tier"], variable_name="caller", value_name="DP"
    ).drop_nulls()
    if melted.height > 50000:
        melted = melted.sample(50000)
    pdf = melted.pipe(_maybe_collect).to_pandas()
    pdf["caller"] = pdf["caller"].str.replace("_DP", "")

    chart = alt.Chart(pdf).mark_boxplot(extent="min-max").encode(
        x=alt.X("caller:N", title="Caller", axis=alt.Axis(labelAngle=-45)),
        y=alt.Y("DP:Q", title="Read Depth"),
        color=alt.Color("caller:N"),
        column=alt.Column("caller_tier:N", title="Caller Tier"),
    ).properties(title="DP Distribution per Caller × Caller Tier")
    _save_chart(chart, "15_dp_per_tier", output_dir)
    return chart


def plot_gt_concordance_per_tier(df, output_dir: str):
    """GT concordance faceted by caller tier."""
    gt_cols = [
        "DNA_mutect2_GT", "RNA_mutect2_GT",
        "DNA_deepsomatic_GT", "RNA_deepsomatic_GT",
    ]
    existing_gt = [c for c in gt_cols if c in df.columns]
    if len(existing_gt) < 2 or "caller_tier" not in df.columns:
        return

    from collections import Counter

    tiers = df["caller_tier"].unique().to_list()
    # Use iter_rows() not to_dicts() — avoid creating Python dicts per row
    agree_counts: dict[tuple, int] = {}
    for row in df.select(existing_gt + ["caller_tier"]).iter_rows():
        gts = [g for g in row[:-1]  # last column is caller_tier
               if g is not None and g not in ("./.", "./.", ".")]
        if len(gts) < 2:
            continue
        best = Counter(gts).most_common(1)[0][1]
        if best >= 2:
            key = (row[-1], best)  # (caller_tier, agreement)
            agree_counts[key] = agree_counts.get(key, 0) + 1

    if not agree_counts:
        return

    pdf = pl.DataFrame(
        [{"caller_tier": k[0], "agreement": k[1], "count": v} for k, v in agree_counts.items()]
    ).pipe(_maybe_collect).to_pandas()
    pdf["agreement"] = pdf["agreement"].astype(str)

    chart = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("agreement:N", title="Number of Callers Agreeing"),
        y=alt.Y("count:Q", title="Count"),
        color=alt.Color("agreement:N"),
        column=alt.Column("caller_tier:N", title="Caller Tier"),
    ).properties(title="GT Concordance per Caller Tier")
    _save_chart(chart, "16_gt_concordance_per_tier", output_dir)
    return chart


def plot_tiered_caller_overlap(df, output_dir: str):
    """N_SUPPORT_CALLERS histogram faceted by caller tier."""
    if "N_SUPPORT_CALLERS" not in df.columns or "caller_tier" not in df.columns:
        return

    counts = df.group_by(["caller_tier", "N_SUPPORT_CALLERS"]).agg(
        pl.len().alias("count")
    )
    total_per_tier = df.group_by("caller_tier").agg(pl.len().alias("total"))
    pdf = counts.join(total_per_tier, on="caller_tier").with_columns(
        (pl.col("count") / pl.col("total") * 100).alias("pct")
    ).pipe(_maybe_collect).to_pandas()

    chart = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("N_SUPPORT_CALLERS:O", title="Number of Supporting Callers"),
        y=alt.Y("pct:Q", title="% of Variants"),
        color=alt.Color("caller_tier:N", title="Caller Tier"),
        column=alt.Column("caller_tier:N", title="Caller Tier"),
    ).properties(title="Caller Support Distribution per Caller Tier")
    _save_chart(chart, "18_caller_overlap_per_tier", output_dir)
    return chart


def plot_tiered_variant_types(df, output_dir: str):
    """Variant type distribution (SNV/INS/DEL/MNV) faceted by caller tier."""
    if "variant_type" not in df.columns or "caller_tier" not in df.columns:
        return

    pdf = df.group_by(["caller_tier", "variant_type"]).agg(
        pl.len().alias("count")
    ).pipe(_maybe_collect).to_pandas()

    chart = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("variant_type:N", title="Variant Type"),
        y=alt.Y("count:Q", title="Number of Variants"),
        color=alt.Color("variant_type:N", title="Variant Type"),
        column=alt.Column("caller_tier:N", title="Caller Tier"),
    ).properties(title="Variant Type Distribution per Caller Tier")
    _save_chart(chart, "19_variant_types_per_tier", output_dir)
    return chart


def plot_ref_alt_dp_scatter(df, output_dir: str):
    """Chart: DNA vs RNA mean REF_DP and ALT_DP scatter plots."""
    needed = ["DNA_REF_DP_mean", "RNA_REF_DP_mean",
              "DNA_ALT_DP_mean", "RNA_ALT_DP_mean"]
    if not all(c in df.columns for c in needed):
        return

    subcharts = []
    for label, x_col, y_col in [
        ("REF_DP", "DNA_REF_DP_mean", "RNA_REF_DP_mean"),
        ("ALT_DP", "DNA_ALT_DP_mean", "RNA_ALT_DP_mean"),
    ]:
        pdf = df.select([x_col, y_col, "set_number"]).drop_nulls()
        if pdf.height > 5000:
            pdf = pdf.sample(5000)
        pdf_pd = pdf.pipe(_maybe_collect).to_pandas()
        pdf_pd["set_number"] = pdf_pd["set_number"].astype(str)
        c = alt.Chart(pdf_pd).mark_circle(opacity=0.4, size=20).encode(
            x=alt.X(f"{x_col}:Q", title=f"DNA Mean {label}"),
            y=alt.Y(f"{y_col}:Q", title=f"RNA Mean {label}"),
            color=alt.Color("set_number:N", title="Set"),
        ).properties(title=f"DNA vs RNA Mean {label}")
        subcharts.append(c)

    chart = alt.hconcat(*subcharts).properties(
        title="DNA vs RNA REF_DP and ALT_DP"
    )
    _save_chart(chart, "17_ref_alt_dp_scatter", output_dir)
    return chart


# ── BAM statistics charts ────────────────────────────────────────────────────


def plot_bam_metrics_bars(bam_stats_df, output_dir: str, top_n: int = 20):
    """BAM metrics: grouped bar chart of per-sample reads for DN/DT/RT (top-N samples)."""
    if bam_stats_df is None or (hasattr(bam_stats_df, 'is_empty') and bam_stats_df.is_empty()):
        return
    needed = ["sample_id", "bam_type", "total_reads"]
    if not all(c in bam_stats_df.columns for c in needed):
        return
    # Get top-N samples by total reads
    top_ids = bam_stats_df.group_by("sample_id").agg(
        pl.col("total_reads").max().alias("max_reads")
    ).sort("max_reads", descending=True).head(top_n)["sample_id"].to_list()
    pdf = bam_stats_df.filter(pl.col("sample_id").is_in(top_ids)).select(
        ["sample_id", "bam_type", "total_reads", "mapped_reads"]
    ).pipe(_maybe_collect).to_pandas()
    chart = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("sample_id:N", title="Sample", axis=alt.Axis(labelAngle=-45)),
        y=alt.Y("total_reads:Q", title="Total Reads"),
        color=alt.Color("bam_type:N", title="BAM Type"),
        xOffset="bam_type:N",
    ).properties(title=f"Per-Sample BAM Read Counts — DN/DT/RT (Top {top_n})")
    _save_chart(chart, "20_bam_metrics", output_dir)
    return chart


def plot_bam_coverage_violin(df, output_dir: str):
    """BAM coverage distribution per BAM type, boxplot by tier."""
    bam_cols = [c for c in df.columns if c.startswith("BAM_DP_")]
    if not bam_cols or "caller_tier" not in df.columns:
        return
    melted = df.select(bam_cols + ["caller_tier"]).unpivot(
        index=["caller_tier"], variable_name="bam_type", value_name="DP"
    ).drop_nulls()
    if melted.height > 50000:
        melted = melted.sample(50000)
    pdf = melted.pipe(_maybe_collect).to_pandas()
    chart = alt.Chart(pdf).mark_boxplot(extent="min-max").encode(
        x=alt.X("bam_type:N", title="BAM Type"),
        y=alt.Y("DP:Q", title="Depth at Variant Position"),
        color=alt.Color("bam_type:N"),
        column=alt.Column("caller_tier:N", title="Caller Tier"),
    ).properties(title="BAM Coverage at Variant Positions per Tier")
    _save_chart(chart, "21_bam_coverage", output_dir)
    return chart


def plot_per_sample_tier_distribution(sample_tier_df, output_dir: str):
    """Stacked bar chart: per-sample per-tier variant counts."""
    if sample_tier_df is None or (hasattr(sample_tier_df, 'is_empty') and sample_tier_df.is_empty()):
        return
    if "sample_id" not in sample_tier_df.columns or "final_tier" not in sample_tier_df.columns:
        return
    pdf = sample_tier_df.select(["sample_id", "final_tier", "n_variants"]).pipe(_maybe_collect).to_pandas()
    chart = alt.Chart(pdf).mark_bar().encode(
        y=alt.Y("sample_id:N", title="Sample ID", sort=None),
        x=alt.X("n_variants:Q", title="Variants"),
        color=alt.Color("final_tier:N", title="Tier"),
    ).properties(title="Per-Sample Per-Tier Variant Distribution")
    _save_chart(chart, "22_sample_tier_dist", output_dir)
    return chart


def plot_per_tier_cross_sample_vaf(sample_tier_df, output_dir: str):
    """Boxplot: per-tier VAF distribution across samples."""
    if sample_tier_df is None or (hasattr(sample_tier_df, 'is_empty') and sample_tier_df.is_empty()):
        return
    vaf_col = None
    for c in ["mean_dna_vaf_mean", "mean_rna_vaf_mean"]:
        if c in sample_tier_df.columns:
            vaf_col = c
            break
    if not vaf_col or "final_tier" not in sample_tier_df.columns:
        return
    pdf = sample_tier_df.select(["final_tier", vaf_col]).pipe(_maybe_collect).to_pandas()
    chart = alt.Chart(pdf).mark_boxplot().encode(
        x=alt.X("final_tier:N", title="Tier"),
        y=alt.Y(f"{vaf_col}:Q", title=f"Mean {vaf_col}"),
        color=alt.Color("final_tier:N"),
    ).properties(title="Per-Tier VAF Distribution Across Samples")
    _save_chart(chart, "23_per_tier_vaf", output_dir)
    return chart


def plot_filter_distribution(df, output_dir: str):
    """FILTER value distribution per set (pie or bar)."""
    if "FILTER" not in df.columns:
        return
    counts = df.group_by("FILTER").agg(pl.len().alias("count")).sort("count", descending=True).pipe(_maybe_collect).to_pandas()
    chart = alt.Chart(counts).mark_bar().encode(
        x=alt.X("count:Q", title="Number of Variants"),
        y=alt.Y("FILTER:N", title="FILTER", sort="-x"),
        color=alt.Color("FILTER:N"),
    ).properties(title="FILTER Distribution")
    _save_chart(chart, "24_filter_distribution", output_dir)
    return chart


def plot_dna_vs_rna_per_caller(df, output_dir: str):
    """DNA vs RNA per-caller VAF scatter at shared positions (cross-modality)."""
    pairs = [("DNA_mutect2", "RNA_mutect2"), ("DNA_deepsomatic", "RNA_deepsomatic"),
             ("DNA_strelka", "RNA_strelka")]
    subcharts = []
    for dna_caller, rna_caller in pairs:
        dna_vaf = f"{dna_caller}_VAF"
        rna_vaf = f"{rna_caller}_VAF"
        if dna_vaf not in df.columns or rna_vaf not in df.columns:
            continue
        cols = [dna_vaf, rna_vaf]
        if "VC" in df.columns:
            cols.append("VC")
        pdf = df.select(cols).drop_nulls(subset=[dna_vaf, rna_vaf])
        if pdf.height > 5000:
            pdf = pdf.sample(5000)
        pdf_pd = pdf.pipe(_maybe_collect).to_pandas()
        caller_label = dna_caller.replace("DNA_", "")
        color_enc = alt.Color("VC:N", scale=alt.Scale(domain=VC_DOMAIN, range=VC_COLORS)) if "VC" in pdf_pd.columns else alt.value("#1f77b4")
        c = alt.Chart(pdf_pd).mark_circle(opacity=0.4, size=20).encode(
            x=alt.X(f"{dna_vaf}:Q", title=f"{caller_label} DNA VAF"),
            y=alt.Y(f"{rna_vaf}:Q", title=f"{caller_label} RNA VAF"),
            color=color_enc,
        ).properties(title=f"{caller_label}: DNA vs RNA VAF")
        subcharts.append(c)
    if not subcharts:
        return
    chart = alt.hconcat(*subcharts).properties(title="DNA vs RNA Per-Caller VAF")
    _save_chart(chart, "26_dna_vs_rna_per_caller", output_dir)
    return chart


def plot_chromosome_density(df, output_dir: str):
    """Variant count per chromosome (Manhattan-style bar chart)."""
    if "CHROM" not in df.columns:
        return
    counts = df.group_by("CHROM").agg(pl.len().alias("count")).sort("count", descending=True).pipe(_maybe_collect).to_pandas()
    chart = alt.Chart(counts).mark_bar().encode(
        x=alt.X("CHROM:N", title="Chromosome", sort="-y"),
        y=alt.Y("count:Q", title="Number of Variants"),
        tooltip=["CHROM", "count"],
    ).properties(title="Variant Density per Chromosome")
    _save_chart(chart, "30_chromosome_density", output_dir)
    return chart


def plot_redi_evidence(df, output_dir: str):
    """REDIportal RNA editing evidence distribution."""
    if "REDI_EVIDENCE" not in df.columns:
        return
    counts = df.group_by("REDI_EVIDENCE").agg(pl.len().alias("count")).sort("count", descending=True).pipe(_maybe_collect).to_pandas()
    chart = alt.Chart(counts).mark_bar().encode(
        x=alt.X("REDI_EVIDENCE:N", title="REDIportal Evidence Level"),
        y=alt.Y("count:Q", title="Number of Variants"),
        color=alt.Color("REDI_EVIDENCE:N"),
    ).properties(title="REDIportal RNA Editing Evidence")
    _save_chart(chart, "29_redi_evidence", output_dir)
    return chart


def plot_tier_quality_distribution(df, output_dir: str):
    """Tier quality score histogram."""
    if "tier_quality" not in df.columns:
        return
    pdf = df.select(["tier_quality"]).pipe(_maybe_collect).to_pandas()
    chart = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("tier_quality:Q", bin=alt.Bin(maxbins=20), title="Tier Quality Score"),
        y=alt.Y("count()", title="Number of Variants"),
    ).properties(title="Tier Quality Score Distribution")
    _save_chart(chart, "27_tier_quality", output_dir)
    return chart


def plot_per_tier_vaf_boxplot(df, output_dir: str):
    """Per-tier DNA VAF boxplot across all variants."""
    if "DNA_VAF_mean" not in df.columns or "final_tier" not in df.columns:
        return
    pdf = df.select(["final_tier", "DNA_VAF_mean", "RNA_VAF_mean"]).drop_nulls(subset=["DNA_VAF_mean"])
    if pdf.height > 50000:
        pdf = pdf.sample(50000)
    pdf = pdf.pipe(_maybe_collect).to_pandas()
    chart = alt.Chart(pdf).mark_boxplot().encode(
        x=alt.X("final_tier:N", title="Tier"),
        y=alt.Y("DNA_VAF_mean:Q", title="DNA Mean VAF"),
        color=alt.Color("final_tier:N"),
    ).properties(title="DNA VAF Distribution per Tier")
    _save_chart(chart, "25_per_tier_vaf", output_dir)
    return chart


def plot_caller_agreement_matrix(df, output_dir: str):
    """6×6 pairwise caller agreement matrix heatmap."""
    callers = ["DNA_mutect2", "DNA_deepsomatic", "DNA_strelka",
               "RNA_mutect2", "RNA_deepsomatic", "RNA_strelka"]
    gt_cols = [f"{c}_GT" for c in callers]
    existing = [c for c in gt_cols if c in df.columns]
    if len(existing) < 2:
        return
    # Compute pairwise agreement: for each pair, fraction where both have non-null GT
    n = df.height
    rows = []
    for c1 in callers:
        col1 = f"{c1}_GT"
        if col1 not in df.columns:
            continue
        for c2 in callers:
            col2 = f"{c2}_GT"
            if col2 not in df.columns:
                continue
            # Both have valid (non-null, non-./.) GT
            both_valid = df.filter(
                pl.col(col1).is_not_null() & pl.col(col2).is_not_null()
                & ~pl.col(col1).is_in(["./.", "./.", "."])
                & ~pl.col(col2).is_in(["./.", "./.", "."])
            ).height
            pct = both_valid / n * 100 if n > 0 else 0
            rows.append({"caller_1": c1.replace("DNA_", "D_").replace("RNA_", "R_"),
                         "caller_2": c2.replace("DNA_", "D_").replace("RNA_", "R_"),
                         "pct": pct})
    pdf = pl.DataFrame(rows).pipe(_maybe_collect).to_pandas()
    chart = alt.Chart(pdf).mark_rect().encode(
        x=alt.X("caller_1:N", title=None),
        y=alt.Y("caller_2:N", title=None),
        color=alt.Color("pct:Q", title="% Both Valid GT", scale=alt.Scale(scheme="blues")),
        tooltip=["caller_1", "caller_2", "pct"],
    ).properties(title="Caller GT Availability Matrix (%)")
    _save_chart(chart, "28_caller_agreement", output_dir)
    return chart


def _chart_section(fig, index: int) -> str:
    """Determine dashboard section for a chart based on its title."""
    title = ""
    try:
        title = fig.title if hasattr(fig, 'title') else ""
        title = str(title).lower()
    except Exception:
        pass
    # Heuristic: match chart title keywords to sections
    if any(kw in title for kw in ["vc", "classification", "variant type", "ti/tv", "cosmic", "filter dist"]):
        return "overview"
    if any(kw in title for kw in ["tier", "caller support", "c1", "gt concordance per"]):
        return "tier"
    if any(kw in title for kw in ["bam", "coverage", "validation", "ref_alt", "vaf v", "dp v", "dp per"]):
        return "bam"
    if any(kw in title for kw in ["per-sample", "sample", "cross"]):
        return "persample"
    return "overview"


def _extract_body_content(html: str) -> str:
    """Extract inner content between <body> and </body> from a full HTML document.

    Each altair chart's to_html() returns a full HTML document. Concatenating
    multiple full documents breaks browser rendering (only first chart shows).
    This extracts just the body content for assembly into a single valid page.
    """
    # Find <body> tag (may have attributes)
    import re
    body_start = re.search(r'<body[^>]*>', html)
    body_end = html.rfind('</body>')
    if body_start and body_end != -1:
        return html[body_start.end():body_end].strip()
    # Fallback: return as-is if body tags not found
    return html


def generate_dashboard(figs: list, output_dir: str):
    """Combine all figures into a single dashboard HTML with 4 sections.

    Sections: Overview, Tier Analysis, BAM & Validation, Per-Sample.
    Each chart gets a unique div ID to avoid vega-embed conflicts.
    """
    dashboard_path = Path(output_dir) / "dashboard.html"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    html_parts = [
        "<!DOCTYPE html><html><head><meta charset='utf-8'>",
        "<title>Seq2neo Variant Statistics Dashboard</title>",
        "<style>",
        "body{font-family:Arial,sans-serif;margin:20px;background:#f5f5f5}",
        ".chart{margin-bottom:40px;background:white;padding:20px;border-radius:8px;box-shadow:0 2px 4px rgba(0,0,0,0.1)}",
        "h1{color:#333}h2{color:#555;margin-top:30px;border-bottom:2px solid #ddd;padding-bottom:8px}",
        ".toc{background:white;padding:15px;border-radius:8px;margin-bottom:30px}",
        ".toc a{color:#1f77b4;text-decoration:none;margin-right:15px}",
        "</style>",
        "</head><body>",
        "<h1>Seq2neo Variant Statistics Dashboard</h1>",
        # Table of contents
        '<div class="toc">',
        '<a href="#overview">Overview</a>',
        '<a href="#tier">Tier Analysis</a>',
        '<a href="#bam">BAM & Validation</a>',
        '<a href="#persample">Per-Sample</a>',
        '</div>',
    ]

    # Section tracking: insert headers when section changes
    # Chart categories: overview, tier, bam, persample, validation
    section_order = [
        ("overview", "Overview"),
        ("tier", "Tier Analysis"),
        ("bam", "BAM & Validation"),
        ("persample", "Per-Sample"),
        ("validation", "Validation"),
    ]
    current_section = None

    for i, fig in enumerate(figs):
        if fig is None:
            continue
        # Determine section from chart filename
        section = _chart_section(fig, i)
        if section != current_section:
            current_section = section
            section_title = dict(section_order).get(section, section.title())
            html_parts.append(f'<h2 id="{section}">{section_title}</h2>')

        html_parts.append(f'<div class="chart" id="chart-{i}">')
        body_content = _extract_body_content(fig.to_html(output_div=f"vis-{i}", inline=True))
        html_parts.append(body_content)
        html_parts.append('</div>')

    html_parts.append("</body></html>")

    with open(dashboard_path, "w") as fh:
        fh.write("\n".join(html_parts))

    print(f"Dashboard written: {dashboard_path}")
