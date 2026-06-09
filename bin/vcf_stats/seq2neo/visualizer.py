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


def _to_pandas(df: pl.DataFrame, cols: list[str]) -> "pd.DataFrame":
    """Convert polars subset to pandas for altair."""
    existing = [c for c in cols if c in df.columns]
    return df.select(existing).to_pandas()


# ── Chart functions ───────────────────────────────────────────────────────


def plot_vc_distribution(df, output_dir: str, group_col: str = "set_number"):
    """Chart 1: Variant counts by VC classification, stacked bar per set."""
    if "VC" not in df.columns:
        return
    counts = (
        df.group_by([group_col, "VC"]).agg(pl.len().alias("count"))
        .sort([group_col, "VC"]).to_pandas()
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
    """Chart 2: Caller support histogram (% by N_SUPPORT_CALLERS)."""
    if "N_SUPPORT_CALLERS" not in df.columns or "set_number" not in df.columns:
        return
    counts = df.group_by(["set_number", "N_SUPPORT_CALLERS"]).agg(pl.len().alias("count"))
    total_per_set = df.group_by("set_number").agg(pl.len().alias("total"))
    counts = counts.join(total_per_set, on="set_number").with_columns(
        (pl.col("count") / pl.col("total") * 100).alias("pct")
    ).to_pandas()
    counts["set_number"] = counts["set_number"].astype(str)
    chart = alt.Chart(counts).mark_bar().encode(
        x=alt.X("N_SUPPORT_CALLERS:O", title="Number of Supporting Callers"),
        y=alt.Y("pct:Q", title="% of Variants"),
        color=alt.Color("set_number:N", title="Set"),
        column=alt.Column("set_number:N", title="Set"),
    ).properties(title="Caller Support Distribution")
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
    melted = melted.to_pandas()
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
    pdf = df.select(["DNA_VAF_mean", "RNA_VAF_mean", "VC"]).drop_nulls()
    if pdf.height > 5000:
        pdf = pdf.sample(5000)
    pdf = pdf.to_pandas()
    chart = alt.Chart(pdf).mark_circle(opacity=0.4, size=20).encode(
        x=alt.X("DNA_VAF_mean:Q", title="DNA Mean VAF"),
        y=alt.Y("RNA_VAF_mean:Q", title="RNA Mean VAF"),
        color=alt.Color("VC:N", scale=alt.Scale(domain=VC_DOMAIN, range=VC_COLORS)),
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
    pdf = pdf.to_pandas()
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

    def count_agreement(row):
        gts = [g for g in [row.get(c) for c in existing]
               if g is not None and g not in ("./.", "./.", ".")]
        if len(gts) < 2:
            return 0
        best = Counter(gts).most_common(1)[0][1]
        return best if best >= 2 else 0

    gt_subset = df.select(existing).to_dicts()
    agreements = [count_agreement(r) for r in gt_subset]
    n_total = len([a for a in agreements if a > 0])
    if n_total == 0:
        return

    counts = pl.DataFrame({
        "agreement_level": ["2", "3", "4"],
        "count": [
            sum(1 for ag in agreements if ag == 2),
            sum(1 for ag in agreements if ag == 3),
            sum(1 for ag in agreements if ag == 4),
        ],
    }).to_pandas()

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
        }).to_pandas()
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
        }).to_pandas()
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
    ).to_pandas()
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
        (pl.col("ti") / pl.col("tv")).alias("ratio")).to_pandas()
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
        ).to_pandas()
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


def plot_per_sample_distribution(sample_stats_df, output_dir: str):
    """Chart 8: Per-sample variant count, horizontal bar chart with sample IDs.

    Sorted by variant count descending. Each bar is one sample, labeled by sample_id.
    """
    if sample_stats_df is None or (hasattr(sample_stats_df, 'is_empty') and sample_stats_df.is_empty()):
        return
    if "total_variants" not in sample_stats_df.columns:
        return
    if "sample_id" not in sample_stats_df.columns:
        return

    pdf = sample_stats_df.select(["sample_id", "total_variants"]).sort(
        "total_variants", descending=True
    ).to_pandas()
    chart = alt.Chart(pdf).mark_bar().encode(
        y=alt.Y("sample_id:N", title="Sample ID", sort=None),
        x=alt.X("total_variants:Q", title="Total Variants per Sample"),
        tooltip=["sample_id", "total_variants"],
    ).properties(title="Per-Sample Variant Counts")
    _save_chart(chart, "08_per_sample_distribution", output_dir)
    return chart


# Alias for backward compatibility
plot_per_sample_violin = plot_per_sample_distribution


def plot_validation_heatmap(report, output_dir: str):
    """Chart 13: Rescue VCF validation heatmap (mismatch % per metric × sample)."""
    if report is None or (hasattr(report, 'is_empty') and report.is_empty()):
        return
    pivot = report.to_pandas().pivot(
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


def plot_vaf_violin_per_tier(df, output_dir: str):
    """Violin-style VAF distribution per caller, faceted by caller tier.

    Uses altair transform_density to create mirrored density (violin) plots.
    Samples data to 50K rows for performance.
    """
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
    pdf = melted.to_pandas()
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


def plot_dp_violin_per_tier(df, output_dir: str):
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
    pdf = melted.to_pandas()
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
    rows = []
    for row in df.select(existing_gt + ["caller_tier"]).to_dicts():
        gts = [g for g in [row.get(c) for c in existing_gt]
               if g is not None and g not in ("./.", "./.", ".")]
        if len(gts) < 2:
            continue
        best = Counter(gts).most_common(1)[0][1]
        if best >= 2:
            rows.append({"caller_tier": row["caller_tier"], "agreement": best})

    if not rows:
        return

    pdf = pl.DataFrame(rows).group_by(["caller_tier", "agreement"]).agg(
        pl.len().alias("count")
    ).to_pandas()
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
    ).to_pandas()

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
    ).to_pandas()

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
        pdf_pd = pdf.to_pandas()
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
    """Combine all figures into a single dashboard HTML using vega-embed.

    Extracts body content from each chart's to_html() (which returns full
    HTML documents) and assembles into a single valid page. Only the first
    chart's <style> block is preserved to avoid CSS conflicts.
    """
    dashboard_path = Path(output_dir) / "dashboard.html"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    html_parts = [
        "<!DOCTYPE html><html><head><meta charset='utf-8'>",
        "<title>Seq2neo Variant Statistics Dashboard</title>",
        "<script src='https://cdn.jsdelivr.net/npm/vega@5'></script>",
        "<script src='https://cdn.jsdelivr.net/npm/vega-lite@5'></script>",
        "<script src='https://cdn.jsdelivr.net/npm/vega-embed@6'></script>",
        "<style>body{font-family:Arial,sans-serif;margin:20px;background:#f5f5f5}",
        ".chart{margin-bottom:40px;background:white;padding:20px;border-radius:8px;box-shadow:0 2px 4px rgba(0,0,0,0.1)}",
        "h1{color:#333}</style>",
        "</head><body>",
        "<h1>Seq2neo Variant Statistics Dashboard</h1>",
    ]

    for i, fig in enumerate(figs):
        if fig is not None:
            html_parts.append(f'<div class="chart" id="chart-{i}">')
            body_content = _extract_body_content(fig.to_html())
            html_parts.append(body_content)
            html_parts.append('</div>')

    html_parts.append("</body></html>")

    with open(dashboard_path, "w") as fh:
        fh.write("\n".join(html_parts))

    print(f"Dashboard written: {dashboard_path}")
