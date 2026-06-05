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
        color=alt.Color("VC:N", scale=alt.Scale(domain=VC_DOMAIN, range=VC_COLORS)),
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
    """Chart 11: Cross-modality & rescue analysis, side-by-side bars."""
    cols_needed = ["CROSS_MODALITY", "RESCUED", "set_number"]
    if not all(c in df.columns for c in cols_needed):
        return

    subcharts = []
    for col in ["CROSS_MODALITY", "RESCUED"]:
        pdf = df.group_by(["set_number", col]).agg(
            pl.len().alias("count")
        ).to_pandas()
        pdf["set_number"] = pdf["set_number"].astype(str)
        c = alt.Chart(pdf).mark_bar().encode(
            x=alt.X("set_number:N", title="Set"),
            y=alt.Y("count:Q", title="Count"),
            color=alt.Color(f"{col}:N"),
        ).properties(title=f"{'Cross-Modality' if col == 'CROSS_MODALITY' else 'Rescued Variants'}")
        subcharts.append(c)

    chart = alt.hconcat(*subcharts).properties(
        title="Cross-Modality and Rescue Analysis"
    )
    _save_chart(chart, "11_cross_modality", output_dir)
    return chart


def plot_ravex_breakdown(ravex_counts: dict[str, int], output_dir: str):
    """Chart 12: RaVeX filter reason breakdown, horizontal bar."""
    if not ravex_counts:
        return
    items = sorted(ravex_counts.items(), key=lambda x: x[1])
    pdf = pl.DataFrame({
        "reason": [it[0] for it in items],
        "count": [it[1] for it in items],
    }).to_pandas()
    chart = alt.Chart(pdf).mark_bar().encode(
        y=alt.Y("reason:N", title="Filter Reason", sort=None),
        x=alt.X("count:Q", title="Count"),
    ).properties(title="RaVeX Filter Reason Breakdown")
    _save_chart(chart, "12_ravex_breakdown", output_dir)
    return chart


def plot_per_sample_violin(sample_stats_df, output_dir: str):
    """Chart 8: Per-sample variant count distribution, boxplot by disease.

    Uses mark_boxplot — altair violin requires complex transform_density
    that is impractical with small sample counts per disease.
    """
    if sample_stats_df is None or (hasattr(sample_stats_df, 'is_empty') and sample_stats_df.is_empty()):
        return
    if "total_variants" not in sample_stats_df.columns:
        return
    if "disease" not in sample_stats_df.columns:
        sample_stats_df = sample_stats_df.with_columns(pl.lit("unknown").alias("disease"))
    pdf = sample_stats_df.select(["total_variants", "disease"]).to_pandas()
    chart = alt.Chart(pdf).mark_boxplot().encode(
        x=alt.X("disease:N", title="Disease", axis=alt.Axis(labelAngle=-45)),
        y=alt.Y("total_variants:Q", title="Total Variants per Sample"),
        color=alt.Color("disease:N"),
    ).properties(title="Per-Sample Variant Count Distribution by Disease")
    _save_chart(chart, "08_per_sample_violin", output_dir)
    return chart


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


def generate_dashboard(figs: list, output_dir: str):
    """Combine all figures into a single dashboard HTML using vega-embed."""
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

    for fig in figs:
        if fig is not None:
            html_parts.append('<div class="chart">')
            html_parts.append(fig.to_html())
            html_parts.append('</div>')

    html_parts.append("</body></html>")

    with open(dashboard_path, "w") as fh:
        fh.write("\n".join(html_parts))

    print(f"Dashboard written: {dashboard_path}")
