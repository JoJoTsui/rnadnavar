"""Generate interactive Altair charts and static image exports.

Creates 22 chart types from seq2neo variant statistics. Each chart is
exported as interactive HTML (embeds vega-embed), PNG (via vl-convert),
and SVG (via vl-convert). No Chrome dependency.

All chart functions accept pl.LazyFrame input (from pl.scan_parquet)
and use explicit .collect() after column projection and sampling.
No function calls .height, .iter_rows(), or subscript access on a
LazyFrame without first calling .collect().
"""

from pathlib import Path

import altair as alt
import pandas as pd
import polars as pl

# Disable altair's 5000-row default limit
alt.data_transformers.disable_max_rows()


# ── Scientific publishing theme (Section 10) ─────────────────────────────
def _register_publishing_theme():
    """Register a 'publishing' Altair theme with clean, journal-ready styling."""
    def _theme():
        return {
            "config": {
                "background": "white",
                "font": "Arial",
                "axis": {
                    "labelFontSize": 11,
                    "titleFontSize": 13,
                    "titleFontWeight": "normal",
                    "gridColor": "#e0e0e0",
                    "gridOpacity": 0.5,
                    "domainColor": "#333",
                    "tickColor": "#333",
                },
                "header": {"labelFontSize": 12, "titleFontSize": 14},
                "legend": {"labelFontSize": 10, "titleFontSize": 11},
                "title": {"fontSize": 15, "fontWeight": "bold"},
                "view": {"strokeWidth": 0},
            }
        }
    alt.themes.register("publishing", _theme)


_register_publishing_theme()

# ── Color palettes ────────────────────────────────────────────────────────
CLASSIFICATION_DOMAIN = ["Somatic", "Germline", "Reference", "Artifact", "RNAedit", "NoConsensus"]
CLASSIFICATION_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]


# ═══════════════════════════════════════════════════════════════════════════
# Helpers — row counting and safe sampling on LazyFrame
# ═══════════════════════════════════════════════════════════════════════════

def _maybe_collect(df):
    """Collect a LazyFrame, pass through an eager DataFrame.

    Used to normalize inputs that may come from pl.scan_parquet() (lazy)
    or test fixtures (eager). Safe to call on either type.
    """
    try:
        if isinstance(df, pl.LazyFrame):
            return df.collect()
    except pl.exceptions.ColumnNotFoundError:
        return pl.DataFrame()
    return df


def _has_column(df, col: str) -> bool:
    """Check if a column exists, safely on both LazyFrame and DataFrame.

    Uses collect_schema().names() for LazyFrame (avoids PerformanceWarning),
    or .columns for eager DataFrame.
    """
    if isinstance(df, pl.LazyFrame):
        return col in df.collect_schema().names()
    return col in df.columns


def _count_rows(df) -> int:
    """Efficient row count on LazyFrame or DataFrame.

    Uses .select(pl.len()).collect().item() which reads no data columns
    — only parquet metadata and row-group statistics.
    """
    if isinstance(df, pl.LazyFrame):
        return df.select(pl.len()).collect().item()
    return df.height


def _sample_if_large(df, max_rows: int = 5000) -> pl.DataFrame:
    """Collect a LazyFrame, sampling if it exceeds max_rows.

    Returns an eager pl.DataFrame suitable for .to_pandas().
    Safe to call on already-eager DataFrames (pass-through with sampling).
    Collects first (LazyFrame.sample() not available in polars < 1.42),
    then samples the eager frame.
    """
    try:
        if isinstance(df, pl.LazyFrame):
            df = df.collect()
    except pl.exceptions.ColumnNotFoundError:
        return pl.DataFrame()
    if df.height > max_rows:
        return df.sample(max_rows)
    return df


def _sort_chromosomes(df: pl.DataFrame, chrom_col: str = "CHROM") -> pl.DataFrame:
    """Sort a DataFrame by natural chromosome order (chr1..22, chrX, chrY, chrM).

    Used by chromosome-grouped charts to replace the default alphanumeric sort
    (which puts chr10 before chr2). Adds a _sort_idx column, sorts, then drops it.

    Chromosomes not in the known order (e.g., alt contigs) sort after known ones,
    alphabetically.
    """
    # Build sort order: chr1..22 → 1..22, chrX→23, chrY→24, chrM/MT→25
    chrom_order = {f"chr{i}": i for i in range(1, 23)}
    chrom_order.update({"chrX": 23, "chrY": 24, "chrM": 25, "chrMT": 25})
    for i in range(1, 23):
        chrom_order[str(i)] = i
    chrom_order.update({"X": 23, "Y": 24, "M": 25, "MT": 25})

    sort_expr = pl.lit(99)  # default: unknown contigs sort last
    for chrom, idx in sorted(chrom_order.items(), key=lambda x: x[1]):
        sort_expr = pl.when(pl.col(chrom_col) == chrom).then(pl.lit(idx)).otherwise(sort_expr)

    return df.with_columns(sort_expr.alias("_sort_idx")).sort(
        ["_sort_idx", chrom_col]
    ).drop("_sort_idx")


def _chrom_sort_list(df: pl.DataFrame, col: str = "CHROM") -> list[str]:
    """Return the sorted chromosome list for passing to altair's sort= parameter.

    Call after _sort_chromosomes() to extract the ordered list of chromosome names.
    Altair's :N nominal type ignores DataFrame row order by default — the sort=
    parameter is the ONLY way to control axis order for nominal axes.
    """
    return _sort_chromosomes(df, col)[col].to_list()

# ═══════════════════════════════════════════════════════════════════════════
# Chromosome natural sort order
# ═══════════════════════════════════════════════════════════════════════════

_CHROMOSOME_ORDER = (
    [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    + [str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"]
)


# ═══════════════════════════════════════════════════════════════════════════
# Chart factories — parameterized by wise dimension
# ═══════════════════════════════════════════════════════════════════════════

def _plot_bar_wise(data, output_dir: str, group_col: str, value_col: str,
                   chart_id: str, title: str, color_col: str = None,
                   pct: bool = False):
    """Generic bar chart factory for any wise dimension.

    Args:
        data: polars DataFrame (eager, pre-aggregated).
        group_col: Column for x-axis categories.
        value_col: Column for y-axis values.
        chart_id: Unique chart filename identifier.
        title: Chart title.
        color_col: Optional column for color encoding (stacked bars).
        pct: If True, add percentage text marks.
    """
    pdf = data.to_pandas()
    if color_col:
        chart = alt.Chart(pdf).mark_bar().encode(
            x=alt.X(f"{group_col}:N", title=group_col.replace("_", " ").title()),
            y=alt.Y(f"{value_col}:Q", title=value_col.replace("_", " ").title()),
            color=alt.Color(f"{color_col}:N", title=color_col.replace("_", " ").title()),
        ).properties(title=title)
    else:
        chart = alt.Chart(pdf).mark_bar().encode(
            x=alt.X(f"{group_col}:N", title=group_col.replace("_", " ").title()),
            y=alt.Y(f"{value_col}:Q", title=value_col.replace("_", " ").title()),
        ).properties(title=title)

    if pct and group_col in pdf.columns and value_col in pdf.columns:
        text = alt.Chart(pdf).mark_text(dy=-8, fontSize=9).encode(
            x=alt.X(f"{group_col}:N"),
            y=alt.Y(f"{value_col}:Q"),
            text=alt.Text(f"{value_col}:Q", format=".1f"),
        )
        chart = chart + text

    _save_chart(chart, chart_id, output_dir)
    return chart


def _plot_box_violin_wise(data, output_dir: str, group_col: str, value_col: str,
                          chart_id: str, title: str, color_col: str = None):
    """Box + violin overlay chart for VAF/DP distributions.

    Args:
        data: polars DataFrame (eager, sampled to 50K max).
        group_col: Column for x-axis categories.
        value_col: Column for y-axis values.
    """
    pdf = data.to_pandas()
    violin = alt.Chart(pdf).transform_density(
        value_col, groupby=[group_col] if group_col in pdf.columns else None,
    ).mark_area(opacity=0.3).encode(
        x=alt.X(f"{group_col}:N", title=group_col.replace("_", " ").title()),
        y=alt.Y(f"{value_col}:Q", title=value_col.replace("_", " ").title()),
        color=alt.Color(f"{color_col or group_col}:N") if color_col else alt.value("#1f77b4"),
    )
    box = alt.Chart(pdf).mark_boxplot(size=30).encode(
        x=alt.X(f"{group_col}:N"),
        y=alt.Y(f"{value_col}:Q"),
        color=alt.Color(f"{color_col or group_col}:N") if color_col else alt.value("#1f77b4"),
    )
    chart = (violin + box).properties(title=title)
    _save_chart(chart, chart_id, output_dir)
    return chart


def _plot_scatter_wise(data, output_dir: str, x_col: str, y_col: str,
                       chart_id: str, title: str, color_col: str = None,
                       opacity: float = 0.4, size: int = 20):
    """Scatter plot factory with configurable coloring.

    Args:
        data: polars DataFrame (eager, sampled to 5K max).
    """
    pdf = data.to_pandas()
    encodings = {
        "x": alt.X(f"{x_col}:Q", title=x_col.replace("_", " ").title()),
        "y": alt.Y(f"{y_col}:Q", title=y_col.replace("_", " ").title()),
    }
    if color_col and color_col in pdf.columns:
        encodings["color"] = alt.Color(f"{color_col}:N", title=color_col)
    chart = alt.Chart(pdf).mark_circle(opacity=opacity, size=size).encode(**encodings).properties(title=title)
    _save_chart(chart, chart_id, output_dir)
    return chart


def _plot_pie_wise(data, output_dir: str, category_col: str, count_col: str,
                   chart_id: str, title: str, color_domain: list = None,
                   color_range: list = None):
    """Pie/donut chart factory.

    Args:
        data: polars DataFrame with category and count columns (eager).
    """
    pdf = data.to_pandas()
    color_kw = {}
    if color_domain and color_range:
        color_kw = {"scale": alt.Scale(domain=color_domain, range=color_range)}
    chart = alt.Chart(pdf).mark_arc(innerRadius=40).encode(
        theta=alt.Theta(f"{count_col}:Q"),
        color=alt.Color(f"{category_col}:N", **color_kw),
    ).properties(title=title)
    _save_chart(chart, chart_id, output_dir)
    return chart


def _save_chart(chart: alt.Chart, name: str, output_dir: str):
    """Save a chart as HTML, PNG, and SVG.

    output_dir is the directory where plots/ is created. If output_dir
    already contains 'plots/' (e.g., 'stats/full/plots/set/'), charts
    are saved directly there.
    """
    plots_dir = Path(output_dir)
    if "plots" not in plots_dir.parts:
        plots_dir = plots_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # HTML always works (no external dependencies)
    html_path = plots_dir / f"{name}.html"
    chart.save(str(html_path))

    # PNG/SVG require vl-convert — best-effort
    png_path = plots_dir / f"{name}.png"
    try:
        chart.save(str(png_path), format="png", scale_factor=3)
    except (ImportError, ModuleNotFoundError):
        pass  # vl-convert not installed, skip PNG

    svg_path = plots_dir / f"{name}.svg"
    try:
        chart.save(str(svg_path), format="svg")
    except (ImportError, ModuleNotFoundError):
        pass  # vl-convert not installed, skip SVG


# ═══════════════════════════════════════════════════════════════════════════
# Aggregate charts — group_by → small result, no sampling needed.
# Memory: < 1 MB per chart (result is tens to hundreds of rows).
# ═══════════════════════════════════════════════════════════════════════════

def plot_vc_distribution(df, output_dir: str, group_col: str = "set_number"):
    """Chart 1: Variant counts by VC classification, stacked bar per set."""
    if not _has_column(df, "FILTER"):
        return
    counts = _maybe_collect(
        df.group_by([group_col, "FILTER"]).agg(pl.len().alias("count"))
        .sort([group_col, "FILTER"])
    )
    if group_col == "CHROM":
        counts = _sort_chromosomes(counts, group_col)
    chrom_order = counts[group_col].to_list() if group_col == "CHROM" else None
    counts = counts.to_pandas()
    counts[group_col] = counts[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()
    chart = alt.Chart(counts).mark_bar().encode(
        x=alt.X(f"{group_col}:N", title=group_title, sort=chrom_order),
        y=alt.Y("count:Q", title="Number of Variants"),
        color=alt.Color(
            "FILTER:N", title="Variant Classification",
            scale=alt.Scale(domain=CLASSIFICATION_DOMAIN, range=CLASSIFICATION_COLORS),
            legend=alt.Legend(orient="right", title="Variant Classification",
                              labelFontSize=11, titleFontSize=12),
        ),
    ).properties(title=f"Variant Classification Distribution by {group_title}")
    _save_chart(chart, "01_vc_distribution", output_dir)
    return chart


def plot_caller_overlap(df, output_dir: str, group_col: str = "set_number"):
    """Chart 2: Variant tier distribution (% by CxDy final tier)."""
    if "final_tier" not in df.columns or group_col not in df.columns:
        return
    counts = df.group_by([group_col, "final_tier"]).agg(pl.len().alias("count"))
    total_per_group = df.group_by(group_col).agg(pl.len().alias("total"))
    counts = counts.join(total_per_group, on=group_col).with_columns(
        (pl.col("count") / pl.col("total") * 100).alias("pct")
    )
    counts = _maybe_collect(counts).to_pandas()
    counts[group_col] = counts[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()
    chart = alt.Chart(counts).mark_bar().encode(
        x=alt.X("final_tier:N", title="Variant Tier (CxDy)"),
        y=alt.Y("pct:Q", title="% of Variants"),
        color=alt.Color(f"{group_col}:N", title=group_title),
        column=alt.Column(f"{group_col}:N", title=group_title),
    ).properties(title=f"Variant Tier Distribution by {group_title}")
    _save_chart(chart, "02_caller_overlap", output_dir)
    return chart


def plot_variant_type_distribution(df, output_dir: str, group_col: str = "set_number"):
    """Chart 9: Variant type distribution (SNV/INS/DEL/MNV), stacked bar per group."""
    if "variant_type" not in df.columns or group_col not in df.columns:
        return
    total_per_group = df.group_by(group_col).agg(pl.len().alias("total"))
    pdf = df.group_by([group_col, "variant_type"]).agg(pl.len().alias("count"))
    pdf = pdf.join(total_per_group, on=group_col).with_columns(
        (pl.col("count") / pl.col("total") * 100).round(1).alias("pct")
    )
    pdf = _maybe_collect(pdf)
    if group_col == "CHROM":
        pdf = _sort_chromosomes(pdf, group_col)
    chrom_order = pdf[group_col].to_list() if group_col == "CHROM" else None
    pdf = pdf.to_pandas()
    pdf[group_col] = pdf[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()
    bars = alt.Chart(pdf).mark_bar().encode(
        x=alt.X(f"{group_col}:N", title=group_title, sort=chrom_order),
        y=alt.Y("count:Q", title="Number of Variants"),
        color=alt.Color("variant_type:N", title="Variant Type"),
    )
    pct_text = alt.Chart(pdf).mark_text(dy=-8, fontSize=9).encode(
        x=alt.X(f"{group_col}:N"), y=alt.Y("count:Q"),
        text=alt.Text("pct:Q", format=".1f"), color=alt.Color("variant_type:N"),
    )
    chart = (bars + pct_text).properties(title=f"Variant Type Distribution by {group_title}")
    _save_chart(chart, "09_variant_type_distribution", output_dir)
    return chart


def plot_ti_tv_ratio(df, output_dir: str, group_col: str = "set_number"):
    """Chart 10: Ti/Tv ratio bar chart per group."""
    if "ti_tv" not in df.columns or group_col not in df.columns:
        return
    ti = df.filter(pl.col("ti_tv") == True).group_by(group_col).agg(pl.len().alias("ti"))
    tv = df.filter(pl.col("ti_tv") == False).group_by(group_col).agg(pl.len().alias("tv"))
    ratio = ti.join(tv, on=group_col).with_columns(
        (pl.col("ti") / pl.col("tv")).alias("ratio"))
    ratio = _maybe_collect(ratio)
    if group_col == "CHROM":
        ratio = _sort_chromosomes(ratio, group_col)
    chrom_order = ratio[group_col].to_list() if group_col == "CHROM" else None
    ratio = ratio.to_pandas()
    ratio[group_col] = ratio[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()
    bars = alt.Chart(ratio).mark_bar().encode(
        x=alt.X(f"{group_col}:N", title=group_title, sort=chrom_order),
        y=alt.Y("ratio:Q", title="Ti/Tv Ratio"),
    )
    text = alt.Chart(ratio).mark_text(dy=-8).encode(
        x=alt.X(f"{group_col}:N"), y=alt.Y("ratio:Q"),
        text=alt.Text("ratio:Q", format=".2f"),
    )
    chart = (bars + text).properties(title=f"Ti/Tv Ratio by {group_title}")
    _save_chart(chart, "10_ti_tv_ratio", output_dir)
    return chart


def plot_cross_modality(df, output_dir: str, group_col: str = "set_number"):
    """Chart 11: Cross-modality & rescue analysis with percentages."""
    cols_needed = ["CROSS_MODALITY", "RESCUED"]
    if not all(c in df.columns for c in cols_needed) or group_col not in df.columns:
        return

    total_per_group = df.group_by(group_col).agg(pl.len().alias("total"))

    subcharts = []
    for col in ["CROSS_MODALITY", "RESCUED"]:
        pdf = df.group_by([group_col, col]).agg(
            pl.len().alias("count")
        ).join(total_per_group, on=group_col).with_columns(
            (pl.col("count") / pl.col("total") * 100).round(1).alias("pct")
        )
        pdf = _maybe_collect(pdf)
        if group_col == "CHROM":
            pdf = _sort_chromosomes(pdf, group_col)
        chrom_order = pdf[group_col].to_list() if group_col == "CHROM" else None
        pdf = pdf.to_pandas()
        pdf[group_col] = pdf[group_col].astype(str)
        group_title = group_col.replace("_", " ").title()

        bars = alt.Chart(pdf).mark_bar().encode(
            x=alt.X(f"{group_col}:N", title=group_title, sort=chrom_order),
            y=alt.Y("count:Q", title="Count"),
            color=alt.Color(f"{col}:N"),
        )
        text = alt.Chart(pdf).mark_text(dy=-8, fontSize=10).encode(
            x=alt.X(f"{group_col}:N"), y=alt.Y("count:Q"),
            text=alt.Text("pct:N", format=".1f"),
            color=alt.Color(f"{col}:N"),
        ).transform_filter(alt.datum["count"] > 0)
        c = (bars + text).properties(
            title=f"{'Cross-Modality' if col == 'CROSS_MODALITY' else 'Rescued Variants'}"
        )
        subcharts.append(c)

    chart = alt.hconcat(*subcharts).properties(title="Cross-Modality and Rescue Analysis (%)")
    _save_chart(chart, "11_cross_modality", output_dir)
    return chart


def plot_filter_distribution(df, output_dir: str, group_col: str = "set_number"):
    """FILTER value distribution per group (bar chart with % marks)."""
    if "FILTER" not in df.columns or group_col not in df.columns:
        return
    total_per_group = df.group_by(group_col).agg(pl.len().alias("total"))
    counts = df.group_by([group_col, "FILTER"]).agg(pl.len().alias("count"))
    counts = counts.join(total_per_group, on=group_col).with_columns(
        (pl.col("count") / pl.col("total") * 100).round(1).alias("pct")
    )
    counts = _maybe_collect(counts)
    if group_col == "CHROM":
        counts = _sort_chromosomes(counts, group_col)
    chrom_order = counts[group_col].to_list() if group_col == "CHROM" else None
    counts = counts.to_pandas()
    counts[group_col] = counts[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()
    bars = alt.Chart(counts).mark_bar().encode(
        x=alt.X(f"{group_col}:N", title=group_title, sort=chrom_order),
        y=alt.Y("count:Q", title="Number of Variants"),
        color=alt.Color("FILTER:N", title="FILTER"),
    )
    pct_text = alt.Chart(counts).mark_text(dy=-8, fontSize=9).encode(
        x=alt.X(f"{group_col}:N"), y=alt.Y("count:Q"),
        text=alt.Text("pct:Q", format=".1f"), color=alt.Color("FILTER:N"),
    )
    chart = (bars + pct_text).properties(title=f"FILTER Distribution by {group_title}")
    _save_chart(chart, "24_filter_distribution", output_dir)
    return chart


def plot_chromosome_density(df, output_dir: str):
    """Variant count per chromosome with natural sort order (chr1..chrX, chrY, chrM)."""
    if "CHROM" not in df.columns:
        return
    counts = df.group_by("CHROM").agg(pl.len().alias("count"))
    counts = _maybe_collect(counts)
    counts = _sort_chromosomes(counts, "CHROM")

    chart = alt.Chart(counts.to_pandas()).mark_bar().encode(
        x=alt.X("CHROM:N", title="Chromosome", sort=counts["CHROM"].to_list()),
        y=alt.Y("count:Q", title="Number of Variants"),
        tooltip=["CHROM", "count"],
    ).properties(title="Variant Density per Chromosome")
    _save_chart(chart, "30_chromosome_density", output_dir)
    return chart


def plot_redi_evidence(df, output_dir: str, group_col: str = "set_number"):
    """REDIportal RNA editing evidence distribution per group."""
    if "REDI_EVIDENCE" not in df.columns or group_col not in df.columns:
        return
    counts = df.group_by([group_col, "REDI_EVIDENCE"]).agg(pl.len().alias("count"))
    counts = _maybe_collect(counts)
    if group_col == "CHROM":
        counts = _sort_chromosomes(counts, group_col)
    chrom_order = counts[group_col].to_list() if group_col == "CHROM" else None
    counts = counts.sort([group_col, "count"], descending=[False, True]).to_pandas()
    counts[group_col] = counts[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()
    chart = alt.Chart(counts).mark_bar().encode(
        x=alt.X(f"{group_col}:N", title=group_title, sort=chrom_order),
        y=alt.Y("count:Q", title="Number of Variants"),
        color=alt.Color("REDI_EVIDENCE:N", title="REDIportal Evidence Level"),
    ).properties(title=f"REDIportal RNA Editing Evidence by {group_title}")
    _save_chart(chart, "29_redi_evidence", output_dir)
    return chart


def plot_tiered_caller_overlap(df, output_dir: str, facet_col: str = None):
    """N_SUPPORT_CALLERS histogram faceted by caller tier, with % marks and log scale."""
    if "N_SUPPORT_CALLERS" not in df.columns or "caller_tier" not in df.columns:
        return
    group_keys = ["caller_tier", "N_SUPPORT_CALLERS"]
    total_keys = ["caller_tier"]
    if facet_col and facet_col in df.columns:
        group_keys.append(facet_col)
        total_keys.append(facet_col)
    counts = df.group_by(group_keys).agg(pl.len().alias("count"))
    total_per_tier = df.group_by(total_keys).agg(pl.len().alias("total"))
    join_on = total_keys
    pdf = counts.join(total_per_tier, on=join_on).with_columns(
        (pl.col("count") / pl.col("total") * 100).round(1).alias("pct")
    )
    pdf = _maybe_collect(pdf).to_pandas()
    col_enc = alt.Column("caller_tier:N", title="Caller Tier")
    if facet_col and facet_col in pdf.columns:
        pdf[facet_col] = pdf[facet_col].astype(str)
        col_enc = alt.Column(f"{facet_col}:N", title=facet_col.replace("_", " ").title())
    bars = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("N_SUPPORT_CALLERS:O", title="Number of Supporting Callers"),
        y=alt.Y("pct:Q", title="% of Variants", scale=alt.Scale(type="log")),
        color=alt.Color("caller_tier:N", title="Caller Tier"),
        column=col_enc,
    )
    chart = bars.properties(title="Caller Support Distribution per Caller Tier")
    _save_chart(chart, "18_caller_overlap_per_tier", output_dir)
    return chart


def plot_tiered_variant_types(df, output_dir: str, facet_col: str = None):
    """Variant type distribution (SNV/INS/DEL/MNV) faceted by caller tier, with % marks."""
    if "variant_type" not in df.columns or "caller_tier" not in df.columns:
        return
    group_keys = ["caller_tier", "variant_type"]
    total_keys = ["caller_tier"]
    if facet_col and facet_col in df.columns:
        group_keys.append(facet_col)
        total_keys.append(facet_col)
    counts = df.group_by(group_keys).agg(pl.len().alias("count"))
    total_per_tier = df.group_by(total_keys).agg(pl.len().alias("total"))
    pdf = counts.join(total_per_tier, on=total_keys).with_columns(
        (pl.col("count") / pl.col("total") * 100).round(1).alias("pct")
    )
    pdf = _maybe_collect(pdf).to_pandas()
    col_enc = alt.Column("caller_tier:N", title="Caller Tier")
    if facet_col and facet_col in pdf.columns:
        pdf[facet_col] = pdf[facet_col].astype(str)
        col_enc = alt.Column(f"{facet_col}:N", title=facet_col.replace("_", " ").title())
    bars = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("variant_type:N", title="Variant Type"),
        y=alt.Y("count:Q", title="Number of Variants"),
        color=alt.Color("variant_type:N", title="Variant Type"),
        column=col_enc,
    )
    chart = bars.properties(title="Variant Type Distribution per Caller Tier")
    _save_chart(chart, "19_variant_types_per_tier", output_dir)
    return chart


def plot_tier_quality_distribution(df, output_dir: str, facet_col: str = None):
    """Tier quality score histogram (sampled — raw data too large for altair PNG)."""
    if "tier_quality" not in df.columns:
        return
    select_cols = ["tier_quality"]
    if facet_col and facet_col in df.columns:
        select_cols.append(facet_col)
    pdf = df.select(select_cols).drop_nulls()
    pdf = _sample_if_large(pdf, max_rows=50000).to_pandas()
    enc = {"x": alt.X("tier_quality:Q", bin=alt.Bin(maxbins=20), title="Tier Quality Score"),
           "y": alt.Y("count()", title="Number of Variants")}
    if facet_col and facet_col in pdf.columns:
        pdf[facet_col] = pdf[facet_col].astype(str)
        enc["color"] = alt.Color(f"{facet_col}:N")
        enc["column"] = alt.Column(f"{facet_col}:N")
    chart = alt.Chart(pdf).mark_bar().encode(**enc).properties(title="Tier Quality Score Distribution (sampled)")
    _save_chart(chart, "27_tier_quality", output_dir)
    return chart


# ═══════════════════════════════════════════════════════════════════════════
# Sampled charts — use _sample_if_large() to limit memory.
# Memory: < 50K rows per chart (sampled).
# ═══════════════════════════════════════════════════════════════════════════

def plot_vaf_distribution(df, output_dir: str, color_col: str = None):
    """Chart 3: Per-caller VAF distribution, box+violin overlay (sampled to 50K)."""
    caller_vaf_cols = [
        f"{c}_VAF" for c in [
            "DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic",
            "RNA_deepsomatic", "DNA_strelka", "RNA_strelka",
        ]
    ]
    existing = [c for c in caller_vaf_cols if c in df.columns]
    if not existing:
        return
    select_cols = existing[:]
    if color_col and color_col in df.columns:
        select_cols.append(color_col)
    melted = df.select(select_cols).unpivot(
        index=[color_col] if color_col and color_col in df.columns else [],
        variable_name="caller", value_name="VAF"
    ).drop_nulls()
    pdf = _sample_if_large(melted, max_rows=50000).to_pandas()
    pdf["caller"] = pdf["caller"].str.replace("_VAF", "")

    # Clamp VAF to [0, 1] for visualization (Strelka VAF uses tier-1 depth
    # denominator and can exceed 1.0 — correct by design, clamped for display).
    vaf_clamp = pdf["VAF"].clip(0.0, 1.0)
    pdf["VAF_display"] = vaf_clamp

    y_scale = alt.Y("VAF_display:Q", title="Variant Allele Frequency (capped at 1.0)",
                    scale=alt.Scale(domain=[0, 1]))
    enc = {"x": alt.X("caller:N", title="Caller", axis=alt.Axis(labelAngle=-45)),
           "y": y_scale}
    if color_col and color_col in pdf.columns:
        enc["color"] = alt.Color(f"{color_col}:N")
        enc["column"] = alt.Column(f"{color_col}:N")

    # Horizontal reference lines at key VAF thresholds
    ref_line_data = pd.DataFrame({"y": [0.005, 0.01, 0.05, 0.10]})
    ref_rules = alt.Chart(ref_line_data).mark_rule(
        strokeDash=[2, 2], opacity=0.4, strokeWidth=1
    ).encode(y=alt.Y("y:Q"))

    is_faceted = "column" in enc
    if not is_faceted:
        enc["color"] = alt.Color("caller:N", scale=alt.Scale(scheme="category10"))
        violin = alt.Chart(pdf).transform_density(
            "VAF_display", groupby=["caller"]
        ).mark_area(opacity=0.3).encode(
            x=alt.X("caller:N", title="Caller", axis=alt.Axis(labelAngle=-45)),
            y=y_scale,
            color=alt.Color("caller:N", scale=alt.Scale(scheme="category10")))
        box = alt.Chart(pdf).mark_boxplot(size=30).encode(**enc)
        chart = (violin + box + ref_rules).properties(title="VAF Distribution per Caller")
    else:
        chart = alt.Chart(pdf).mark_boxplot(size=30).encode(**enc).properties(
            title="VAF Distribution per Caller")
    _save_chart(chart, "03_vaf_distribution", output_dir)
    return chart


def plot_dna_vs_rna_vaf(df, output_dir: str, color_col: str = None):
    """Chart 4: DNA vs RNA mean VAF scatter (sampled to 5K)."""
    if "DNA_VAF_mean" not in df.columns or "RNA_VAF_mean" not in df.columns:
        return
    cols = ["DNA_VAF_mean", "RNA_VAF_mean"]
    if color_col and color_col in df.columns:
        cols.append(color_col)
    elif "FILTER" in df.columns:
        cols.append("FILTER")  # fallback: color by classification
    pdf = df.select(cols).drop_nulls(subset=["DNA_VAF_mean", "RNA_VAF_mean"])
    pdf = _sample_if_large(pdf, max_rows=5000).to_pandas()
    # Determine color encoding
    actual_color = None
    if color_col and color_col in pdf.columns:
        actual_color = alt.Color(f"{color_col}:N", title=color_col.replace("_", " ").title())
    elif "FILTER" in pdf.columns:
        actual_color = alt.Color("FILTER:N", scale=alt.Scale(domain=CLASSIFICATION_DOMAIN, range=CLASSIFICATION_COLORS))
    else:
        actual_color = alt.value("#1f77b4")
    chart = alt.Chart(pdf).mark_circle(opacity=0.4, size=20).encode(
        x=alt.X("DNA_VAF_mean:Q", title="DNA Mean VAF"),
        y=alt.Y("RNA_VAF_mean:Q", title="RNA Mean VAF"),
        color=actual_color,
    ).properties(title="DNA vs RNA Mean VAF")
    _save_chart(chart, "04_dna_vs_rna_vaf", output_dir)
    return chart


def plot_dna_vs_rna_dp(df, output_dir: str, group_col: str = "set_number"):
    """Chart 5: DNA vs RNA mean DP scatter (sampled to 5K)."""
    if "DNA_DP_mean" not in df.columns or "RNA_DP_mean" not in df.columns:
        return
    cols = ["DNA_DP_mean", "RNA_DP_mean"]
    if group_col in df.columns:
        cols.append(group_col)
    pdf = df.select(cols).drop_nulls()
    pdf = _sample_if_large(pdf, max_rows=5000).to_pandas()
    if group_col in pdf.columns:
        pdf[group_col] = pdf[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()
    enc = {"x": alt.X("DNA_DP_mean:Q", title="DNA Mean Depth (capped at 2000)", scale=alt.Scale(domain=[0, 2000])),
           "y": alt.Y("RNA_DP_mean:Q", title="RNA Mean Depth (capped at 2000)", scale=alt.Scale(domain=[0, 2000]))}
    if group_col in pdf.columns:
        enc["color"] = alt.Color(f"{group_col}:N", title=group_title)
    chart = alt.Chart(pdf).mark_circle(opacity=0.4, size=20).encode(**enc).properties(
        title=f"DNA vs RNA Mean Depth (capped at 2000)")
    _save_chart(chart, "05_dna_vs_rna_dp", output_dir)
    return chart


def plot_vaf_boxplot_per_tier(df, output_dir: str, facet_col: str = None):
    """VAF distribution per caller, boxplot faceted by caller tier (sampled to 50K)."""
    caller_vaf_cols = [
        f"{c}_VAF" for c in [
            "DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic",
            "RNA_deepsomatic", "DNA_strelka", "RNA_strelka",
        ]
    ]
    existing = [c for c in caller_vaf_cols if c in df.columns]
    if not existing or "caller_tier" not in df.columns:
        return
    index_cols = ["caller_tier"]
    if facet_col and facet_col in df.columns:
        index_cols.append(facet_col)
    melted = df.select(existing + index_cols).unpivot(
        index=index_cols, variable_name="caller", value_name="VAF"
    ).drop_nulls()
    pdf = _sample_if_large(melted, max_rows=50000).to_pandas()
    pdf["caller"] = pdf["caller"].str.replace("_VAF", "")
    col_enc = alt.Column("caller_tier:N", title="Caller Tier")
    if facet_col and facet_col in pdf.columns:
        pdf[facet_col] = pdf[facet_col].astype(str)
        col_enc = alt.Column(f"{facet_col}:N", title=facet_col.replace("_", " ").title())
    pdf["VAF_display"] = pdf["VAF"].clip(0.0, 1.0)
    y_scale = alt.Y("VAF_display:Q", title="Variant Allele Frequency (capped at 1.0)",
                    scale=alt.Scale(domain=[0, 1]))
    base_enc = {"x": alt.X("caller:N", title="Caller", axis=alt.Axis(labelAngle=-45)),
                "color": alt.Color("caller:N", scale=alt.Scale(scheme="category10")),
                "column": col_enc,
                "y": y_scale}
    chart = alt.Chart(pdf).mark_boxplot(size=30).encode(**base_enc).properties(
        title="VAF Distribution per Caller × Caller Tier")
    _save_chart(chart, "14_vaf_violin_per_tier", output_dir)
    return chart


def plot_dp_boxplot_per_tier(df, output_dir: str, facet_col: str = None):
    """DP distribution per caller, box+violin faceted by caller tier (sampled to 50K)."""
    caller_dp_cols = [
        f"{c}_DP" for c in [
            "DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic",
            "RNA_deepsomatic", "DNA_strelka", "RNA_strelka",
        ]
    ]
    existing = [c for c in caller_dp_cols if c in df.columns]
    if not existing or "caller_tier" not in df.columns:
        return
    index_cols = ["caller_tier"]
    if facet_col and facet_col in df.columns:
        index_cols.append(facet_col)
    melted = df.select(existing + index_cols).unpivot(
        index=index_cols, variable_name="caller", value_name="DP"
    ).drop_nulls()
    pdf = _sample_if_large(melted, max_rows=50000).to_pandas()
    pdf["caller"] = pdf["caller"].str.replace("_DP", "")
    col_enc = alt.Column("caller_tier:N", title="Caller Tier")
    if facet_col and facet_col in pdf.columns:
        pdf[facet_col] = pdf[facet_col].astype(str)
        col_enc = alt.Column(f"{facet_col}:N", title=facet_col.replace("_", " ").title())
    base_enc = {"x": alt.X("caller:N", title="Caller", axis=alt.Axis(labelAngle=-45)),
                "color": alt.Color("caller:N"), "column": col_enc,
                "y": alt.Y("DP:Q", title="Read Depth (capped at 2000)", scale=alt.Scale(domain=[0, 2000]))}
    chart = alt.Chart(pdf).mark_boxplot(size=30).encode(**base_enc).properties(
        title="DP Distribution per Caller × Caller Tier (capped at 2000)")
    _save_chart(chart, "15_dp_per_tier", output_dir)
    return chart


def plot_ref_alt_dp_scatter(df, output_dir: str, group_col: str = "set_number"):
    """DNA vs RNA mean REF_DP and ALT_DP scatter plots (sampled to 5K each)."""
    needed = ["DNA_REF_DP_mean", "RNA_REF_DP_mean",
              "DNA_ALT_DP_mean", "RNA_ALT_DP_mean"]
    if not all(c in df.columns for c in needed):
        return
    subcharts = []
    for label, x_col, y_col in [
        ("REF_DP", "DNA_REF_DP_mean", "RNA_REF_DP_mean"),
        ("ALT_DP", "DNA_ALT_DP_mean", "RNA_ALT_DP_mean"),
    ]:
        cols = [x_col, y_col]
        if group_col in df.columns:
            cols.append(group_col)
        pdf = df.select(cols).drop_nulls()
        pdf = _sample_if_large(pdf, max_rows=5000).to_pandas()
        group_title = group_col.replace("_", " ").title()
        enc = {"x": alt.X(f"{x_col}:Q", title=f"DNA Mean {label} (capped at 2000)", scale=alt.Scale(domain=[0, 2000])),
               "y": alt.Y(f"{y_col}:Q", title=f"RNA Mean {label} (capped at 2000)", scale=alt.Scale(domain=[0, 2000]))}
        if group_col in pdf.columns:
            pdf[group_col] = pdf[group_col].astype(str)
            enc["color"] = alt.Color(f"{group_col}:N", title=group_title)
        c = alt.Chart(pdf).mark_circle(opacity=0.4, size=20).encode(**enc).properties(
            title=f"DNA vs RNA Mean {label} (capped at 2000)")
        subcharts.append(c)
    chart = alt.hconcat(*subcharts).properties(title="DNA vs RNA REF_DP and ALT_DP (capped at 2000)")
    _save_chart(chart, "17_ref_alt_dp_scatter", output_dir)
    return chart


def plot_per_tier_vaf_boxplot(df, output_dir: str):
    """Per-tier DNA VAF boxplot across all variants (sampled to 50K)."""
    if "DNA_VAF_mean" not in df.columns or "final_tier" not in df.columns:
        return
    pdf = df.select(["final_tier", "DNA_VAF_mean", "RNA_VAF_mean"]).drop_nulls(subset=["DNA_VAF_mean"])
    pdf = _sample_if_large(pdf, max_rows=50000).to_pandas()
    pdf["VAF_display"] = pdf["DNA_VAF_mean"].clip(0.0, 1.0)
    chart = alt.Chart(pdf).mark_boxplot().encode(
        x=alt.X("final_tier:N", title="Tier"),
        y=alt.Y("VAF_display:Q", title="DNA Mean VAF (capped at 1.0)",
                scale=alt.Scale(domain=[0, 1])),
        color=alt.Color("final_tier:N", scale=alt.Scale(scheme="category10")),
    ).properties(title="DNA VAF Distribution per Tier")
    _save_chart(chart, "25_per_tier_vaf", output_dir)
    return chart


def plot_dna_vs_rna_per_caller(df, output_dir: str, color_col: str = None):
    """DNA vs RNA per-caller VAF scatter at shared positions (sampled to 5K per pair)."""
    pairs = [("DNA_mutect2", "RNA_mutect2"), ("DNA_deepsomatic", "RNA_deepsomatic"),
             ("DNA_strelka", "RNA_strelka")]
    subcharts = []
    for dna_caller, rna_caller in pairs:
        dna_vaf = f"{dna_caller}_VAF"
        rna_vaf = f"{rna_caller}_VAF"
        if dna_vaf not in df.columns or rna_vaf not in df.columns:
            continue
        cols = [dna_vaf, rna_vaf]
        active_color = None
        if color_col and color_col in df.columns:
            cols.append(color_col)
            active_color = color_col
        elif "FILTER" in df.columns:
            cols.append("FILTER")
            active_color = "FILTER"
        pdf = df.select(cols).drop_nulls(subset=[dna_vaf, rna_vaf])
        pdf = _sample_if_large(pdf, max_rows=5000).to_pandas()
        caller_label = dna_caller.replace("DNA_", "")
        if active_color and active_color in pdf.columns:
            if active_color == "FILTER":
                color_enc = alt.Color("FILTER:N", scale=alt.Scale(domain=CLASSIFICATION_DOMAIN, range=CLASSIFICATION_COLORS))
            else:
                color_enc = alt.Color(f"{active_color}:N")
        else:
            color_enc = alt.value("#1f77b4")
        c = alt.Chart(pdf).mark_circle(opacity=0.4, size=20).encode(
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


# ═══════════════════════════════════════════════════════════════════════════
# Count-based charts — use _count_rows(), no data columns loaded.
# Memory: < 1 KB per chart (just integer counts).
# ═══════════════════════════════════════════════════════════════════════════

def plot_cosmic_gnomad_annotation(df, output_dir: str, group_col: str = "set_number"):
    """Chart 7: COSMIC/gnomAD annotation coverage as side-by-side pies, per group."""
    has_cosmic_col = "COSMIC_ID" in df.columns
    has_gnomad_col = "GNOMAD_AF" in df.columns
    if not has_cosmic_col and not has_gnomad_col:
        return

    has_groups = group_col in df.columns

    if has_groups:
        # Per-group annotation bars (single group_by scan, not per-group _count_rows)
        agg_exprs = [pl.len().alias("n_total")]
        if has_cosmic_col:
            agg_exprs.append(pl.col("COSMIC_ID").is_not_null().sum().alias("n_cosmic"))
        if has_gnomad_col:
            agg_exprs.append(pl.col("GNOMAD_AF").is_not_null().sum().alias("n_gnomad"))
        counts_df = _maybe_collect(df.group_by(group_col).agg(agg_exprs))
        if counts_df.is_empty():
            return
        rows = []
        for row in counts_df.iter_rows():
            g = row[0]
            g_n = row[1]  # n_total
            r = {"group": str(g), "n_total": g_n}
            idx = 2
            if has_cosmic_col:
                r["cosmic_pct"] = round(row[idx] / g_n * 100, 2) if g_n > 0 else 0
                idx += 1
            if has_gnomad_col:
                r["gnomad_pct"] = round(row[idx] / g_n * 100, 2) if g_n > 0 else 0
                idx += 1
            rows.append(r)
        if not rows:
            return
        pdf_pl = pl.DataFrame(rows)
        if group_col == "CHROM":
            pdf_pl = _sort_chromosomes(pdf_pl, "group")
        chrom_order = pdf_pl["group"].to_list() if group_col == "CHROM" else None
        pdf = pdf_pl.to_pandas()
        group_title = group_col.replace("_", " ").title()

        subcharts = []
        for db, col, color in [("COSMIC", "cosmic_pct", "#1f77b4"), ("gnomAD", "gnomad_pct", "#ff7f0e")]:
            if col not in pdf.columns:
                continue
            c = alt.Chart(pdf).mark_bar().encode(
                x=alt.X("group:N", title=group_title, sort=chrom_order),
                y=alt.Y(f"{col}:Q", title=f"% with {db} Annotation"),
                color=alt.value(color),
            ).properties(title=f"{db} Annotation by {group_title}")
            subcharts.append(c)
        chart = alt.hconcat(*subcharts) if len(subcharts) == 2 else subcharts[0]
        chart = chart.properties(title=f"COSMIC and gnomAD Annotation by {group_title}")
    else:
        # Global pie charts
        n = _count_rows(df)
        has_cosmic = _count_rows(df.filter(pl.col("COSMIC_ID").is_not_null())) if has_cosmic_col else -1
        has_gnomad = _count_rows(df.filter(pl.col("GNOMAD_AF").is_not_null())) if has_gnomad_col else -1
        if has_cosmic < 0 and has_gnomad < 0:
            return
        charts = []
        if has_cosmic >= 0:
            cd = pl.DataFrame({
                "category": ["In COSMIC", "Not in COSMIC"],
                "count": [has_cosmic, n - has_cosmic],
            }).to_pandas()
            charts.append(alt.Chart(cd).mark_arc(innerRadius=40).encode(
                theta=alt.Theta("count:Q"),
                color=alt.Color("category:N", scale=alt.Scale(
                    domain=["In COSMIC", "Not in COSMIC"], range=["#1f77b4", "#d3d3d3"]
                )),
            ).properties(title="COSMIC Annotation"))
        if has_gnomad >= 0:
            gd = pl.DataFrame({
                "category": ["Has gnomAD AF", "No gnomAD AF"],
                "count": [has_gnomad, n - has_gnomad],
            }).to_pandas()
            charts.append(alt.Chart(gd).mark_arc(innerRadius=40).encode(
                theta=alt.Theta("count:Q"),
                color=alt.Color("category:N", scale=alt.Scale(
                    domain=["Has gnomAD AF", "No gnomAD AF"], range=["#ff7f0e", "#d3d3d3"]
                )),
            ).properties(title="gnomAD Annotation"))
        chart = alt.hconcat(*charts).properties(title="COSMIC and gnomAD Annotation Coverage") if len(charts) == 2 else charts[0]
    _save_chart(chart, "07_cosmic_gnomad", output_dir)
    return chart


def plot_caller_agreement_matrix(df, output_dir: str):
    """6×6 pairwise caller agreement matrix heatmap (count queries only)."""
    callers = ["DNA_mutect2", "DNA_deepsomatic", "DNA_strelka",
               "RNA_mutect2", "RNA_deepsomatic", "RNA_strelka"]
    gt_cols = [f"{c}_GT" for c in callers]
    existing = [c for c in gt_cols if c in df.columns]
    if len(existing) < 2:
        return

    n = _count_rows(df)
    rows = []
    for c1 in callers:
        col1 = f"{c1}_GT"
        if col1 not in df.columns:
            continue
        for c2 in callers:
            col2 = f"{c2}_GT"
            if col2 not in df.columns:
                continue
            both_valid = _count_rows(
                df.filter(
                    pl.col(col1).is_not_null() & pl.col(col2).is_not_null()
                    & ~pl.col(col1).is_in(["./.", "./.", "."])
                    & ~pl.col(col2).is_in(["./.", "./.", "."])
                )
            )
            pct = both_valid / n * 100 if n > 0 else 0
            rows.append({"caller_1": c1.replace("DNA_", "D_").replace("RNA_", "R_"),
                         "caller_2": c2.replace("DNA_", "D_").replace("RNA_", "R_"),
                         "pct": pct})
    pdf = pl.DataFrame(rows).to_pandas()
    chart = alt.Chart(pdf).mark_rect().encode(
        x=alt.X("caller_1:N", title=None),
        y=alt.Y("caller_2:N", title=None),
        color=alt.Color("pct:Q", title="% Both Valid GT", scale=alt.Scale(scheme="blues")),
        tooltip=["caller_1", "caller_2", "pct"],
    ).properties(title="Caller GT Availability Matrix (%)")
    _save_chart(chart, "28_caller_agreement", output_dir)
    return chart


# ═══════════════════════════════════════════════════════════════════════════
# iter_rows charts — collect needed columns, then iterate in Python.
# Memory: ~2 GB for GT columns (58M rows × 4-5 short strings).
# ═══════════════════════════════════════════════════════════════════════════

def plot_gt_concordance(df, output_dir: str, group_col: str = "set_number"):
    """Chart 6: GT concordance among 4 callers with GT fields, per group."""
    from collections import Counter

    gt_cols = [
        "DNA_mutect2_GT", "RNA_mutect2_GT",
        "DNA_deepsomatic_GT", "RNA_deepsomatic_GT",
    ]
    existing = [c for c in gt_cols if c in df.columns]
    if len(existing) < 2:
        return

    # Determine groups: if group_col exists, compute per-group; else global
    has_groups = group_col in df.columns
    cols_to_get = existing + ([group_col] if has_groups else [])
    pdf = _maybe_collect(df.select(cols_to_get))

    if has_groups:
        # Per-group concordance
        rows = []
        for group_val in pdf[group_col].unique().to_list():
            group_df = pdf.filter(pl.col(group_col) == group_val)
            agree_counts = {2: 0, 3: 0, 4: 0}
            for row in group_df.iter_rows():
                gts = [g for g in row[:-1] if g is not None and g not in ("./.", "./.", ".")]
                if len(gts) >= 2:
                    best = Counter(gts).most_common(1)[0][1]
                    if best >= 2:
                        agree_counts[best] = agree_counts.get(best, 0) + 1
            for level in [2, 3, 4]:
                rows.append({"group": str(group_val), "agreement_level": str(level),
                             "count": agree_counts[level]})
        if not rows:
            return
        counts_df = pl.DataFrame(rows).to_pandas()
        group_title = group_col.replace("_", " ").title()
        chart = alt.Chart(counts_df).mark_bar().encode(
            x=alt.X("agreement_level:N", title="Number of Callers Agreeing"),
            y=alt.Y("count:Q", title="Number of Variants"),
            color=alt.Color("agreement_level:N"),
            column=alt.Column("group:N", title=group_title),
        ).properties(title=f"GT Concordance by {group_title}")
    else:
        agree_counts = {2: 0, 3: 0, 4: 0}
        n_total = 0
        for row in pdf.iter_rows():
            gts = [g for g in row if g is not None and g not in ("./.", "./.", ".")]
            if len(gts) >= 2:
                best = Counter(gts).most_common(1)[0][1]
                if best >= 2:
                    agree_counts[best] = agree_counts.get(best, 0) + 1
                    n_total += 1
        if n_total == 0:
            return
        counts_df = pl.DataFrame({
            "agreement_level": ["2", "3", "4"],
            "count": [agree_counts[2], agree_counts[3], agree_counts[4]],
        }).to_pandas()
        chart = alt.Chart(counts_df).mark_bar().encode(
            x=alt.X("agreement_level:N", title="Number of Callers Agreeing on GT"),
            y=alt.Y("count:Q", title="Number of Variants"),
        ).properties(title=f"GT Concordance Among 4 Callers (n={n_total} with ≥2 agreement)")
    _save_chart(chart, "06_gt_concordance", output_dir)
    return chart


def plot_gt_concordance_per_tier(df, output_dir: str, facet_col: str = None):
    """GT concordance faceted by caller tier."""
    gt_cols = [
        "DNA_mutect2_GT", "RNA_mutect2_GT",
        "DNA_deepsomatic_GT", "RNA_deepsomatic_GT",
    ]
    existing_gt = [c for c in gt_cols if c in df.columns]
    if len(existing_gt) < 2 or "caller_tier" not in df.columns:
        return

    from collections import Counter

    cols_to_collect = existing_gt + ["caller_tier"]
    if facet_col and facet_col in df.columns:
        cols_to_collect.append(facet_col)
    pdf = _maybe_collect(df.select(cols_to_collect))

    n_gt = len(existing_gt)
    agree_counts: dict[tuple, int] = {}
    for row in pdf.iter_rows():
        # row[:n_gt] are GT columns; row[n_gt] is caller_tier; row[n_gt+1] (if present) is facet_col
        gts = [g for g in row[:n_gt]
               if g is not None and g not in ("./.", "./.", ".")]
        if len(gts) < 2:
            continue
        best = Counter(gts).most_common(1)[0][1]
        if best >= 2:
            key = (row[n_gt], best)  # caller_tier at index n_gt
            agree_counts[key] = agree_counts.get(key, 0) + 1

    if not agree_counts:
        return

    result = pl.DataFrame(
        [{"caller_tier": k[0], "agreement": k[1], "count": v} for k, v in agree_counts.items()]
    ).to_pandas()
    result["agreement"] = result["agreement"].astype(str)

    col_enc = alt.Column("caller_tier:N", title="Caller Tier")
    if facet_col and facet_col in result.columns:
        result[facet_col] = result[facet_col].astype(str)
        col_enc = alt.Column(f"{facet_col}:N", title=facet_col.replace("_", " ").title())
    chart = alt.Chart(result).mark_bar().encode(
        x=alt.X("agreement:N", title="Number of Callers Agreeing"),
        y=alt.Y("count:Q", title="Count"),
        color=alt.Color("agreement:N"),
        column=col_enc,
    ).properties(title="GT Concordance per Caller Tier")
    _save_chart(chart, "16_gt_concordance_per_tier", output_dir)
    return chart


# ═══════════════════════════════════════════════════════════════════════════
# Small-data charts — receive eager DataFrames (not lazy scans).
# These use sample_stats_df, sample_tier_df, bam_stats_df, or validation
# report — all small (hundreds of rows). No sampling needed.
# ═══════════════════════════════════════════════════════════════════════════

def plot_per_sample_distribution(sample_stats_df, output_dir: str, top_n: int = None):
    """Chart 8: Per-sample variant count, horizontal bar chart. All samples.

    If set_number column exists, facets by set using row encoding.
    """
    if sample_stats_df is None or (hasattr(sample_stats_df, 'is_empty') and sample_stats_df.is_empty()):
        return
    if "total_variants" not in sample_stats_df.columns or "sample_id" not in sample_stats_df.columns:
        return
    select_cols = ["sample_id", "total_variants"]
    has_set = "set_number" in sample_stats_df.columns
    if has_set:
        select_cols.append("set_number")
    pdf = sample_stats_df.select(select_cols).sort(
        "total_variants", descending=True
    )
    if top_n:
        pdf = pdf.head(top_n)
    if has_set:
        pdf = pdf.with_columns(pl.col("set_number").cast(pl.Utf8))
    pdf = pdf.to_pandas()
    n_samples = len(pdf)
    enc = {
        "y": alt.Y("sample_id:N", title="Sample ID", sort=None),
        "x": alt.X("total_variants:Q", title="Total Variants per Sample"),
        "tooltip": ["sample_id", "total_variants"],
    }
    if has_set:
        enc["row"] = alt.Row("set_number:N", title="Set")
    chart = alt.Chart(pdf).mark_bar().encode(**enc).properties(
        title=f"Per-Sample Variant Counts (n={n_samples})",
        height=max(300, n_samples * 12)  # dynamic height for scroll
    )
    _save_chart(chart, "08_per_sample_distribution", output_dir)
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


def plot_bam_metrics_bars(bam_stats_df, output_dir: str, top_n: int = 20):
    """BAM metrics: grouped bar chart of per-sample reads for DN/DT/RT, faceted by set."""
    if bam_stats_df is None or (hasattr(bam_stats_df, 'is_empty') and bam_stats_df.is_empty()):
        return
    needed = ["sample_id", "bam_type", "total_reads"]
    if not all(c in bam_stats_df.columns for c in needed):
        return
    # Show all samples, faceted by set_number if available
    select_cols = ["sample_id", "bam_type", "total_reads", "mapped_reads"]
    has_set = "set_number" in bam_stats_df.columns
    if has_set:
        select_cols.append("set_number")
    pdf = bam_stats_df.select(select_cols).to_pandas()
    enc = {
        "x": alt.X("sample_id:N", title="Sample", axis=alt.Axis(labelAngle=-45, labelLimit=120)),
        "y": alt.Y("total_reads:Q", title="Total Reads"),
        "color": alt.Color("bam_type:N", title="BAM Type"),
        "xOffset": "bam_type:N",
    }
    if has_set:
        pdf["set_number"] = pdf["set_number"].astype(str)
        enc["row"] = alt.Row("set_number:N", title="Set")
    chart = alt.Chart(pdf).mark_bar().encode(**enc).properties(
        title="Per-Sample BAM Read Counts — DN/DT/RT")
    _save_chart(chart, "20_bam_metrics", output_dir)
    return chart


def plot_per_sample_tier_distribution(sample_tier_df, output_dir: str):
    """Stacked bar chart: per-sample per-tier variant counts.

    If set_number column exists, facets by set using row encoding.
    """
    if sample_tier_df is None or (hasattr(sample_tier_df, 'is_empty') and sample_tier_df.is_empty()):
        return
    if "sample_id" not in sample_tier_df.columns or "final_tier" not in sample_tier_df.columns:
        return
    select_cols = ["sample_id", "final_tier", "n_variants"]
    has_set = "set_number" in sample_tier_df.columns
    if has_set:
        select_cols.append("set_number")
    pdf = sample_tier_df.select(select_cols)
    if has_set:
        pdf = pdf.with_columns(pl.col("set_number").cast(pl.Utf8))
    pdf = pdf.to_pandas()
    enc = {
        "y": alt.Y("sample_id:N", title="Sample ID", sort=None),
        "x": alt.X("n_variants:Q", title="Variants"),
        "color": alt.Color("final_tier:N", title="Tier"),
    }
    if has_set:
        enc["row"] = alt.Row("set_number:N", title="Set")
    chart = alt.Chart(pdf).mark_bar().encode(**enc).properties(
        title="Per-Sample Per-Tier Variant Distribution")
    _save_chart(chart, "22_sample_tier_dist", output_dir)
    return chart


def plot_bam_coverage_violin(df, output_dir: str, color_col: str = None):
    """BAM pileup coverage depth distribution across variant positions.

    Shows per-position DP (total, REF, ALT) for each BAM type (DN/DT/RT)
    as violin plots. Handles partial BAM type availability (e.g., DT+RT
    present but DN absent). Data is sampled to 10K rows to keep memory bounded.
    """
    dp_cols = []
    for bt in ["DN", "DT", "RT"]:
        for suffix in ["DP", "REF_DP", "ALT_DP"]:
            col = f"BAM_{bt}_{suffix}"
            if col in df.columns:
                dp_cols.append(col)

    if len(dp_cols) < 1:
        return

    select_cols = dp_cols[:]
    if color_col and color_col in df.columns:
        select_cols.append(color_col)
    melted = df.select(select_cols).unpivot(
        index=[color_col] if color_col and color_col in df.columns else [],
        variable_name="metric", value_name="depth"
    ).drop_nulls()
    melted = melted.filter(pl.col("depth") <= 2000)
    sampled = _sample_if_large(melted, max_rows=10000)

    if sampled is None or (hasattr(sampled, 'is_empty') and sampled.is_empty()):
        return
    if hasattr(sampled, 'height') and sampled.height == 0:
        return
    sampled = sampled.to_pandas()
    if len(sampled) == 0:
        return

    enc = {"x": alt.X("depth:Q", title="Depth at Variant Position"),
           "y": alt.Y("density:Q", title="Density"),
           "color": alt.Color("metric:N", title="BAM Metric")}
    if color_col and color_col in sampled.columns:
        enc["column"] = alt.Column(f"{color_col}:N")
    chart = alt.Chart(sampled).transform_density(
        "depth", groupby=["metric"]
    ).mark_area(opacity=0.5).encode(**enc).properties(
        title="BAM Pileup Depth Distribution at Variant Positions (depth capped at 2000, sampled)")
    _save_chart(chart, "21_bam_coverage_violin", output_dir)
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
    pdf = sample_tier_df.select(["final_tier", vaf_col]).to_pandas()
    chart = alt.Chart(pdf).mark_boxplot().encode(
        x=alt.X("final_tier:N", title="Tier"),
        y=alt.Y(f"{vaf_col}:Q", title=f"Mean {vaf_col}"),
        color=alt.Color("final_tier:N"),
    ).properties(title="Per-Tier VAF Distribution Across Samples")
    _save_chart(chart, "23_per_tier_vaf", output_dir)
    return chart


# ═══════════════════════════════════════════════════════════════════════════
# Threshold analysis charts (4.3)
# ═══════════════════════════════════════════════════════════════════════════

def plot_vaf_threshold_sweep(sweep_df, output_dir: str):
    """Chart 31: VAF threshold sweep — multi-line retention% vs threshold per caller.

    Args:
        sweep_df: Output of compute_vaf_threshold_sweep() — long-form DataFrame
                  with caller, threshold, pct_retained columns.
    """
    if sweep_df is None or (hasattr(sweep_df, 'is_empty') and sweep_df.is_empty()):
        return
    needed = ["caller", "threshold", "pct_retained"]
    if not all(c in sweep_df.columns for c in needed):
        return
    # Filter to overall (non-classification) rows
    pdf = sweep_df.filter(pl.col("classification").is_null() if "classification" in sweep_df.columns else pl.lit(True))
    if isinstance(pdf, pl.LazyFrame):
        pdf = pdf.collect()
    if pdf.is_empty():
        return
    pdf = pdf.to_pandas()
    chart = alt.Chart(pdf).mark_line(point=True).encode(
        x=alt.X("threshold:Q", title="VAF Threshold"),
        y=alt.Y("pct_retained:Q", title="% Variants Retained"),
        color=alt.Color("caller:N", title="Caller"),
    ).properties(title="VAF Threshold Sweep: Retention % vs Threshold per Caller")
    _save_chart(chart, "31_vaf_threshold_sweep", output_dir)
    return chart


def plot_dp_threshold_sweep(sweep_df, output_dir: str):
    """Chart 35: DP threshold sweep — multi-line retention% vs threshold per caller/metric.

    Args:
        sweep_df: Output of compute_dp_threshold_sweep() — long-form DataFrame
                  with caller, metric, threshold, pct_retained columns.
    """
    if sweep_df is None or (hasattr(sweep_df, 'is_empty') and sweep_df.is_empty()):
        return
    needed = ["caller", "metric", "threshold", "pct_retained"]
    if not all(c in sweep_df.columns for c in needed):
        return
    if hasattr(sweep_df, 'is_empty') and sweep_df.is_empty():
        return
    pdf = sweep_df.to_pandas()
    chart = alt.Chart(pdf).mark_line(point=True).encode(
        x=alt.X("threshold:Q", title="Depth Threshold"),
        y=alt.Y("pct_retained:Q", title="% Variants Retained"),
        color=alt.Color("caller:N", title="Caller", scale=alt.Scale(scheme="category10")),
        column=alt.Column("metric:N", title="DP Metric"),
    ).properties(title="DP Threshold Sweep: Retention % vs Threshold per Caller")
    _save_chart(chart, "35_dp_threshold_sweep", output_dir)
    return chart


def plot_caller_concordance_vs_vaf(df, output_dir: str, color_col: str = None):
    """Chart 32: Caller concordance vs VAF — box plot by # supporting callers.

    Shows how the number of supporting callers relates to VAF.
    """
    if "N_SUPPORT_CALLERS" not in df.columns or "DNA_VAF_mean" not in df.columns:
        return
    select_cols = ["N_SUPPORT_CALLERS", "DNA_VAF_mean"]
    if color_col and color_col in df.columns:
        select_cols.append(color_col)
    pdf = df.select(select_cols).drop_nulls()
    pdf = _sample_if_large(pdf, max_rows=50000).to_pandas()
    pdf["N_SUPPORT_CALLERS"] = pdf["N_SUPPORT_CALLERS"].astype(int).astype(str)
    enc = {"x": alt.X("N_SUPPORT_CALLERS:N", title="Number of Supporting Callers"),
           "y": alt.Y("DNA_VAF_mean:Q", title="DNA Mean VAF"),
           "color": alt.Color("N_SUPPORT_CALLERS:N")}
    if color_col and color_col in pdf.columns:
        enc["column"] = alt.Column(f"{color_col}:N")
    chart = alt.Chart(pdf).mark_boxplot().encode(**enc).properties(title="Caller Concordance vs VAF")
    _save_chart(chart, "32_caller_concordance_vs_vaf", output_dir)
    return chart


def plot_filter_effectiveness_heatmap(matrix_df, output_dir: str):
    """Chart 33: Filter effectiveness heatmap — FILTER × Classification.

    Args:
        matrix_df: Output of compute_filter_effectiveness_matrix() with
                   classification, filter_flag, pct_flagged columns.
    """
    if matrix_df is None or (hasattr(matrix_df, 'is_empty') and matrix_df.is_empty()):
        return
    needed = ["classification", "filter_flag", "pct_flagged"]
    if not all(c in matrix_df.columns for c in needed):
        return
    pdf = matrix_df.to_pandas()
    chart = alt.Chart(pdf).mark_rect().encode(
        x=alt.X("filter_flag:N", title="Filter Flag"),
        y=alt.Y("classification:N", title="Classification"),
        color=alt.Color("pct_flagged:Q", title="% Flagged",
                        scale=alt.Scale(scheme="redyellowgreen", reverse=True)),
        tooltip=["classification", "filter_flag", "pct_flagged"],
    ).properties(title="Filter Effectiveness: % Flagged by Classification × Filter")
    _save_chart(chart, "33_filter_effectiveness_heatmap", output_dir)
    return chart


def plot_database_enrichment_by_tier(df, output_dir: str):
    """Chart 34: Database enrichment by tier — COSMIC/gnomAD % per CxDy tier.

    Shows what fraction of variants in each tier have COSMIC or gnomAD annotations.
    """
    if "final_tier" not in df.columns:
        return
    n = _count_rows(df)
    if n == 0:
        return

    # Build enrichment table per tier (single group_by scan, not per-tier _count_rows)
    agg_exprs = [pl.len().alias("n_total")]
    if "COSMIC_ID" in df.columns:
        agg_exprs.append(pl.col("COSMIC_ID").is_not_null().sum().alias("n_cosmic"))
    if "GNOMAD_AF" in df.columns:
        agg_exprs.append(pl.col("GNOMAD_AF").is_not_null().sum().alias("n_gnomad"))
    counts_df = _maybe_collect(df.group_by("final_tier").agg(agg_exprs))
    if counts_df.is_empty():
        return
    rows = []
    for row in counts_df.sort("final_tier").iter_rows():
        tier = row[0]
        tier_n = row[1]  # n_total
        if tier_n == 0:
            continue
        # Skip D=0 tiers: D=0 means "no database evidence" by design.
        # Including them inflates the chart with zero-value bars that
        # confuse readers into thinking data is missing.
        if tier.endswith("D0"):
            continue
        r = {"tier": tier, "n_total": tier_n}
        idx = 2
        if "COSMIC_ID" in df.columns:
            r["cosmic_pct"] = round(row[idx] / tier_n * 100, 2)
            idx += 1
        if "GNOMAD_AF" in df.columns:
            r["gnomad_pct"] = round(row[idx] / tier_n * 100, 2)
            idx += 1
        rows.append(r)

    if not rows:
        return
    pdf = pl.DataFrame(rows).to_pandas()

    charts = []
    for db, col in [("COSMIC", "cosmic_pct"), ("gnomAD", "gnomad_pct")]:
        if col not in pdf.columns:
            continue
        c = alt.Chart(pdf).mark_bar().encode(
            x=alt.X("tier:N", title="Tier (CxDy)"),
            y=alt.Y(f"{col}:Q", title=f"% with {db} Annotation"),
            color=alt.Color("tier:N"),
        ).properties(title=f"{db} Annotation by Tier")
        charts.append(c)

    if len(charts) == 2:
        chart = alt.hconcat(*charts).properties(
            title="Database Annotation Enrichment by Tier (D=0 tiers excluded — no DB evidence by design)")
    elif charts:
        chart = charts[0].properties(
            title=f"Database Annotation Enrichment by Tier (D=0 tiers excluded)")
    else:
        return
    _save_chart(chart, "34_database_enrichment_by_tier", output_dir)
    return chart


# ═══════════════════════════════════════════════════════════════════════════
# New chart functions (Phase 4.4)
# ═══════════════════════════════════════════════════════════════════════════

def plot_dp_distribution(df, output_dir: str, color_col: str = None):
    """Chart 12: Per-caller DP distribution, box+violin overlay (sampled to 50K)."""
    caller_dp_cols = [
        f"{c}_DP" for c in [
            "DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic",
            "RNA_deepsomatic", "DNA_strelka", "RNA_strelka",
        ]
    ]
    existing = [c for c in caller_dp_cols if c in df.columns]
    if not existing:
        return
    select_cols = existing[:]
    if color_col and color_col in df.columns:
        select_cols.append(color_col)
    melted = df.select(select_cols).unpivot(
        index=[color_col] if color_col and color_col in df.columns else [],
        variable_name="caller", value_name="DP"
    ).drop_nulls()
    pdf = _sample_if_large(melted, max_rows=50000).to_pandas()
    pdf["caller"] = pdf["caller"].str.replace("_DP", "")
    enc = {"x": alt.X("caller:N", title="Caller", axis=alt.Axis(labelAngle=-45)),
           "y": alt.Y("DP:Q", title="Read Depth (capped at 2000)", scale=alt.Scale(domain=[0, 2000]))}
    if color_col and color_col in pdf.columns:
        enc["color"] = alt.Color(f"{color_col}:N")
        enc["column"] = alt.Column(f"{color_col}:N")
    is_faceted = "column" in enc
    if not is_faceted:
        enc["color"] = alt.Color("caller:N")
        violin = alt.Chart(pdf).transform_density("DP", groupby=["caller"]
            ).mark_area(opacity=0.3).encode(
                x=alt.X("caller:N", title="Caller", axis=alt.Axis(labelAngle=-45)),
                y=alt.Y("DP:Q", title="Read Depth (capped at 2000)", scale=alt.Scale(domain=[0, 2000])),
                color=alt.Color("caller:N"))
        box = alt.Chart(pdf).mark_boxplot(size=30).encode(**enc)
        chart = (violin + box).properties(title="DP Distribution per Caller (capped at 2000)")
    else:
        chart = alt.Chart(pdf).mark_boxplot(size=30).encode(**enc).properties(
            title="DP Distribution per Caller (capped at 2000)")
    _save_chart(chart, "12_dp_distribution", output_dir)
    return chart


def plot_mean_vaf_per_group(summary_df, output_dir: str, group_col: str):
    """Chart 35: Mean VAF per group from wise summary data (bar chart)."""
    if summary_df is None or (hasattr(summary_df, 'is_empty') and summary_df.is_empty()):
        return
    needed_dna = ["mean_dna_vaf", group_col]
    needed_rna = ["mean_rna_vaf", group_col]
    has_dna = all(c in summary_df.columns for c in needed_dna)
    has_rna = all(c in summary_df.columns for c in needed_rna)
    if not has_dna and not has_rna:
        return
    pdf = summary_df.select(
        [group_col] + ([c for c in ["mean_dna_vaf", "mean_rna_vaf"] if c in summary_df.columns])
    ).to_pandas()
    pdf[group_col] = pdf[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()
    charts = []
    for col, label, color in [("mean_dna_vaf", "DNA Mean VAF", "#1f77b4"),
                               ("mean_rna_vaf", "RNA Mean VAF", "#ff7f0e")]:
        if col not in pdf.columns:
            continue
        c = alt.Chart(pdf).mark_bar().encode(
            x=alt.X(f"{group_col}:N", title=group_title),
            y=alt.Y(f"{col}:Q", title=label),
            color=alt.value(color),
        ).properties(title=f"{label} by {group_title}")
        charts.append(c)
    if not charts:
        return
    chart = alt.hconcat(*charts).properties(title=f"Mean VAF by {group_title}") if len(charts) == 2 else charts[0]
    _save_chart(chart, "40_mean_vaf_per_group", output_dir)
    return chart


def plot_mean_dp_per_group(summary_df, output_dir: str, group_col: str):
    """Chart 36: Mean DP per group from wise summary data (bar chart)."""
    if summary_df is None or (hasattr(summary_df, 'is_empty') and summary_df.is_empty()):
        return
    needed_dna = ["mean_dna_dp", group_col]
    needed_rna = ["mean_rna_dp", group_col]
    has_dna = all(c in summary_df.columns for c in needed_dna)
    has_rna = all(c in summary_df.columns for c in needed_rna)
    if not has_dna and not has_rna:
        return
    pdf = summary_df.select(
        [group_col] + ([c for c in ["mean_dna_dp", "mean_rna_dp"] if c in summary_df.columns])
    ).to_pandas()
    pdf[group_col] = pdf[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()
    charts = []
    for col, label, color in [("mean_dna_dp", "DNA Mean DP", "#1f77b4"),
                               ("mean_rna_dp", "RNA Mean DP", "#ff7f0e")]:
        if col not in pdf.columns:
            continue
        c = alt.Chart(pdf).mark_bar().encode(
            x=alt.X(f"{group_col}:N", title=group_title),
            y=alt.Y(f"{col}:Q", title=label),
            color=alt.value(color),
        ).properties(title=f"{label} by {group_title}")
        charts.append(c)
    if not charts:
        return
    chart = alt.hconcat(*charts).properties(title=f"Mean DP by {group_title}") if len(charts) == 2 else charts[0]
    _save_chart(chart, "41_mean_dp_per_group", output_dir)
    return chart


def plot_n_support_callers_dist(df, output_dir: str, group_col: str = "set_number"):
    """Chart 37: N_SUPPORT_CALLERS distribution histogram per group."""
    if "N_SUPPORT_CALLERS" not in df.columns or group_col not in df.columns:
        return
    pdf = df.select([group_col, "N_SUPPORT_CALLERS"]).drop_nulls()
    pdf = _sample_if_large(pdf, max_rows=50000).to_pandas()
    pdf[group_col] = pdf[group_col].astype(str)
    pdf["N_SUPPORT_CALLERS"] = pdf["N_SUPPORT_CALLERS"].astype(int)
    group_title = group_col.replace("_", " ").title()
    chart = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("N_SUPPORT_CALLERS:O", title="Number of Supporting Callers"),
        y=alt.Y("count()", title="Number of Variants"),
        color=alt.Color(f"{group_col}:N"),
        column=alt.Column(f"{group_col}:N", title=group_title),
    ).properties(title=f"Caller Support Distribution by {group_title}")
    _save_chart(chart, "42_n_support_callers_dist", output_dir)
    return chart


def plot_caller_tier_heatmap(df, output_dir: str, facet_col: str = None):
    """Chart 38: Caller × Tier detection matrix heatmap."""
    caller_vaf_cols = [
        f"{c}_VAF" for c in [
            "DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic",
            "RNA_deepsomatic", "DNA_strelka", "RNA_strelka",
        ]
    ]
    existing = [c for c in caller_vaf_cols if c in df.columns]
    if not existing or "final_tier" not in df.columns:
        return

    rows = []
    for vaf_col in existing:
        caller = vaf_col.replace("_VAF", "")
        tiers = _maybe_collect(df.select("final_tier").unique())["final_tier"].to_list()
        for tier in sorted(tiers):
            n_total = _count_rows(df.filter(pl.col("final_tier") == tier))
            n_detected = _count_rows(df.filter(
                (pl.col("final_tier") == tier) & pl.col(vaf_col).is_not_null()
            ))
            row = {"caller": caller, "tier": tier, "pct_detected": round(
                n_detected / n_total * 100, 2) if n_total > 0 else 0}
            if facet_col and facet_col in df.columns:
                for fv in _maybe_collect(df.select(facet_col).unique())[facet_col].to_list():
                    fv_total = _count_rows(df.filter(
                        (pl.col("final_tier") == tier) & (pl.col(facet_col) == fv)))
                    fv_detected = _count_rows(df.filter(
                        (pl.col("final_tier") == tier) & (pl.col(facet_col) == fv)
                        & pl.col(vaf_col).is_not_null()))
                    row_fv = {"caller": caller, "tier": tier, "facet_val": str(fv),
                              "pct_detected": round(fv_detected / fv_total * 100, 2) if fv_total > 0 else 0}
                    rows.append(row_fv)
            else:
                rows.append(row)
        if facet_col and facet_col in df.columns:
            break  # Already computed with facet

    if not rows:
        return
    pdf = pl.DataFrame(rows).to_pandas()
    chart = alt.Chart(pdf).mark_rect().encode(
        x=alt.X("tier:N", title="Tier (CxDy)"),
        y=alt.Y("caller:N", title="Caller"),
        color=alt.Color("pct_detected:Q", title="% Detected", scale=alt.Scale(scheme="blues")),
        tooltip=["caller", "tier", "pct_detected"],
    )
    if facet_col and "facet_val" in pdf.columns:
        chart = chart.encode(column=alt.Column("facet_val:N"))
    chart = chart.properties(title="Caller Detection Rate by Tier")
    _save_chart(chart, "38_caller_tier_heatmap", output_dir)
    return chart


def plot_sample_overview_scatter(sample_stats_df, output_dir: str):
    """Chart 39: Per-sample overview — mean VAF vs mean DP scatter."""
    if sample_stats_df is None or (hasattr(sample_stats_df, 'is_empty') and sample_stats_df.is_empty()):
        return
    needed = ["sample_id", "total_variants"]
    vaf_col = "mean_dna_vaf_mean" if "mean_dna_vaf_mean" in sample_stats_df.columns else None
    dp_col = "mean_dna_dp_mean" if "mean_dna_dp_mean" in sample_stats_df.columns else None
    if not vaf_col or not dp_col:
        return
    pdf = sample_stats_df.select(
        [c for c in ["sample_id", "total_variants", "set_number", "disease", vaf_col, dp_col]
         if c in sample_stats_df.columns]
    ).to_pandas()
    color_col = "disease" if "disease" in pdf.columns else ("set_number" if "set_number" in pdf.columns else None)
    enc = {"x": alt.X(f"{vaf_col}:Q", title="Mean DNA VAF"),
           "y": alt.Y(f"{dp_col}:Q", title="Mean DNA Depth"),
           "size": alt.Size("total_variants:Q", title="Total Variants"),
           "tooltip": ["sample_id", "total_variants", vaf_col, dp_col]}
    if color_col:
        pdf[color_col] = pdf[color_col].astype(str)
        enc["color"] = alt.Color(f"{color_col}:N")
    chart = alt.Chart(pdf).mark_circle(opacity=0.7).encode(**enc).properties(
        title="Sample Overview: Mean VAF vs Mean Depth")
    _save_chart(chart, "39_sample_overview_scatter", output_dir)
    return chart


# ═══════════════════════════════════════════════════════════════════════════
# ML Threshold Guidance Charts (Section 7)
# ═══════════════════════════════════════════════════════════════════════════

def plot_filter_vaf_dp_heatmap(cross_tab_df, output_dir: str):
    """Chart 43: FILTER x VAF_bin x DP_bin heatmap for ML threshold guidance.

    Shows the joint distribution of classification, VAF bin, and DP bin
    as a heatmap faceted by classification.
    """
    if cross_tab_df is None or (hasattr(cross_tab_df, 'is_empty') and cross_tab_df.is_empty()):
        return
    needed = ["classification", "vaf_bin", "dp_bin", "count"]
    if not all(c in cross_tab_df.columns for c in needed):
        return
    pdf = cross_tab_df.to_pandas()
    # Order VAF and DP bins sensibly
    vaf_order = ["<0.01", "0.01-0.05", "0.05-0.10", "0.10-0.25", "0.25-0.50", "0.50-1.0"]
    dp_order = ["<10", "10-50", "50-100", "100-200", "200-500", "500+"]
    chart = alt.Chart(pdf).mark_rect().encode(
        x=alt.X("vaf_bin:N", title="VAF Bin", sort=vaf_order),
        y=alt.Y("dp_bin:N", title="DP Bin", sort=dp_order),
        color=alt.Color("count:Q", title="Count", scale=alt.Scale(scheme="blues", type="log")),
        facet=alt.Facet("classification:N", columns=3, title="Classification"),
        tooltip=["classification", "vaf_bin", "dp_bin", "count"],
    ).properties(title="FILTER x VAF x DP Cross-Tabulation", width=200, height=180)
    _save_chart(chart, "43_filter_vaf_dp_heatmap", output_dir)
    return chart


def plot_low_vaf_rna_support(low_vaf_df, output_dir: str):
    """Chart 44: Low VAF variants with strong RNA support — stacked bar by FILTER."""
    if low_vaf_df is None or (hasattr(low_vaf_df, 'is_empty') and low_vaf_df.is_empty()):
        return
    if "FILTER" not in low_vaf_df.columns or "n_variants" not in low_vaf_df.columns:
        return
    pdf = low_vaf_df.to_pandas()
    chart = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("FILTER:N", title="Classification"),
        y=alt.Y("n_variants:Q", title="Number of Variants"),
        color=alt.Color("FILTER:N", title="Classification",
                        scale=alt.Scale(domain=CLASSIFICATION_DOMAIN, range=CLASSIFICATION_COLORS)),
        tooltip=["FILTER", "n_variants", "mean_vaf"],
    ).properties(title="Low VAF (< 0.05) Variants with RNA Support (N_RNA_CALLERS >= 2)")
    _save_chart(chart, "44_low_vaf_rna_support", output_dir)
    return chart


# ═══════════════════════════════════════════════════════════════════════════
# FP Cross-Tabulation Charts (Section 8)
# ═══════════════════════════════════════════════════════════════════════════

def plot_fp_cross_tab_heatmap(fp_df, output_dir: str):
    """Chart 45: FP cross-tabulation heatmap — FILTER x N_SUPPORT_CALLERS for non-Somatic."""
    if fp_df is None or (hasattr(fp_df, 'is_empty') and fp_df.is_empty()):
        return
    needed = ["FILTER", "N_SUPPORT_CALLERS", "count"]
    if not all(c in fp_df.columns for c in needed):
        return
    pdf = fp_df.to_pandas()
    pdf["N_SUPPORT_CALLERS"] = pdf["N_SUPPORT_CALLERS"].astype(int).astype(str)
    chart = alt.Chart(pdf).mark_rect().encode(
        x=alt.X("N_SUPPORT_CALLERS:N", title="Number of Supporting Callers"),
        y=alt.Y("FILTER:N", title="Classification"),
        color=alt.Color("count:Q", title="Count", scale=alt.Scale(scheme="orangered", type="log")),
        tooltip=["FILTER", "N_SUPPORT_CALLERS", "count"],
    ).properties(title="FP Cross-Tabulation: Non-Somatic FILTER x Caller Support")
    _save_chart(chart, "45_fp_cross_tab_heatmap", output_dir)
    return chart


# ═══════════════════════════════════════════════════════════════════════════
# Somatic Modality Sub-Classification Charts (Section 9)
# ═══════════════════════════════════════════════════════════════════════════

def plot_somatic_modality_pie(modality_df, output_dir: str):
    """Chart 46: Somatic modality sub-classification — pie chart."""
    if modality_df is None or (hasattr(modality_df, 'is_empty') and modality_df.is_empty()):
        return
    if "somatic_modality" not in modality_df.columns or "n_variants" not in modality_df.columns:
        return
    modality_domain = ["MultiModality", "DNA_only", "RNA_only", "Weak", "Unknown"]
    modality_colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#8c564b"]
    pdf = modality_df.to_pandas()
    chart = alt.Chart(pdf).mark_arc(innerRadius=40).encode(
        theta=alt.Theta("n_variants:Q"),
        color=alt.Color("somatic_modality:N", title="Somatic Modality",
                        scale=alt.Scale(domain=modality_domain, range=modality_colors)),
        tooltip=["somatic_modality", "n_variants"],
    ).properties(title="Somatic Modality Sub-Classification")
    _save_chart(chart, "46_somatic_modality_pie", output_dir)
    return chart


def plot_somatic_modality_bars(modality_df, output_dir: str):
    """Chart 47: Somatic modality — bar chart with mean VAF/DP."""
    if modality_df is None or (hasattr(modality_df, 'is_empty') and modality_df.is_empty()):
        return
    if "somatic_modality" not in modality_df.columns or "n_variants" not in modality_df.columns:
        return
    pdf = modality_df.to_pandas()
    # Variants bar
    bars = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("somatic_modality:N", title="Somatic Modality"),
        y=alt.Y("n_variants:Q", title="Number of Variants"),
        color=alt.Color("somatic_modality:N", title="Modality"),
        tooltip=list(pdf.columns),
    ).properties(title="Somatic Modality Sub-Classification: Variant Counts")
    # VAF overlay if available
    subcharts = [bars]
    if "mean_dna_vaf" in pdf.columns:
        vaf_chart = alt.Chart(pdf).mark_bar().encode(
            x=alt.X("somatic_modality:N", title="Somatic Modality"),
            y=alt.Y("mean_dna_vaf:Q", title="Mean DNA VAF"),
            color=alt.Color("somatic_modality:N"),
        ).properties(title="Mean DNA VAF by Modality")
        subcharts.append(vaf_chart)
    chart = alt.hconcat(*subcharts) if len(subcharts) > 1 else bars
    _save_chart(chart, "47_somatic_modality_bars", output_dir)
    return chart

def _chart_section(fig, index: int) -> str:
    """Determine dashboard section for a chart based on its title."""
    title = ""
    try:
        title = fig.title if hasattr(fig, 'title') else ""
        title = str(title).lower()
    except Exception:
        pass
    if any(kw in title for kw in ["vc", "classification", "variant type", "ti/tv", "cosmic", "filter dist"]):
        return "overview"
    if any(kw in title for kw in ["tier", "caller support", "c1", "gt concordance per"]):
        return "tier"
    if any(kw in title for kw in ["bam", "coverage", "validation", "ref_alt", "vaf v", "dp v", "dp per"]):
        return "bam"
    if any(kw in title for kw in ["per-sample", "sample", "cross"]):
        return "persample"
    if any(kw in title for kw in ["threshold", "vaf sweep", "concordance vs", "effectiveness", "enrichment",
                                    "cross-tabulation", "low vaf", "fp cross", "partition", "modality"]):
        return "threshold"
    return "overview"


def _extract_body_content(html: str) -> str:
    """Extract inner content between <body> and </body> from a full HTML document."""
    import re
    body_start = re.search(r'<body[^>]*>', html)
    body_end = html.rfind('</body>')
    if body_start and body_end != -1:
        return html[body_start.end():body_end].strip()
    return html


def plot_bam_dp_distribution(df, output_dir: str, color_col: str = None):
    """Chart 36: BAM pileup DP distribution — violin + box per BAM type and metric.

    Shows per-position DP (total, REF, ALT) for DT and RT BAM types
    as box+violin overlay. Data is sampled to 10K rows.
    """
    dp_cols = []
    for bt in ["DT", "RT"]:
        for suffix, label in [("DP", "BAM DP"), ("REF_DP", "BAM REF DP"), ("ALT_DP", "BAM ALT DP")]:
            col = f"BAM_{bt}_{suffix}"
            if col in df.columns:
                dp_cols.append((col, f"{bt} {label}"))

    if not dp_cols:
        return

    cols, labels = zip(*dp_cols)
    select_cols = list(cols)
    melted = df.select(select_cols).unpivot(
        variable_name="metric", value_name="depth"
    ).drop_nulls()
    # Map column names to readable labels
    label_map = dict(dp_cols)
    melted = melted.with_columns(
        pl.col("metric").replace_strict(label_map, default=pl.col("metric"))
    )
    sampled = _sample_if_large(melted, max_rows=10000)

    if sampled is None or (hasattr(sampled, 'is_empty') and sampled.is_empty()):
        return
    if hasattr(sampled, 'height') and sampled.height == 0:
        return
    sampled = sampled.to_pandas()
    if len(sampled) == 0:
        return

    enc = {"x": alt.X("metric:N", title="BAM Metric", axis=alt.Axis(labelAngle=-45)),
           "y": alt.Y("depth:Q", title="Depth at Variant Position", scale=alt.Scale(type="symlog"))}
    if color_col and color_col in sampled.columns:
        enc["color"] = alt.Color(f"{color_col}:N")
        enc["column"] = alt.Column(f"{color_col}:N")
    is_faceted = "column" in enc
    if not is_faceted:
        # Boxplot-only approach with symlog scale (violin + log produces distorted density)
        enc["color"] = alt.Color("metric:N", scale=alt.Scale(scheme="category10"))
        box = alt.Chart(sampled).mark_boxplot(size=30).encode(**enc)
        chart = box.properties(title="BAM Pileup DP Distribution per BAM Type (symlog scale)")
    else:
        chart = alt.Chart(sampled).mark_boxplot(size=30).encode(**enc).properties(
            title="BAM Pileup DP Distribution per BAM Type (symlog scale)")
    _save_chart(chart, "36_bam_dp_distribution", output_dir)
    return chart


def plot_per_tier_dp_boxplot(df, output_dir: str):
    """Chart 37: Per-tier DNA DP boxplot (mirrors plot_per_tier_vaf_boxplot)."""
    if "DNA_DP_mean" not in df.columns or "final_tier" not in df.columns:
        return
    pdf = df.select(["final_tier", "DNA_DP_mean"]).drop_nulls(subset=["DNA_DP_mean"])
    pdf = _sample_if_large(pdf, max_rows=50000).to_pandas()
    chart = alt.Chart(pdf).mark_boxplot().encode(
        x=alt.X("final_tier:N", title="Tier"),
        y=alt.Y("DNA_DP_mean:Q", title="DNA Mean DP (capped at 2000)", scale=alt.Scale(domain=[0, 2000])),
        color=alt.Color("final_tier:N", scale=alt.Scale(scheme="category10")),
    ).properties(title="DNA DP Distribution per Tier (capped at 2000)")
    _save_chart(chart, "37_per_tier_dp", output_dir)
    return chart


def generate_dashboard(figs: list, output_dir: str):
    """Combine all figures into a single dashboard HTML with 4 sections."""
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
        '<div class="toc">',
        '<a href="#overview">Overview</a>',
        '<a href="#tier">Tier Analysis</a>',
        '<a href="#bam">BAM & Validation</a>',
        '<a href="#persample">Per-Sample</a>',
        '</div>',
    ]

    section_order = [
        ("overview", "Overview"),
        ("tier", "Tier Analysis"),
        ("threshold", "Threshold Analysis"),
        ("bam", "BAM & Validation"),
        ("persample", "Per-Sample"),
        ("validation", "Validation"),
    ]
    current_section = None

    for i, fig in enumerate(figs):
        if fig is None:
            continue
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
