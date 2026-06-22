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

CALLER_DOMAIN = ["DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic", "RNA_deepsomatic", "DNA_strelka", "RNA_strelka"]
CALLER_COLORS = ["#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a"]

RESCUE_DOMAIN = ["YES", "NO"]
RESCUE_COLORS = ["#2ca02c", "#d62728"]

VARIANT_TYPE_DOMAIN = ["SNV", "INS", "DEL", "MNV"]
VARIANT_TYPE_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

SOMATIC_MODALITY_DOMAIN = ["MultiModality", "DNA_only", "RNA_only", "Weak", "Unknown"]
SOMATIC_MODALITY_COLORS = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#8c564b"]

BAM_TYPE_DOMAIN = ["DN", "DT", "RT"]
BAM_TYPE_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c"]

MODALITY_DOMAIN = ["DNA", "RNA"]
MODALITY_COLORS = ["#1f77b4", "#ff7f0e"]

AGREEMENT_DOMAIN = ["2", "3", "4"]
AGREEMENT_COLORS = ["#ff7f0e", "#2ca02c", "#1f77b4"]

# Registry: entity name → (domain, colors). Dynamic entities use scheme fallback.
_COLOR_REGISTRY = {
    "caller": (CALLER_DOMAIN, CALLER_COLORS),
    "FILTER": (CLASSIFICATION_DOMAIN, CLASSIFICATION_COLORS),
    "RESCUED": (RESCUE_DOMAIN, RESCUE_COLORS),
    "variant_type": (VARIANT_TYPE_DOMAIN, VARIANT_TYPE_COLORS),
    "somatic_modality": (SOMATIC_MODALITY_DOMAIN, SOMATIC_MODALITY_COLORS),
    "bam_type": (BAM_TYPE_DOMAIN, BAM_TYPE_COLORS),
    "modality": (MODALITY_DOMAIN, MODALITY_COLORS),
    "agreement_level": (AGREEMENT_DOMAIN, AGREEMENT_COLORS),
    "agreement": (AGREEMENT_DOMAIN, AGREEMENT_COLORS),
    "modality_evidence_caller": (
        ["cross_modality", "dna_confident", "rna_rescued", "low_confidence"],
        ["#2ca02c", "#1f77b4", "#ff7f0e", "#d62728"]
    ),
}


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
# Unified helpers — colors, faceting, text, clipping, scales
# ═══════════════════════════════════════════════════════════════════════════

def _color_scale(entity: str) -> alt.Scale:
    """Return a unified alt.Scale for the given entity type.

    Uses _COLOR_REGISTRY for known entities, falls back to category10 scheme.
    """
    if entity in _COLOR_REGISTRY:
        domain, colors = _COLOR_REGISTRY[entity]
        return alt.Scale(domain=domain, range=colors)
    return alt.Scale(scheme="category10")


def _apply_faceting(chart, group_col: str, columns: int = None):
    """Apply the right faceting pattern based on group_col.

    Rules:
      set_number         → 2 columns, x-independent (per-set samples), y-shared
      disease_normalized → 4 columns, y-shared
      sample_id          → 2 columns, x-independent, y-shared
      FILTER             → 3 columns, y-shared
      others             → 3 columns default
    """
    _FACET_RULES = {
        "set_number":         {"columns": 2, "resolve": {"x": "independent", "y": "shared"}},
        "disease_normalized": {"columns": 4, "resolve": {"y": "shared"}},
        "sample_id":          {"columns": 2, "resolve": {"x": "independent", "y": "shared"}},
        "FILTER":             {"columns": 3, "resolve": {"y": "shared"}},
    }
    rule = _FACET_RULES.get(group_col, {"columns": 3, "resolve": {}})
    n_cols = columns if columns is not None else rule["columns"]
    result = chart.facet(
        facet=alt.Facet(f"{group_col}:N"),
        columns=n_cols,
    )
    resolve = rule.get("resolve", {})
    if resolve:
        result = result.resolve_scale(**resolve)
    return result


def _add_heatmap_text(base, x_enc, y_enc, text_col: str, pdf, fontSize: int = 8,
                      fmt: str = ",d"):
    """Add auto-colored text overlay to a heatmap chart.

    Text is white on dark cells (count > median) and black on light cells.
    """
    median_val = max(float(pdf[text_col].median()), 1) if len(pdf) > 0 else 1
    text = base.mark_text(baseline="middle", fontSize=fontSize).encode(
        x=x_enc,
        y=y_enc,
        text=alt.Text(f"{text_col}:Q", format=fmt),
        color=alt.condition(
            f"datum.{text_col} > {median_val}",
            alt.value("white"),
            alt.value("black"),
        ),
    )
    return text


def _si(v):
    """Format a numeric count with SI prefixes: 1M→"1.5M", 1K→"1.5K", else integer."""
    if v >= 1e6:
        return f"{v / 1e6:.1f}M"
    if v >= 1e3:
        return f"{v / 1e3:.0f}K"
    return str(int(v))


def _luminance(hex_color: str) -> float:
    """Compute luminance of a hex color (0.299*R + 0.587*G + 0.521*B).

    Returns a float in [0, 1]; < 0.5 is considered dark.
    """
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i:i + 2], 16) / 255.0 for i in (0, 2, 4))
    return 0.299 * r + 0.587 * g + 0.521 * b


def _make_bar_text(pdf, x_field, count_col, stack=None, x_offset=None,
                   x_sort=None, show_pct=False, group_col=None,
                   dy=-8, fontSize=8, base=None, color_col=None,
                   overlap="hide"):
    """Create count text labels for a bar chart with auto-contrast and responsive sizing.

    Args:
        pdf:        pandas DataFrame with chart data.
        x_field:    Column name for x-axis (nominal).
        count_col:  Column name for count values (quantitative).
        stack:      Stack mode ("zero") for stacked bars, None for simple.
        x_offset:   Column name for grouped bars (xOffset).
        x_sort:     Sort order for x-axis (list or None).
        show_pct:   If True, show "count (pct%)" labels — requires group_col.
        group_col:  Column to group by when computing percentages.
        dy:         Vertical offset for text marks (overridden by overlap).
        fontSize:   Font size for text marks (overridden by auto-scale when None/0).
        base:       Optional alt.Chart(pdf) base — pass the SAME base used for
                    bars when the layered (bars + text) chart will be .facet()'d.
                    Altair requires all layers to share one data source for faceting.
        color_col:  Column name used for color encoding (bar fill). When provided,
                    luminance of the fill color is used for auto-contrast text color.
        overlap:    Label overlap handling: "hide" (default, skip when >30 cats),
                    "rotate" (angle=90, align="left", dy=-4),
                    "stagger" (alternate dy=-8 vs dy=-16).

    Returns:
        alt.Chart text layer ready to be layered with (bars + text).
    """
    n_cats = pdf[x_field].nunique()

    # ── Responsive fontSize based on x-axis category count ────────────────
    if fontSize is None or fontSize <= 0 or fontSize == 8:
        if n_cats <= 10:
            fontSize = 10
        elif n_cats <= 25:
            fontSize = 9
        elif n_cats <= 50:
            fontSize = 8
        else:
            fontSize = 7

    # ── Overlap: always render text, use responsive fontSize for density ───
    # Skipping labels (empty/invisible chart) breaks both .facet() on layered
    # charts (data source mismatch) and vl-convert PNG export (child_width
    # signal error). Responsive fontSize (above) handles label density.

    # ── Overlap: rotation or stagger overrides ────────────────────────────
    if overlap == "rotate":
        dy = -4
    elif overlap == "stagger":
        # dy will be set per-row below via alternation column
        pass

    # ── Build text label ──────────────────────────────────────────────────
    if show_pct and group_col and group_col in pdf.columns:
        total = pdf.groupby(group_col)[count_col].transform("sum")
        pct = (pdf[count_col] / total * 100).round(0).astype(int)
        pdf["_lbl"] = pdf[count_col].apply(_si) + " (" + pct.astype(str) + "%)"
        text_enc = {"text": "_lbl:N"}
    else:
        pdf["_lbl"] = pdf[count_col].apply(_si)
        text_enc = {"text": "_lbl:N"}

    # ── Auto-contrast text color from bar fill luminance ──────────────────
    if color_col and color_col in pdf.columns:
        # Map color_scale domains to hex values for luminance check.
        # Use the _COLOR_REGISTRY if present for this color_col.
        entity = color_col
        domain, colors_range = _COLOR_REGISTRY.get(entity, (None, None))
        if domain and colors_range:
            color_map = dict(zip(domain, colors_range))
            def _text_color(val):
                hex_c = str(color_map.get(val, "#333"))
                return "white" if _luminance(hex_c) < 0.5 else "black"
            pdf["_text_color"] = pdf[color_col].apply(_text_color)
        else:
            pdf["_text_color"] = "#333"
        text_color = alt.Color("_text_color:N", scale=None)
    else:
        text_color = alt.value("#333")

    # ── Build encoding ────────────────────────────────────────────────────
    enc = {"x": alt.X(f"{x_field}:N", sort=x_sort)}
    if stack:
        enc["y"] = alt.Y(f"{count_col}:Q", stack=stack)
    else:
        enc["y"] = alt.Y(f"{count_col}:Q")
    if x_offset:
        enc["xOffset"] = f"{x_offset}:N"
    enc.update(text_enc)
    enc["color"] = text_color

    # ── Overlap encodings ─────────────────────────────────────────────────
    text_kw = {"fontSize": fontSize}
    if overlap == "rotate":
        text_kw.update(angle=90, align="left", dy=dy)
    elif overlap == "stagger":
        # Alternating dy: even-indexed rows get dy=-8, odd-indexed get dy=-16
        pdf["_dy"] = [(-8 if i % 2 == 0 else -16) for i in range(len(pdf))]
        enc["yOffset"] = "_dy:Q"
    else:
        text_kw["dy"] = dy

    chart_source = base if base is not None else alt.Chart(pdf)
    return chart_source.mark_text(**text_kw).encode(**enc)


def _clip_dp(pdf, col: str, cap: int = 2000):
    """Clip DP values to [0, cap] and return (clipped_pdf, n_clipped)."""
    n_over = int((pdf[col] > cap).sum())
    pdf[col] = pdf[col].clip(0, cap)
    return pdf, n_over


def _count_scale():
    """Linear scale for count axes — raw values, no log transform."""
    return alt.Scale()


def _count_axis(**kwargs):
    """SI-formatted axis labels for count data (1K, 10K, 100K, 1M, etc.)."""
    defaults = dict(format="~s")
    defaults.update(kwargs)
    return alt.Axis(**defaults)


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
    box = alt.Chart(pdf).mark_boxplot().encode(
        x=alt.X(f"{group_col}:N"),
        y=alt.Y(f"{value_col}:Q"),
        color=alt.Color(f"{color_col or group_col}:N") if color_col else alt.value("#1f77b4"),
    )
    chart = box.properties(title=title)
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

def plot_vc_distribution(df, output_dir: str, group_col: str = "set_number", facet_col: str = None):
    """Chart 1: Variant counts by VC classification, stacked bar per set."""
    if not _has_column(df, "FILTER"):
        return
    # Deduplicate group_by columns (handles variant-category-wise where group_col="FILTER")
    group_cols = list(dict.fromkeys([group_col, "FILTER"]))
    if facet_col and facet_col in df.columns and facet_col not in group_cols:
        group_cols.insert(0, facet_col)
    counts = _maybe_collect(
        df.group_by(group_cols).agg(pl.len().alias("count"))
        .sort(group_cols)
    )
    if group_col == "CHROM":
        counts = _sort_chromosomes(counts, group_col)
    chrom_order = counts[group_col].to_list() if group_col == "CHROM" else None
    counts = counts.to_pandas()
    counts[group_col] = counts[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()
    base = alt.Chart(counts)
    bars = base.mark_bar().encode(
        x=alt.X(f"{group_col}:N", title=group_title, sort=chrom_order),
        y=alt.Y("count:Q", title="Number of Variants", scale=_count_scale(), axis=_count_axis()),
        color=alt.Color(
            "FILTER:N", title="Variant Classification",
            scale=_color_scale("FILTER"),
            legend=alt.Legend(orient="right", title="Variant Classification",
                              labelFontSize=11, titleFontSize=12),
        ),
    )
    text = _make_bar_text(counts, x_field=group_col, count_col="count",
                          stack="zero", show_pct=True, group_col=group_col,
                          x_sort=chrom_order, base=base)
    chart = (bars + text).properties(title=f"Variant Classification Distribution by {group_title}")
    if facet_col:
        chart = chart.facet(facet=alt.Facet(f"{facet_col}:N"), columns=2).resolve_scale(x="independent")
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
    ).properties(title=f"Variant Tier Distribution by {group_title}")
    chart = _apply_faceting(chart, group_col)
    _save_chart(chart, "02_caller_overlap", output_dir)
    return chart


def plot_variant_type_distribution(df, output_dir: str, group_col: str = "set_number", facet_col: str = None):
    """Chart 9: Variant type distribution (SNV/INS/DEL/MNV), stacked bar per group."""
    if "variant_type" not in df.columns or group_col not in df.columns:
        return
    facet_groups = [facet_col] if facet_col and facet_col in df.columns else []
    total_per_group = df.group_by([*facet_groups, group_col]).agg(pl.len().alias("total"))
    pdf = df.group_by([*facet_groups, group_col, "variant_type"]).agg(pl.len().alias("count"))
    pdf = pdf.join(total_per_group, on=[*facet_groups, group_col]).with_columns(
        (pl.col("count") / pl.col("total") * 100).round(1).alias("pct")
    )
    pdf = _maybe_collect(pdf)
    if group_col == "CHROM":
        pdf = _sort_chromosomes(pdf, group_col)
    chrom_order = pdf[group_col].to_list() if group_col == "CHROM" else None
    pdf = pdf.to_pandas()
    pdf[group_col] = pdf[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()
    base = alt.Chart(pdf)
    bars = base.mark_bar().encode(
        x=alt.X(f"{group_col}:N", title=group_title, sort=chrom_order),
        y=alt.Y("count:Q", title="Number of Variants", scale=_count_scale(), axis=_count_axis()),
        color=alt.Color("variant_type:N", title="Variant Type", scale=_color_scale("variant_type")),
    )
    text = _make_bar_text(pdf, x_field=group_col, count_col="count",
                          stack="zero", show_pct=True, group_col=group_col,
                          x_sort=chrom_order, base=base)
    chart = (bars + text).properties(title=f"Variant Type Distribution by {group_title}")
    if facet_col:
        chart = chart.facet(facet=alt.Facet(f"{facet_col}:N"), columns=2).resolve_scale(x="independent")
    _save_chart(chart, "09_variant_type_distribution", output_dir)
    return chart


def plot_ti_tv_ratio(df, output_dir: str, group_col: str = "set_number", facet_col: str = None):
    """Chart 10: Ti/Tv ratio bar chart per group."""
    if "ti_tv" not in df.columns or group_col not in df.columns:
        return
    facet_groups = [facet_col] if facet_col and facet_col in df.columns else []
    ti = df.filter(pl.col("ti_tv") == True).group_by([*facet_groups, group_col]).agg(pl.len().alias("ti"))
    tv = df.filter(pl.col("ti_tv") == False).group_by([*facet_groups, group_col]).agg(pl.len().alias("tv"))
    ratio = ti.join(tv, on=[*facet_groups, group_col]).with_columns(
        (pl.col("ti") / pl.col("tv")).alias("ratio"))
    ratio = _maybe_collect(ratio)
    if group_col == "CHROM":
        ratio = _sort_chromosomes(ratio, group_col)
    chrom_order = ratio[group_col].to_list() if group_col == "CHROM" else None
    ratio = ratio.to_pandas()
    ratio[group_col] = ratio[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()
    base = alt.Chart(ratio)
    bars = base.mark_bar().encode(
        x=alt.X(f"{group_col}:N", title=group_title, sort=chrom_order),
        y=alt.Y("ratio:Q", title="Ti/Tv Ratio"),
    )
    text = base.mark_text(dy=-8).encode(
        x=alt.X(f"{group_col}:N"), y=alt.Y("ratio:Q"),
        text=alt.Text("ratio:Q", format=".2f"),
    )
    chart = (bars + text).properties(title=f"Ti/Tv Ratio by {group_title}")
    if facet_col:
        chart = chart.facet(facet=alt.Facet(f"{facet_col}:N"), columns=2).resolve_scale(x="independent")
    _save_chart(chart, "10_ti_tv_ratio", output_dir)
    return chart


def plot_cross_modality(df, output_dir: str, group_col: str = "set_number", facet_col: str = None):
    """Chart 11: Modality evidence & cross-modality rescue analysis with percentages.

    If modality_evidence_caller column is available, uses 4-category classification
    (cross_modality, dna_confident, rna_rescued, low_confidence). Falls back to
    RESCUED column for backward compatibility with older parquet files.
    """
    if group_col not in df.columns:
        return

    # Determine evidence column: prefer modality_evidence_caller, fall back to RESCUED
    evidence_col = "modality_evidence_caller"
    if not _has_column(df, evidence_col):
        has_cross_modality = _has_column(df, "CROSS_MODALITY")
        has_rescued = _has_column(df, "RESCUED")
        if not has_cross_modality and not has_rescued:
            return
        # Backward compat: use CROSS_MODALITY and RESCUED as before
        cols_to_plot = [c for c in ["CROSS_MODALITY", "RESCUED"] if _has_column(df, c)]
        title_prefix_map = {"CROSS_MODALITY": "Cross-Modality", "RESCUED": "Rescued Variants"}
    else:
        cols_to_plot = [evidence_col]
        title_prefix_map = {evidence_col: "Modality Evidence"}

    facet_groups = [facet_col] if facet_col and facet_col in df.columns else []
    total_per_group = df.group_by([*facet_groups, group_col]).agg(pl.len().alias("total"))

    subcharts = []
    for col in cols_to_plot:
        pdf = df.group_by([*facet_groups, group_col, col]).agg(
            pl.len().alias("count")
        ).join(total_per_group, on=[*facet_groups, group_col]).with_columns(
            (pl.col("count") / pl.col("total") * 100).round(1).alias("pct")
        )
        pdf = _maybe_collect(pdf)
        if group_col == "CHROM":
            pdf = _sort_chromosomes(pdf, group_col)
        chrom_order = pdf[group_col].to_list() if group_col == "CHROM" else None
        pdf = pdf.to_pandas()
        pdf[group_col] = pdf[group_col].astype(str)
        group_title = group_col.replace("_", " ").title()

        base = alt.Chart(pdf)
        bars = base.mark_bar().encode(
            x=alt.X(f"{group_col}:N", title=group_title, sort=chrom_order),
            y=alt.Y("count:Q", title="Count", scale=_count_scale(), axis=_count_axis()),
            color=alt.Color(f"{col}:N", scale=_color_scale(col) if col in _COLOR_REGISTRY else alt.Scale()),
        )
        text = _make_bar_text(pdf, x_field=group_col, count_col="count",
                              stack="zero", show_pct=True, group_col=group_col,
                              x_sort=chrom_order, base=base)
        c = (bars + text).properties(
            title=title_prefix_map.get(col, col.replace("_", " ").title())
        )
        if facet_col:
            c = c.facet(facet=alt.Facet(f"{facet_col}:N"), columns=2).resolve_scale(x="independent")
        subcharts.append(c)

    chart_title = "Modality Evidence Analysis (%)" if evidence_col in cols_to_plot else "Cross-Modality and Rescue Analysis (%)"
    chart = alt.hconcat(*subcharts).properties(title=chart_title) if len(subcharts) > 1 else subcharts[0]
    _save_chart(chart, "11_modality_evidence", output_dir)
    return chart


def plot_filter_distribution(df, output_dir: str, group_col: str = "set_number"):
    """FILTER value distribution per group (bar chart with % marks)."""
    if "FILTER" not in df.columns or group_col not in df.columns:
        return
    total_per_group = df.group_by(group_col).agg(pl.len().alias("total"))
    # Deduplicate group_by columns (defensive — handles group_col="FILTER")
    group_cols = list(dict.fromkeys([group_col, "FILTER"]))
    counts = df.group_by(group_cols).agg(pl.len().alias("count"))
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
        y=alt.Y("count:Q", title="Number of Variants", scale=_count_scale(), axis=_count_axis()),
        color=alt.Color("FILTER:N", title="FILTER"),
    )
    text = _make_bar_text(counts, x_field=group_col, count_col="count",
                          stack="zero", show_pct=True, group_col=group_col,
                          x_sort=chrom_order)
    chart = (bars + text).properties(title=f"FILTER Distribution by {group_title}")
    _save_chart(chart, "24_filter_distribution", output_dir)
    return chart


def plot_chromosome_density(df, output_dir: str):
    """Variant count per chromosome with natural sort order (chr1..chrX, chrY, chrM)."""
    if "CHROM" not in df.columns:
        return
    counts = df.group_by("CHROM").agg(pl.len().alias("count"))
    counts = _maybe_collect(counts)
    counts = _sort_chromosomes(counts, "CHROM")

    chrom_order = counts["CHROM"].to_list()
    pdf = counts.to_pandas()
    bars = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("CHROM:N", title="Chromosome", sort=chrom_order),
        y=alt.Y("count:Q", title="Number of Variants", scale=_count_scale(), axis=_count_axis()),
        tooltip=["CHROM", "count"],
    )
    text = _make_bar_text(pdf, x_field="CHROM", count_col="count", x_sort=chrom_order)
    chart = (bars + text).properties(title="Variant Density per Chromosome")
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
    bars = alt.Chart(counts).mark_bar().encode(
        x=alt.X(f"{group_col}:N", title=group_title, sort=chrom_order),
        y=alt.Y("count:Q", title="Number of Variants", scale=_count_scale(), axis=_count_axis()),
        color=alt.Color("REDI_EVIDENCE:N", title="REDIportal Evidence Level"),
    )
    text = _make_bar_text(counts, x_field=group_col, count_col="count",
                          stack="zero", show_pct=True, group_col=group_col,
                          x_sort=chrom_order)
    chart = (bars + text).properties(title=f"REDIportal RNA Editing Evidence by {group_title}")
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
        y=alt.Y("count:Q", title="Number of Variants", scale=_count_scale(), axis=_count_axis()),
        color=alt.Color("variant_type:N", title="Variant Type", scale=_color_scale("variant_type")),
        column=col_enc,
    )
    text = _make_bar_text(pdf, x_field="variant_type", count_col="count")
    text = text.encode(column=col_enc)
    chart = (bars + text).properties(title="Variant Type Distribution per Caller Tier")
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
        enc["color"] = alt.Color("caller:N", scale=_color_scale("caller"))
        box = alt.Chart(pdf).mark_boxplot().encode(**enc)
        chart = (box + ref_rules).properties(
            title="VAF Distribution per Caller",
            width=alt.Step(60))
    else:
        # Faceted: boxplot with column encoding (NOT .facet() — boxplot is composite mark)
        facet_col_name = color_col
        base_enc = {"x": alt.X("caller:N", title="Caller", axis=alt.Axis(labelAngle=-45)),
                    "y": y_scale,
                    "color": alt.Color("caller:N", scale=_color_scale("caller")),
                    "column": alt.Column(f"{facet_col_name}:N")}
        box = alt.Chart(pdf).mark_boxplot().encode(**base_enc)
        chart = box.properties(
            title="VAF Distribution per Caller",
            width=alt.Step(60)
        )
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
        actual_color = alt.Color("FILTER:N", scale=_color_scale("FILTER"))
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
    n_over = int((pdf["DNA_DP_mean"] > 2000).sum() + (pdf["RNA_DP_mean"] > 2000).sum())
    pdf["DNA_DP_mean"] = pdf["DNA_DP_mean"].clip(0, 2000)
    pdf["RNA_DP_mean"] = pdf["RNA_DP_mean"].clip(0, 2000)
    group_title = group_col.replace("_", " ").title()
    subtitle = f"{n_over} values > 2000 clipped" if n_over else ""
    enc = {"x": alt.X("DNA_DP_mean:Q", title="DNA Mean Depth (capped at 2000)", scale=alt.Scale(domain=[0, 2000])),
           "y": alt.Y("RNA_DP_mean:Q", title="RNA Mean Depth (capped at 2000)", scale=alt.Scale(domain=[0, 2000]))}
    if group_col in pdf.columns:
        enc["color"] = alt.Color(f"{group_col}:N", title=group_title)
    chart = alt.Chart(pdf).mark_circle(opacity=0.4, size=20).encode(**enc).properties(
        title=alt.Title("DNA vs RNA Mean Depth (capped at 2000)", subtitle=subtitle))
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
                "color": alt.Color("caller:N", scale=_color_scale("caller")),
                "column": col_enc,
                "y": y_scale}
    chart = alt.Chart(pdf).mark_boxplot().encode(**base_enc).properties(
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
    n_over = int((pdf["DP"] > 2000).sum())
    pdf["DP"] = pdf["DP"].clip(0, 2000)
    pdf["caller"] = pdf["caller"].str.replace("_DP", "")
    col_enc = alt.Column("caller_tier:N", title="Caller Tier")
    if facet_col and facet_col in pdf.columns:
        pdf[facet_col] = pdf[facet_col].astype(str)
        col_enc = alt.Column(f"{facet_col}:N", title=facet_col.replace("_", " ").title())
    subtitle = f"{n_over} values > 2000 clipped" if n_over else ""
    base_enc = {"x": alt.X("caller:N", title="Caller", axis=alt.Axis(labelAngle=-45)),
                "color": alt.Color("caller:N", scale=_color_scale("caller")), "column": col_enc,
                "y": alt.Y("DP:Q", title="Read Depth (capped at 2000)", scale=alt.Scale(domain=[0, 2000]))}
    chart = alt.Chart(pdf).mark_boxplot().encode(**base_enc).properties(
        title=alt.Title("DP Distribution per Caller × Caller Tier (capped at 2000)", subtitle=subtitle))
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
        n_over = int((pdf[x_col] > 2000).sum() + (pdf[y_col] > 2000).sum())
        pdf[x_col] = pdf[x_col].clip(0, 2000)
        pdf[y_col] = pdf[y_col].clip(0, 2000)
        group_title = group_col.replace("_", " ").title()
        subtitle = f"{n_over} values > 2000 clipped" if n_over else ""
        enc = {"x": alt.X(f"{x_col}:Q", title=f"DNA Mean {label} (capped at 2000)", scale=alt.Scale(domain=[0, 2000])),
               "y": alt.Y(f"{y_col}:Q", title=f"RNA Mean {label} (capped at 2000)", scale=alt.Scale(domain=[0, 2000]))}
        if group_col in pdf.columns:
            pdf[group_col] = pdf[group_col].astype(str)
            enc["color"] = alt.Color(f"{group_col}:N", title=group_title)
        c = alt.Chart(pdf).mark_circle(opacity=0.4, size=20).encode(**enc).properties(
            title=alt.Title(f"DNA vs RNA Mean {label} (capped at 2000)", subtitle=subtitle))
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
        color=alt.Color("final_tier:N", scale=_color_scale("final_tier")),
    ).properties(title="DNA VAF Distribution per Tier")
    _save_chart(chart, "25_per_tier_vaf", output_dir)
    return chart


def plot_dna_vs_rna_per_caller(df, output_dir: str, color_col: str = None,
                              metrics=("VAF", "DP")):
    """DNA vs RNA per-caller VAF/DP scatter at shared positions.

    Layout: vconcat(hconcat(VAF_row), hconcat(DP_row)) — 2 rows x 3 columns.
    DP subcharts are clipped to [0, 2000]. All panels share a color scale
    when color_col is provided. Sampled to 5K per pair.
    """
    pairs = [("DNA_mutect2", "RNA_mutect2"), ("DNA_deepsomatic", "RNA_deepsomatic"),
             ("DNA_strelka", "RNA_strelka")]

    # Build subcharts organized by row: {metric: [charts]}
    metric_rows = {}
    # Shared color domain for all panels (built from first pair)
    shared_color_domain = None
    shared_color_range = None

    for metric in metrics:
        subcharts = []
        for dna_caller, rna_caller in pairs:
            dna_col = f"{dna_caller}_{metric}"
            rna_col = f"{rna_caller}_{metric}"
            if dna_col not in df.columns or rna_col not in df.columns:
                continue
            cols = [dna_col, rna_col]
            active_color = None
            if color_col and color_col in df.columns:
                cols.append(color_col)
                active_color = color_col
            elif "FILTER" in df.columns:
                cols.append("FILTER")
                active_color = "FILTER"
            pdf = df.select(cols).drop_nulls(subset=[dna_col, rna_col])
            pdf = _sample_if_large(pdf, max_rows=5000).to_pandas()

            # DP: clip to [0, 2000]
            if metric == "DP":
                n_over = int((pdf[dna_col] > 2000).sum() + (pdf[rna_col] > 2000).sum())
                pdf[dna_col] = pdf[dna_col].clip(0, 2000)
                pdf[rna_col] = pdf[rna_col].clip(0, 2000)

            caller_label = dna_caller.replace("DNA_", "")
            if active_color and active_color in pdf.columns:
                if active_color == "FILTER":
                    color_enc = alt.Color("FILTER:N", scale=_color_scale("FILTER"),
                                          legend=alt.Legend(title="Classification"))
                else:
                    color_enc = alt.Color(f"{active_color}:N",
                                          legend=alt.Legend(title=active_color.replace("_", " ").title()))
            else:
                color_enc = alt.value("#1f77b4")

            x_title = f"{caller_label} DNA {metric}"
            y_title = f"{caller_label} RNA {metric}"
            if metric == "DP":
                x_title += " (capped 2000)"
                y_title += " (capped 2000)"

            c = alt.Chart(pdf).mark_circle(opacity=0.4, size=20).encode(
                x=alt.X(f"{dna_col}:Q", title=x_title),
                y=alt.Y(f"{rna_col}:Q", title=y_title),
                color=color_enc,
            ).properties(title=f"{caller_label}", width=150, height=150)
            subcharts.append(c)
        if subcharts:
            metric_rows[metric] = alt.hconcat(*subcharts)

    if not metric_rows:
        return

    # vconcat rows: VAF first, DP second
    row_charts = [metric_rows[m] for m in metrics if m in metric_rows]
    chart = alt.vconcat(*row_charts).properties(
        title=alt.Title("DNA vs RNA Per-Caller VAF and DP", subtitle="VAF row (top), DP row (bottom, capped 2000)"))
    _save_chart(chart, "26_dna_vs_rna_per_caller", output_dir)
    return chart


# ═══════════════════════════════════════════════════════════════════════════
# Count-based charts — use _count_rows(), no data columns loaded.
# Memory: < 1 KB per chart (just integer counts).
# ═══════════════════════════════════════════════════════════════════════════

def plot_cosmic_gnomad_annotation(df, output_dir: str, group_col: str = "set_number", facet_col: str = None):
    """Chart 7: COSMIC/gnomAD annotation coverage as side-by-side pies, per group."""
    has_cosmic_col = "COSMIC_ID" in df.columns
    has_gnomad_col = "GNOMAD_AF" in df.columns
    if not has_cosmic_col and not has_gnomad_col:
        return

    has_groups = group_col in df.columns
    facet_groups = [facet_col] if facet_col and facet_col in df.columns else []

    if has_groups:
        # Per-group annotation bars (single group_by scan, not per-group _count_rows)
        agg_exprs = [pl.len().alias("n_total")]
        if has_cosmic_col:
            agg_exprs.append(pl.col("COSMIC_ID").is_not_null().sum().alias("n_cosmic"))
        if has_gnomad_col:
            agg_exprs.append(pl.col("GNOMAD_AF").is_not_null().sum().alias("n_gnomad"))
        counts_df = _maybe_collect(df.group_by([*facet_groups, group_col]).agg(agg_exprs))
        if counts_df.is_empty():
            return
        rows = []
        for rd in counts_df.iter_rows(named=True):
            g_n = rd["n_total"]
            r = {"group": str(rd[group_col]), "n_total": g_n}
            if facet_col:
                r[facet_col] = str(rd[facet_col])
            if has_cosmic_col:
                r["cosmic_pct"] = round(rd["n_cosmic"] / g_n * 100, 2) if g_n > 0 else 0
            if has_gnomad_col:
                r["gnomad_pct"] = round(rd["n_gnomad"] / g_n * 100, 2) if g_n > 0 else 0
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
            if facet_col:
                c = c.facet(facet=alt.Facet(f"{facet_col}:N"), columns=2).resolve_scale(x="independent")
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
    base = alt.Chart(pdf)
    rect = base.mark_rect().encode(
        x=alt.X("caller_1:N", title=None),
        y=alt.Y("caller_2:N", title=None),
        color=alt.Color("pct:Q", title="% Both Valid GT", scale=alt.Scale(scheme="blues")),
        tooltip=["caller_1", "caller_2", "pct"],
    )
    text = _add_heatmap_text(base, alt.X("caller_1:N"), alt.Y("caller_2:N"), "pct", pdf, fontSize=8, fmt=".1f")
    chart = (rect + text).properties(title="Caller GT Availability Matrix (%)")
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
        # Fill null group values to prevent "undefined" facet titles
        pdf = pdf.with_columns(pl.col(group_col).cast(pl.Utf8).fill_null("Unknown"))
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
        base = alt.Chart(counts_df)
        bars = base.mark_bar().encode(
            x=alt.X("agreement_level:N", title="Number of Callers Agreeing"),
            y=alt.Y("count:Q", title="Number of Variants", scale=_count_scale(), axis=_count_axis()),
            color=alt.Color("agreement_level:N"),
        )
        text = _make_bar_text(counts_df, x_field="agreement_level", count_col="count", base=base)
        chart = (bars + text).properties(title=f"GT Concordance by {group_title}")
        chart = _apply_faceting(chart, group_col)
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
        bars = alt.Chart(counts_df).mark_bar().encode(
            x=alt.X("agreement_level:N", title="Number of Callers Agreeing on GT"),
            y=alt.Y("count:Q", title="Number of Variants", scale=_count_scale(), axis=_count_axis()),
        )
        text = _make_bar_text(counts_df, x_field="agreement_level", count_col="count")
        chart = (bars + text).properties(title=f"GT Concordance Among 4 Callers (n={n_total} with ≥2 agreement)")
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
    bars = alt.Chart(result).mark_bar().encode(
        x=alt.X("agreement:N", title="Number of Callers Agreeing"),
        y=alt.Y("count:Q", title="Count", scale=_count_scale(), axis=_count_axis()),
        color=alt.Color("agreement:N"),
        column=col_enc,
    )
    text = _make_bar_text(result, x_field="agreement", count_col="count")
    text = text.encode(column=col_enc)
    chart = (bars + text).properties(title="GT Concordance per Caller Tier")
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
        "x": alt.X("sample_id:N", title="Sample ID", sort=None,
                    axis=alt.Axis(labelAngle=-45, labelLimit=100)),
        "y": alt.Y("total_variants:Q", title="Total Variants per Sample", scale=_count_scale(), axis=_count_axis()),
        "tooltip": ["sample_id", "total_variants"],
    }
    base = alt.Chart(pdf)
    bars = base.mark_bar().encode(**enc)
    text = _make_bar_text(pdf, x_field="sample_id", count_col="total_variants",
                          base=base, overlap="hide")
    chart = (bars + text).properties(
        title=f"Per-Sample Variant Counts (n={n_samples})",
        height=300,
    )
    if has_set:
        chart = chart.facet(
            facet=alt.Facet("set_number:N", title="Set"),
            columns=2,
        ).resolve_scale(x="independent", y="shared")
    _save_chart(chart, "08_per_sample_distribution", output_dir)
    return chart


def plot_validation_heatmap(report, output_dir: str):
    """Chart 13: Rescue VCF validation heatmap (mismatch % per sample × metric).

    Transposed layout: samples on x-axis, metrics on y-axis.
    Faceted by set_number when available.
    """
    if report is None or (hasattr(report, 'is_empty') and report.is_empty()):
        return
    pdf = report.to_pandas()
    has_set = "set_number" in pdf.columns
    # Build a sample→set lookup before pivoting (pivot drops extra columns)
    set_lookup = pdf[["sample_id", "set_number"]].drop_duplicates() if has_set else None
    pivot = pdf.pivot(
        index="sample_id", columns="metric", values="mismatch_pct"
    ).reset_index().melt(id_vars="sample_id", var_name="metric", value_name="mismatch_pct")
    # Re-join set_number after melt
    if set_lookup is not None:
        pivot = pivot.merge(set_lookup, on="sample_id", how="left")
    base = alt.Chart(pivot)
    rect = base.mark_rect().encode(
        x=alt.X("sample_id:N", title="Sample", axis=alt.Axis(labelAngle=-45)),
        y=alt.Y("metric:N", title="Validation Metric"),
        color=alt.Color("mismatch_pct:Q", title="Mismatch %",
                        scale=alt.Scale(scheme="redyellowgreen", reverse=True)),
    )
    text = _add_heatmap_text(base, alt.X("sample_id:N"), alt.Y("metric:N"),
                             "mismatch_pct", pivot, fontSize=7, fmt=".1f")
    chart = (rect + text).properties(title="Rescue VCF Validation: Mismatch % by Sample × Metric")
    if has_set:
        chart = chart.facet(
            facet=alt.Facet("set_number:N", title="Set"), columns=2,
        ).resolve_scale(x="independent")
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
        "color": alt.Color("bam_type:N", title="BAM Type", scale=_color_scale("bam_type")),
        "xOffset": "bam_type:N",
    }
    if has_set:
        pdf["set_number"] = pdf["set_number"].astype(str)
    chart = alt.Chart(pdf).mark_bar().encode(**enc).properties(
        title="Per-Sample BAM Read Counts — DN/DT/RT",
        width=350)
    if has_set:
        chart = chart.facet(
            facet=alt.Facet("set_number:N", title="Set"),
            columns=2,
        ).resolve_scale(x="independent", y="shared")
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
        "x": alt.X("sample_id:N", title="Sample ID", sort=None,
                    axis=alt.Axis(labelAngle=-45)),
        "y": alt.Y("n_variants:Q", title="Variants", scale=_count_scale(), axis=_count_axis()),
        "color": alt.Color("final_tier:N", title="Tier", scale=_color_scale("final_tier")),
    }
    base = alt.Chart(pdf)
    bars = base.mark_bar().encode(**enc)
    text = _make_bar_text(pdf, x_field="sample_id", count_col="n_variants",
                          stack="zero", show_pct=True, group_col="sample_id",
                          base=base, overlap="hide")
    chart = (bars + text).properties(
        title="Per-Sample Per-Tier Variant Distribution")
    if has_set:
        chart = chart.facet(
            facet=alt.Facet("set_number:N", title="Set"),
            columns=2,
        ).resolve_scale(x="independent")
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

    enc = {"x": alt.X("metric:N", title="BAM Metric"),
           "y": alt.Y("depth:Q", title="Depth at Variant Position"),
           "color": alt.Color("metric:N", title="BAM Metric", scale=_color_scale("metric"), legend=None)}
    if color_col and color_col in sampled.columns:
        enc["column"] = alt.Column(f"{color_col}:N")
    chart = alt.Chart(sampled).mark_boxplot().encode(**enc).properties(
        title="BAM Pileup Depth Distribution at Variant Positions (depth capped at 2000, sampled)")
    _save_chart(chart, "21_bam_coverage_boxplot", output_dir)
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
        color=alt.Color("final_tier:N", scale=_color_scale("final_tier")),
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
    # Include both overall and per-classification rows
    pdf = sweep_df
    if isinstance(pdf, pl.LazyFrame):
        pdf = pdf.collect()
    if pdf.is_empty():
        return
    if "classification" in pdf.columns:
        pdf = pdf.with_columns(pl.col("classification").fill_null("Overall"))
    else:
        pdf = pdf.with_columns(pl.lit("Overall").alias("classification"))
    pdf = pdf.to_pandas()
    chart = alt.Chart(pdf).mark_line(point=True).encode(
        x=alt.X("threshold:Q", title="VAF Threshold"),
        y=alt.Y("pct_retained:Q", title="% Variants Retained"),
        color=alt.Color("caller:N", title="Caller", scale=_color_scale("caller")),
    ).properties(
        title="VAF Threshold Sweep: Retention % vs Threshold per Caller",
        width=250, height=200,
    ).facet(
        facet=alt.Facet("classification:N", title="Classification"),
        columns=3,
    )
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
    # Include both overall and per-classification rows
    pdf = sweep_df
    if isinstance(pdf, pl.LazyFrame):
        pdf = pdf.collect()
    if pdf.is_empty():
        return
    if "classification" in pdf.columns:
        pdf = pdf.with_columns(pl.col("classification").fill_null("Overall"))
    else:
        pdf = pdf.with_columns(pl.lit("Overall").alias("classification"))
    pdf = pdf.to_pandas()
    # Create combined facet label for grid layout (metric × classification)
    pdf["panel"] = pdf["metric"] + " — " + pdf["classification"]
    chart = alt.Chart(pdf).mark_line(point=True).encode(
        x=alt.X("threshold:Q", title="Depth Threshold"),
        y=alt.Y("pct_retained:Q", title="% Variants Retained"),
        color=alt.Color("caller:N", title="Caller", scale=_color_scale("caller")),
    ).properties(
        title="DP Threshold Sweep: Retention % vs Threshold per Caller",
        width=220, height=180,
    ).facet(
        facet=alt.Facet("panel:N", title="DP Metric — Classification"),
        columns=3,
    )
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
        # Use column encoding (NOT .facet()) — boxplot is composite mark
        enc["column"] = alt.Column(f"{color_col}:N")
    chart = alt.Chart(pdf).mark_boxplot().encode(**enc).properties(
        title="Caller Concordance vs VAF", width=200)
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
    base = alt.Chart(pdf)
    rect = base.mark_rect().encode(
        x=alt.X("filter_flag:N", title="Filter Flag"),
        y=alt.Y("classification:N", title="Classification"),
        color=alt.Color("pct_flagged:Q", title="% Flagged",
                        scale=alt.Scale(scheme="redyellowgreen", reverse=True)),
        tooltip=["classification", "filter_flag", "pct_flagged"],
    )
    text = _add_heatmap_text(base, alt.X("filter_flag:N"), alt.Y("classification:N"), "pct_flagged", pdf, fontSize=7, fmt=".1f")
    chart = (rect + text).properties(title="Filter Effectiveness: % Flagged by Classification × Filter")
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
    n_over = int((pdf["DP"] > 2000).sum())
    pdf["DP"] = pdf["DP"].clip(0, 2000)
    pdf["caller"] = pdf["caller"].str.replace("_DP", "")
    enc = {"x": alt.X("caller:N", title="Caller", axis=alt.Axis(labelAngle=-45)),
           "y": alt.Y("DP:Q", title="Read Depth (capped at 2000)", scale=alt.Scale(domain=[0, 2000]))}
    if color_col and color_col in pdf.columns:
        # Use column encoding (NOT .facet()) — boxplot is composite mark
        enc["color"] = alt.Color(f"{color_col}:N")
        enc["column"] = alt.Column(f"{color_col}:N")
        chart = alt.Chart(pdf).mark_boxplot().encode(**enc).properties(
            title="DP Distribution per Caller (capped at 2000)")
    else:
        enc["color"] = alt.Color("caller:N", scale=_color_scale("caller"))
        box = alt.Chart(pdf).mark_boxplot().encode(**enc)
        chart = box.properties(title="DP Distribution per Caller (capped at 2000)")
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
    ).properties(title=f"Caller Support Distribution by {group_title}")
    chart = _apply_faceting(chart, group_col)
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
    base = alt.Chart(pdf)
    rect = base.mark_rect().encode(
        x=alt.X("tier:N", title="Tier (CxDy)"),
        y=alt.Y("caller:N", title="Caller"),
        color=alt.Color("pct_detected:Q", title="% Detected", scale=alt.Scale(scheme="blues")),
        tooltip=["caller", "tier", "pct_detected"],
    )
    text = _add_heatmap_text(base, alt.X("tier:N"), alt.Y("caller:N"), "pct_detected", pdf, fontSize=7, fmt=".1f")
    chart = (rect + text)
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
    # Collapse partition dimension (train/val/test) if present — prevents
    # 3 overlapping text labels per cell.
    agg_cols = ["classification", "vaf_bin", "dp_bin"]
    if any(c not in agg_cols for c in cross_tab_df.columns if c != "count"):
        cross_tab_df = cross_tab_df.group_by(agg_cols).agg(pl.col("count").sum())
    pdf = cross_tab_df.to_pandas()
    # Order VAF and DP bins sensibly
    vaf_order = ["<0.01", "0.01-0.05", "0.05-0.10", "0.10-0.25", "0.25-0.50", "0.50-1.0"]
    dp_order = ["<10", "10-50", "50-100", "100-200", "200-500", "500+"]
    base = alt.Chart(pdf)
    rect = base.mark_rect().encode(
        x=alt.X("vaf_bin:N", title="VAF Bin", sort=vaf_order),
        y=alt.Y("dp_bin:N", title="DP Bin", sort=dp_order),
        color=alt.Color("count:Q", title="Count", scale=alt.Scale(scheme="blues", type="log")),
        tooltip=["classification", "vaf_bin", "dp_bin", "count"],
    )
    median_count = max(float(pdf["count"].median()), 1)
    text = base.mark_text(baseline="middle", fontSize=8).encode(
        x=alt.X("vaf_bin:N", sort=vaf_order),
        y=alt.Y("dp_bin:N", sort=dp_order),
        text=alt.Text("count:Q", format=",d"),
        color=alt.condition(
            f"datum.count > {median_count}",
            alt.value("white"),
            alt.value("black"),
        ),
    )
    chart = (rect + text).properties(
        title="FILTER x VAF x DP Cross-Tabulation", width=200, height=180
    ).facet(
        facet=alt.Facet("classification:N", title="Classification"),
        columns=3,
    )
    _save_chart(chart, "43_filter_vaf_dp_heatmap", output_dir)
    return chart


def plot_low_vaf_rna_support(low_vaf_df, output_dir: str):
    """Chart 44: Low VAF variants with strong RNA support — stacked bar by FILTER."""
    if low_vaf_df is None or (hasattr(low_vaf_df, 'is_empty') and low_vaf_df.is_empty()):
        return
    if "FILTER" not in low_vaf_df.columns or "n_variants" not in low_vaf_df.columns:
        return
    pdf = low_vaf_df.to_pandas()
    bars = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("FILTER:N", title="Classification"),
        y=alt.Y("n_variants:Q", title="Number of Variants", scale=_count_scale(), axis=_count_axis()),
        color=alt.Color("FILTER:N", title="Classification",
                        scale=_color_scale("FILTER")),
        tooltip=["FILTER", "n_variants", "mean_vaf"],
    )
    text = _make_bar_text(pdf, x_field="FILTER", count_col="n_variants")
    chart = (bars + text).properties(title="Low VAF (< 0.05) Variants with RNA Support (N_RNA_CALLERS >= 2)")
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
    base = alt.Chart(pdf)
    rect = base.mark_rect().encode(
        x=alt.X("N_SUPPORT_CALLERS:N", title="Number of Supporting Callers"),
        y=alt.Y("FILTER:N", title="Classification"),
        color=alt.Color("count:Q", title="Count", scale=alt.Scale(scheme="orangered", type="log")),
        tooltip=["FILTER", "N_SUPPORT_CALLERS", "count"],
    )
    median_count = max(float(pdf["count"].median()), 1)
    text = base.mark_text(baseline="middle", fontSize=8).encode(
        x=alt.X("N_SUPPORT_CALLERS:N"),
        y=alt.Y("FILTER:N"),
        text=alt.Text("count:Q", format=".2s"),
        color=alt.condition(
            f"datum.count > {median_count}",
            alt.value("white"),
            alt.value("black"),
        ),
    )
    chart = (rect + text).properties(title="Classification × Caller Support Cross-Tabulation",
                                     width=210, height=210)
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
    pdf = modality_df.to_pandas()
    # Filter out zero-count categories to avoid phantom legend entries
    pdf = pdf[pdf["n_variants"] > 0].copy()
    if pdf.empty:
        return
    total = pdf["n_variants"].sum()
    pdf["pct"] = (pdf["n_variants"] / total * 100).round(1)
    pdf["label"] = pdf.apply(lambda r: f"{r['n_variants']:,}\n({r['pct']}%)", axis=1)
    # Build color scale from only present categories
    present = pdf["somatic_modality"].tolist()
    domain = [d for d in SOMATIC_MODALITY_DOMAIN if d in present]
    colors = [SOMATIC_MODALITY_COLORS[SOMATIC_MODALITY_DOMAIN.index(d)] for d in domain]
    arc = alt.Chart(pdf).mark_arc(innerRadius=40).encode(
        theta=alt.Theta("n_variants:Q", stack=True),
        color=alt.Color("somatic_modality:N", title="Somatic Modality",
                        scale=alt.Scale(domain=domain, range=colors)),
        tooltip=["somatic_modality", "n_variants", "pct"],
    )
    text = alt.Chart(pdf).mark_text(size=10, radiusOffset=20, radius=90).encode(
        theta=alt.Theta("n_variants:Q", stack=True),
        text="label:N",
    )
    chart = (arc + text).properties(title="Somatic Modality Sub-Classification")
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
        y=alt.Y("n_variants:Q", title="Number of Variants", scale=_count_scale(), axis=_count_axis()),
        color=alt.Color("somatic_modality:N", title="Modality", scale=_color_scale("somatic_modality")),
        tooltip=list(pdf.columns),
    )
    text = _make_bar_text(pdf, x_field="somatic_modality", count_col="n_variants")
    bars_with_text = (bars + text).properties(title="Somatic Modality Sub-Classification: Variant Counts")
    # VAF overlay if available
    subcharts = [bars_with_text]
    if "mean_dna_vaf" in pdf.columns:
        vaf_chart = alt.Chart(pdf).mark_bar().encode(
            x=alt.X("somatic_modality:N", title="Somatic Modality"),
            y=alt.Y("mean_dna_vaf:Q", title="Mean DNA VAF"),
            color=alt.Color("somatic_modality:N", scale=_color_scale("somatic_modality")),
        ).properties(title="Mean DNA VAF by Modality")
        subcharts.append(vaf_chart)
    chart = alt.hconcat(*subcharts) if len(subcharts) > 1 else bars_with_text
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
                                    "cross-tabulation", "low vaf", "fp cross", "partition", "modality",
                                    "rescue", "rescued"]):
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
    sampled, n_over = _clip_dp(sampled, "depth", cap=2000)
    subtitle = f"{n_over} values > 2000 clipped" if n_over else ""

    enc = {"x": alt.X("metric:N", title="BAM Metric", axis=alt.Axis(labelAngle=-45)),
           "y": alt.Y("depth:Q", title="Depth at Variant Position (capped 2000)", scale=alt.Scale(domain=[0, 2000]))}
    if color_col and color_col in sampled.columns:
        # Use column encoding (NOT .facet()) — boxplot is composite mark
        enc["color"] = alt.Color(f"{color_col}:N")
        enc["column"] = alt.Column(f"{color_col}:N")
        chart = alt.Chart(sampled).mark_boxplot().encode(**enc).properties(
            title=alt.Title("BAM Pileup DP Distribution per BAM Type", subtitle=subtitle))
    else:
        enc["color"] = alt.Color("metric:N", scale=_color_scale("metric"), legend=None)
        box = alt.Chart(sampled).mark_boxplot().encode(**enc)
        chart = box.properties(title=alt.Title("BAM Pileup DP Distribution per BAM Type", subtitle=subtitle))
    _save_chart(chart, "36_bam_dp_distribution", output_dir)
    return chart


def plot_per_tier_dp_boxplot(df, output_dir: str):
    """Chart 37: Per-tier DNA DP boxplot (mirrors plot_per_tier_vaf_boxplot)."""
    if "DNA_DP_mean" not in df.columns or "final_tier" not in df.columns:
        return
    pdf = df.select(["final_tier", "DNA_DP_mean"]).drop_nulls(subset=["DNA_DP_mean"])
    pdf = _sample_if_large(pdf, max_rows=50000).to_pandas()
    n_over = int((pdf["DNA_DP_mean"] > 2000).sum())
    pdf["DNA_DP_mean"] = pdf["DNA_DP_mean"].clip(0, 2000)
    chart = alt.Chart(pdf).mark_boxplot().encode(
        x=alt.X("final_tier:N", title="Tier"),
        y=alt.Y("DNA_DP_mean:Q", title="DNA Mean DP (capped at 2000)", scale=alt.Scale(domain=[0, 2000])),
        color=alt.Color("final_tier:N", scale=_color_scale("final_tier")),
    ).properties(title=alt.Title("DNA DP Distribution per Tier (capped at 2000)",
                                 subtitle=f"{n_over} values > 2000 clipped" if n_over else ""))
    _save_chart(chart, "37_per_tier_dp", output_dir)
    return chart


# ═══════════════════════════════════════════════════════════════════════════
# Rescue Analytics Charts (Section 10)
# ═══════════════════════════════════════════════════════════════════════════


def plot_rescue_breakdown(breakdown_df, output_dir: str,
                          evidence_col: str = "modality_evidence_caller"):
    """Chart 48: Per-set rescued vs non-rescued counts with % labels."""
    if breakdown_df is None or (hasattr(breakdown_df, 'is_empty') and breakdown_df.is_empty()):
        return
    # Fall back to RESCUED if evidence_col not in columns
    if evidence_col not in breakdown_df.columns:
        evidence_col = "RESCUED"
    pdf = breakdown_df.to_pandas()
    non_evidence = [c for c in pdf.columns if c not in (evidence_col, "count", "total", "pct")]
    group_col = non_evidence[0] if non_evidence else None
    x_enc = alt.X(f"{group_col}:N", title=group_col.replace("_", " ").title()) if group_col else alt.X(f"{evidence_col}:N")
    x_field = group_col if group_col else evidence_col
    evidence_title = evidence_col.replace("_", " ").title()
    bars = alt.Chart(pdf).mark_bar().encode(
        x=x_enc,
        y=alt.Y("count:Q", title="Variant Count", scale=_count_scale(), axis=_count_axis(), stack="zero"),
        color=alt.Color(f"{evidence_col}:N", title=evidence_title,
                        scale=_color_scale(evidence_col) if evidence_col in _COLOR_REGISTRY else alt.Scale()),
        tooltip=list(pdf.columns),
    )
    text = _make_bar_text(pdf, x_field=x_field, count_col="count",
                          stack="zero", show_pct=True, group_col=x_field)
    chart = (bars + text).properties(title=f"Rescue Breakdown by {evidence_title}")
    _save_chart(chart, "48_rescue_breakdown", output_dir)
    return chart


def plot_rescue_by_filter(rescue_filter_df, output_dir: str,
                          evidence_col: str = "modality_evidence_caller"):
    """Chart 49: Rescued/non-rescued counts per FILTER category."""
    if rescue_filter_df is None or (hasattr(rescue_filter_df, 'is_empty') and rescue_filter_df.is_empty()):
        return
    # Fall back to RESCUED if evidence_col not in columns
    if evidence_col not in rescue_filter_df.columns:
        evidence_col = "RESCUED"
    pdf = rescue_filter_df.to_pandas()
    evidence_title = evidence_col.replace("_", " ").title()
    bars = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("FILTER:N", title="Classification"),
        y=alt.Y("count:Q", title="Variant Count", scale=_count_scale(), axis=_count_axis()),
        color=alt.Color(f"{evidence_col}:N", title=evidence_title,
                        scale=_color_scale(evidence_col) if evidence_col in _COLOR_REGISTRY else alt.Scale()),
        xOffset=f"{evidence_col}:N",
        tooltip=["FILTER", evidence_col, "count", "pct"],
    )
    text = _make_bar_text(pdf, x_field="FILTER", count_col="count",
                          x_offset=evidence_col, show_pct=True, group_col="FILTER")
    chart = (bars + text).properties(title=f"Rescue Status by Classification (FILTER)")
    _save_chart(chart, "49_rescue_by_filter", output_dir)
    return chart


def plot_rescue_cross_tab_heatmap(cross_tab_df, output_dir: str,
                                  evidence_col: str = "modality_evidence_caller"):
    """Chart 50: Modality evidence × FILTER × set_number heatmap with text labels."""
    if cross_tab_df is None or (hasattr(cross_tab_df, 'is_empty') and cross_tab_df.is_empty()):
        return
    # Fall back to RESCUED if evidence_col not in columns
    if evidence_col not in cross_tab_df.columns:
        evidence_col = "RESCUED"
    pdf = cross_tab_df.to_pandas()
    has_set = "set_number" in pdf.columns
    y_col = "FILTER"
    evidence_title = evidence_col.replace("_", " ").title()
    base = alt.Chart(pdf)
    rect = base.mark_rect().encode(
        x=alt.X(f"{evidence_col}:N", title=evidence_title),
        y=alt.Y(f"{y_col}:N", title="Classification"),
        color=alt.Color("count:Q", title="Count", scale=alt.Scale(scheme="blues", type="log")),
        tooltip=list(pdf.columns),
    )
    median_count = max(float(pdf["count"].median()), 1)
    text = base.mark_text(baseline="middle", fontSize=8).encode(
        x=alt.X(f"{evidence_col}:N"),
        y=alt.Y(f"{y_col}:N"),
        text=alt.Text("count:Q", format=",d"),
        color=alt.condition(f"datum.count > {median_count}", alt.value("white"), alt.value("black")),
    )
    chart = (rect + text).properties(title=f"Modality Evidence Cross-Tabulation", width=150, height=200)
    if has_set:
        chart = chart.facet(facet=alt.Facet("set_number:N", title="Set"), columns=2)
    _save_chart(chart, "50_rescue_cross_tab_heatmap", output_dir)
    return chart


def plot_rescue_vaf_boxplot(df, output_dir: str,
                            evidence_col: str = "modality_evidence_caller"):
    """Chart 51: DNA/RNA VAF distributions by modality evidence category."""
    if df is None or (hasattr(df, 'is_empty') and df.is_empty()):
        return
    # Fall back to RESCUED if evidence_col not in columns
    if not _has_column(df, evidence_col):
        if not _has_column(df, "RESCUED"):
            return
        evidence_col = "RESCUED"
    vaf_cols = [c for c in ["DNA_VAF_mean", "RNA_VAF_mean"] if _has_column(df, c)]
    if not vaf_cols:
        return
    select_cols = vaf_cols + [evidence_col]
    pdf = df.select(select_cols)
    # Fill null evidence values
    if isinstance(pdf, pl.LazyFrame):
        pdf = _sample_if_large(pdf, max_rows=50000)
    else:
        pdf = pdf.drop_nulls(subset=vaf_cols)
        pdf = _sample_if_large(pdf, max_rows=50000)
    pdf = pdf.to_pandas()
    melted = pdf.melt(id_vars=[evidence_col], var_name="modality", value_name="VAF")
    melted["modality"] = melted["modality"].str.replace("_VAF_mean", "")
    melted["VAF"] = melted["VAF"].clip(0, 1)
    evidence_title = evidence_col.replace("_", " ").title()
    color_scale = _color_scale(evidence_col) if evidence_col in _COLOR_REGISTRY else alt.Scale()
    chart = alt.Chart(melted).mark_boxplot(size=40).encode(
        x=alt.X(f"{evidence_col}:N", title=evidence_title),
        y=alt.Y("VAF:Q", title="Mean VAF", scale=alt.Scale(domain=[0, 1])),
        color=alt.Color(f"{evidence_col}:N", scale=color_scale),
        column=alt.Column("modality:N", title="Modality"),
    ).properties(title=f"VAF Distribution by {evidence_title}", width=200)
    _save_chart(chart, "51_rescue_vaf_boxplot", output_dir)
    return chart


def plot_rescue_dp_boxplot(df, output_dir: str,
                            evidence_col: str = "modality_evidence_caller"):
    """Chart 52: DNA/RNA DP distributions by modality evidence category."""
    if df is None or (hasattr(df, 'is_empty') and df.is_empty()):
        return
    # Fall back to RESCUED if evidence_col not in columns
    if not _has_column(df, evidence_col):
        if not _has_column(df, "RESCUED"):
            return
        evidence_col = "RESCUED"
    dp_cols = [c for c in ["DNA_DP_mean", "RNA_DP_mean"] if _has_column(df, c)]
    if not dp_cols:
        return
    select_cols = dp_cols + [evidence_col]
    pdf = df.select(select_cols)
    if isinstance(pdf, pl.LazyFrame):
        pdf = _sample_if_large(pdf, max_rows=50000)
    else:
        pdf = pdf.drop_nulls(subset=dp_cols)
        pdf = _sample_if_large(pdf, max_rows=50000)
    pdf = pdf.to_pandas()
    melted = pdf.melt(id_vars=[evidence_col], var_name="modality", value_name="DP")
    melted["modality"] = melted["modality"].str.replace("_DP_mean", "")
    melted["DP"] = melted["DP"].clip(0, 2000)
    evidence_title = evidence_col.replace("_", " ").title()
    color_scale = _color_scale(evidence_col) if evidence_col in _COLOR_REGISTRY else alt.Scale()
    chart = alt.Chart(melted).mark_boxplot(size=40).encode(
        x=alt.X(f"{evidence_col}:N", title=evidence_title),
        y=alt.Y("DP:Q", title="Mean DP (capped 2000)", scale=alt.Scale(domain=[0, 2000])),
        color=alt.Color(f"{evidence_col}:N", scale=color_scale),
        column=alt.Column("modality:N", title="Modality"),
    ).properties(title=f"DP Distribution by {evidence_title}", width=200)
    _save_chart(chart, "52_rescue_dp_boxplot", output_dir)
    return chart


def plot_rescue_sample_distribution(sample_rescue_df, output_dir: str,
                                    evidence_col: str = "modality_evidence_caller"):
    """Chart 53: Per-sample rescue counts, faceted by set using 2×2 grid."""
    if sample_rescue_df is None or (hasattr(sample_rescue_df, 'is_empty') and sample_rescue_df.is_empty()):
        return
    # Fall back to RESCUED if evidence_col not in columns
    if evidence_col not in sample_rescue_df.columns:
        evidence_col = "RESCUED"
    # Aggregate to sample × evidence level (drop FILTER detail)
    group_cols = ["sample_id", evidence_col]
    has_set = "set_number" in sample_rescue_df.columns
    if has_set:
        group_cols.append("set_number")
    pdf = sample_rescue_df.group_by(group_cols).agg(pl.col("count").sum()).to_pandas()
    evidence_title = evidence_col.replace("_", " ").title()
    base = alt.Chart(pdf)
    bars = base.mark_bar().encode(
        x=alt.X("sample_id:N", title="Sample", sort=None,
                 axis=alt.Axis(labelAngle=-45, labelLimit=100)),
        y=alt.Y("count:Q", title="Variant Count", scale=_count_scale(), axis=_count_axis(), stack="zero"),
        color=alt.Color(f"{evidence_col}:N", title=evidence_title,
                        scale=_color_scale(evidence_col) if evidence_col in _COLOR_REGISTRY else alt.Scale()),
        tooltip=["sample_id", evidence_col, "count"],
    )
    text = _make_bar_text(pdf, x_field="sample_id", count_col="count",
                          stack="zero", show_pct=True, group_col="sample_id",
                          base=base, overlap="hide")
    chart = (bars + text).properties(title=f"Per-Sample {evidence_title} Counts", height=300)
    if has_set:
        pdf["set_number"] = pdf["set_number"].astype(str)
        chart = chart.facet(
            facet=alt.Facet("set_number:N", title="Set"),
            columns=2,
        ).resolve_scale(x="independent", y="shared")
    _save_chart(chart, "53_rescue_sample_distribution", output_dir)
    return chart


def plot_rescue_by_tier(rescue_tier_df, output_dir: str,
                        evidence_col: str = "modality_evidence_caller"):
    """Chart 54: Rescue rate per CxDy tier, stacked bar."""
    if rescue_tier_df is None or (hasattr(rescue_tier_df, 'is_empty') and rescue_tier_df.is_empty()):
        return
    # Fall back to RESCUED if evidence_col not in columns
    if evidence_col not in rescue_tier_df.columns:
        evidence_col = "RESCUED"
    pdf = rescue_tier_df.to_pandas()
    evidence_title = evidence_col.replace("_", " ").title()
    bars = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("final_tier:N", title="Caller Tier"),
        y=alt.Y("count:Q", title="Variant Count", scale=_count_scale(), axis=_count_axis(), stack="zero"),
        color=alt.Color(f"{evidence_col}:N", title=evidence_title,
                        scale=_color_scale(evidence_col) if evidence_col in _COLOR_REGISTRY else alt.Scale()),
        tooltip=["final_tier", evidence_col, "count", "pct"],
    )
    text = _make_bar_text(pdf, x_field="final_tier", count_col="count",
                          stack="zero", show_pct=True, group_col="final_tier")
    chart = (bars + text).properties(title=f"Modality Evidence by Caller Tier")
    _save_chart(chart, "54_rescue_by_tier", output_dir)
    return chart


def plot_rescue_rate_trend(breakdown_df, output_dir: str,
                            evidence_col: str = "modality_evidence_caller"):
    """Chart 55: Rescue proportion (%) across sets, line chart."""
    if breakdown_df is None or (hasattr(breakdown_df, 'is_empty') and breakdown_df.is_empty()):
        return
    # Fall back to RESCUED if evidence_col not in columns
    if evidence_col not in breakdown_df.columns:
        evidence_col = "RESCUED"
    pdf = breakdown_df.to_pandas()
    if "set_number" not in pdf.columns:
        return
    evidence_title = evidence_col.replace("_", " ").title()
    # For backward compat with RESCUED, filter to "YES"; for modality_evidence_caller,
    # filter to "cross_modality" as the "rescued" analog
    filter_value = "cross_modality" if evidence_col == "modality_evidence_caller" else "YES"
    rescued_only = pdf[pdf[evidence_col] == filter_value].copy()
    if rescued_only.empty:
        return
    chart = alt.Chart(rescued_only).mark_line(point=True, strokeWidth=2).encode(
        x=alt.X("set_number:N", title="Set"),
        y=alt.Y("pct:Q", title="Rescue Rate (%)"),
        tooltip=["set_number", "count", "total", "pct"],
    ).properties(title="Rescue Rate Trend Across Sets")
    _save_chart(chart, "55_rescue_rate_trend", output_dir)
    return chart


def plot_rescue_caller_support(caller_support_df, output_dir: str,
                               evidence_col: str = "modality_evidence_caller"):
    """Chart 56: N_SUPPORT_CALLERS distribution by modality evidence status."""
    if caller_support_df is None or (hasattr(caller_support_df, 'is_empty') and caller_support_df.is_empty()):
        return
    # Fall back to RESCUED if evidence_col not in columns
    if evidence_col not in caller_support_df.columns:
        evidence_col = "RESCUED"
    pdf = caller_support_df.to_pandas()
    pdf["N_SUPPORT_CALLERS"] = pdf["N_SUPPORT_CALLERS"].astype(int).astype(str)
    evidence_title = evidence_col.replace("_", " ").title()
    bars = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("N_SUPPORT_CALLERS:N", title="Number of Supporting Callers"),
        y=alt.Y("count:Q", title="Variant Count", scale=_count_scale(), axis=_count_axis()),
        color=alt.Color(f"{evidence_col}:N", title=evidence_title,
                        scale=_color_scale(evidence_col) if evidence_col in _COLOR_REGISTRY else alt.Scale()),
        xOffset=f"{evidence_col}:N",
        tooltip=["N_SUPPORT_CALLERS", evidence_col, "count"],
    )
    text = _make_bar_text(pdf, x_field="N_SUPPORT_CALLERS", count_col="count",
                          x_offset=evidence_col, show_pct=True,
                          group_col="N_SUPPORT_CALLERS")
    chart = (bars + text).properties(title=f"Caller Support Distribution by {evidence_title}")
    _save_chart(chart, "56_rescue_caller_support", output_dir)
    return chart


# ═══════════════════════════════════════════════════════════════════════════
# Multi-Allelic Visualization Charts (Section 2)
# ═══════════════════════════════════════════════════════════════════════════

MULTIALLELIC_CLASS_DOMAIN = ["normalization_artifact", "noise", "true_multi_allelic"]
MULTIALLELIC_CLASS_COLORS = ["#d62728", "#ff7f0e", "#2ca02c"]


def plot_multiallelic_classification(df, output_dir: str):
    """Multi-allelic classification: stacked bar chart by set_number/disease.

    Groups by multiallelic_class and set_number (or disease_normalized),
    filters to multiallelic_class != "single", and shows the 3 multi-allelic
    categories as a stacked bar chart.
    """
    if not _has_column(df, "multiallelic_class"):
        return

    group_col = "disease_normalized" if _has_column(df, "disease_normalized") else "set_number"
    if not _has_column(df, group_col):
        return

    counts = _maybe_collect(
        df.filter(pl.col("multiallelic_class") != "single")
        .group_by([group_col, "multiallelic_class"])
        .agg(pl.len().alias("count"))
        .sort(group_col)
    )
    if counts.is_empty():
        return
    pdf = counts.to_pandas()
    pdf[group_col] = pdf[group_col].astype(str)
    group_title = group_col.replace("_", " ").title()

    base = alt.Chart(pdf)
    bars = base.mark_bar().encode(
        x=alt.X(f"{group_col}:N", title=group_title),
        y=alt.Y("count:Q", title="Number of Multi-Allelic Variants",
                scale=_count_scale(), axis=_count_axis()),
        color=alt.Color("multiallelic_class:N", title="Multi-Allelic Class",
                        scale=alt.Scale(domain=MULTIALLELIC_CLASS_DOMAIN,
                                        range=MULTIALLELIC_CLASS_COLORS)),
    )
    text = _make_bar_text(pdf, x_field=group_col, count_col="count",
                          stack="zero", show_pct=True, group_col=group_col,
                          base=base, color_col="multiallelic_class")
    chart = (bars + text).properties(title=f"Multi-Allelic Classification by {group_title}")
    _save_chart(chart, "classification", output_dir)
    return chart


def plot_allele_balance_scatter(df, output_dir: str):
    """Multi-allelic allele balance scatter: major_vaf vs minor_vaf.

    Filters to n_alleles_at_site >= 2, computes per-site major and minor VAF,
    colors by multiallelic_class. Adds y=x reference line. Samples to 5000 points.
    """
    if not _has_column(df, "n_alleles_at_site"):
        return

    vaf_cols = [c for c in df.columns if c.endswith("_VAF") and c.startswith(("DNA_", "RNA_"))]
    if len(vaf_cols) < 2:
        return

    select_cols = ["n_alleles_at_site", "multiallelic_class"] + vaf_cols
    color_col = "multiallelic_class" if _has_column(df, "multiallelic_class") else None
    if color_col:
        select_cols = list(dict.fromkeys(select_cols))  # deduplicate

    pdf = df.select(select_cols).filter(pl.col("n_alleles_at_site") >= 2)
    pdf = _sample_if_large(pdf, max_rows=5000)
    if pdf.is_empty():
        return
    pdf = pdf.to_pandas()

    # Compute major and minor VAF per site from per-caller VAF columns
    vaf_existing = [c for c in vaf_cols if c in pdf.columns]
    if not vaf_existing:
        return

    def row_vaf_stats(row):
        vafs = [row[c] for c in vaf_existing if pd.notna(row[c])]
        if len(vafs) >= 2:
            return pd.Series({"major_vaf": max(vafs), "minor_vaf": min(vafs)})
        return pd.Series({"major_vaf": float("nan"), "minor_vaf": float("nan")})

    vaf_df = pdf[vaf_existing].apply(row_vaf_stats, axis=1)
    pdf["major_vaf"] = vaf_df["major_vaf"].clip(0, 1)
    pdf["minor_vaf"] = vaf_df["minor_vaf"].clip(0, 1)
    pdf = pdf.dropna(subset=["major_vaf", "minor_vaf"])

    if pdf.empty:
        return

    # y=x reference line
    ref_line = alt.Chart(pd.DataFrame({"x": [0, 1], "y": [0, 1]})).mark_rule(
        strokeDash=[4, 4], opacity=0.5, color="gray"
    ).encode(x=alt.X("x:Q"), y=alt.Y("y:Q"))

    enc = {
        "x": alt.X("major_vaf:Q", title="Major Allele VAF", scale=alt.Scale(domain=[0, 1])),
        "y": alt.Y("minor_vaf:Q", title="Minor Allele VAF", scale=alt.Scale(domain=[0, 1])),
    }
    if color_col and color_col in pdf.columns:
        enc["color"] = alt.Color(f"{color_col}:N", title="Multi-Allelic Class",
                                 scale=alt.Scale(domain=MULTIALLELIC_CLASS_DOMAIN,
                                                 range=MULTIALLELIC_CLASS_COLORS))

    scatter = alt.Chart(pdf).mark_circle(opacity=0.5, size=30).encode(**enc)
    chart = (scatter + ref_line).properties(title="Multi-Allelic Allele Balance (major vs minor VAF)")
    _save_chart(chart, "allele_balance", output_dir)
    return chart


def plot_vaf_sum_histogram(df, output_dir: str):
    """Multi-allelic VAF sum histogram.

    Filters to n_alleles_at_site >= 2, histograms vaf_sum with bin=20.
    Adds vertical rule at x=1.0. Color bins: >1.0 red, <=1.0 blue.
    """
    if not _has_column(df, "n_alleles_at_site") or not _has_column(df, "vaf_sum"):
        return

    pdf = df.select(["n_alleles_at_site", "vaf_sum"]).filter(
        pl.col("n_alleles_at_site") >= 2
    ).drop_nulls(subset=["vaf_sum"])
    pdf = _maybe_collect(pdf)
    if pdf.is_empty():
        return
    pdf = pdf.to_pandas()

    # Bin and color
    pdf["vaf_bucket"] = pdf["vaf_sum"].apply(lambda x: ">1.0" if x > 1.0 else "≤1.0")

    ref_rule = alt.Chart(pd.DataFrame({"x": [1.0]})).mark_rule(
        strokeDash=[6, 4], color="red", strokeWidth=2
    ).encode(x=alt.X("x:Q"))

    hist = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("vaf_sum:Q", bin=alt.Bin(maxbins=20), title="VAF Sum Across Alleles"),
        y=alt.Y("count()", title="Number of Sites"),
        color=alt.Color("vaf_bucket:N", title="VAF Sum",
                        scale=alt.Scale(domain=[">1.0", "≤1.0"],
                                        range=["#d62728", "#1f77b4"])),
    )

    chart = (hist + ref_rule).properties(title="Multi-Allelic VAF Sum Distribution (sites with >=2 alleles)")
    _save_chart(chart, "vaf_sum", output_dir)
    return chart


def plot_category_conflict_summary(df, output_dir: str):
    """Multi-allelic category conflict summary: horizontal bar of top 10 conflict types.

    Filters to category_conflict == True, groups by (CHROM, POS), collects unique
    FILTER values as sorted string. Top 10 conflict patterns shown as horizontal bars.
    """
    if not _has_column(df, "category_conflict") or not _has_column(df, "FILTER"):
        return

    pdf = _maybe_collect(
        df.filter(pl.col("category_conflict") == True)
        .select(["CHROM", "POS", "FILTER"])
    )
    if pdf.is_empty():
        return

    # Group by (CHROM, POS), collect sorted unique FILTER values
    conflicts = (
        pdf.group_by(["CHROM", "POS"])
        .agg(pl.col("FILTER").unique().sort().str.join("/").alias("conflict_type"))
    )
    # Count per conflict type
    counts = (
        conflicts.group_by("conflict_type")
        .agg(pl.len().alias("count"))
        .sort("count", descending=True)
        .head(10)
    )
    if counts.is_empty():
        return
    counts = counts.to_pandas()

    bars = alt.Chart(counts).mark_bar().encode(
        y=alt.Y("conflict_type:N", title="Conflict Type (FILTER values at site)",
                sort="-x"),
        x=alt.X("count:Q", title="Number of Sites", scale=_count_scale(), axis=_count_axis()),
    )
    text = _make_bar_text(counts, x_field="count", count_col="count",
                          x_sort=None, show_pct=False)
    # For horizontal bar, swap x/y in text encoding
    text = alt.Chart(counts).mark_text(dy=-8, fontSize=9).encode(
        y=alt.Y("conflict_type:N", sort="-x"),
        x=alt.X("count:Q"),
        text=alt.Text("count:Q", format=",d"),
        color=alt.value("#333"),
    )
    chart = (bars + text).properties(
        title="Top 10 Multi-Allelic Category Conflicts (FILTER disagreement at same site)")
    _save_chart(chart, "category_conflict", output_dir)
    return chart


# ═══════════════════════════════════════════════════════════════════════════
# BAM Visualization Charts (Section 5)
# ═══════════════════════════════════════════════════════════════════════════

def plot_bam_metrics_sample_wise(bam_stats_df, output_dir: str):
    """Faceted bar chart of BAM metrics per sample, colored by bam_type.

    Takes bam_stats_df (not combined_df). Melts metrics into long form
    and creates a faceted bar chart with one panel per metric. Samples
    are ordered by set_number then sample_id.
    """
    if bam_stats_df is None or (hasattr(bam_stats_df, 'is_empty') and bam_stats_df.is_empty()):
        return
    if "sample_id" not in bam_stats_df.columns or "bam_type" not in bam_stats_df.columns:
        return

    # Identify metric columns (numeric, not sample_id/bam_type/set_number)
    skip_cols = {"sample_id", "bam_type", "set_number", "disease", "disease_normalized"}
    metric_cols = [c for c in bam_stats_df.columns if c not in skip_cols
                   and bam_stats_df[c].dtype in (pl.Int64, pl.Float64, pl.Int32, pl.Float32)]

    if not metric_cols:
        return

    # Build melted long-form DataFrame
    id_cols = ["sample_id", "bam_type"]
    has_set = "set_number" in bam_stats_df.columns
    if has_set:
        id_cols.append("set_number")

    melted = bam_stats_df.select(id_cols + metric_cols).unpivot(
        index=id_cols, variable_name="metric", value_name="value"
    )

    # Order samples by set_number then sample_id
    if has_set:
        melted = melted.sort(["set_number", "sample_id"])
    else:
        melted = melted.sort("sample_id")

    pdf = melted.to_pandas()
    pdf["sample_id"] = pdf["sample_id"].astype(str)
    if has_set:
        pdf["set_number"] = pdf["set_number"].astype(str)

    # Faceted bar chart: one panel per metric
    n_metrics = len(metric_cols)
    n_cols = min(3, n_metrics)

    chart = alt.Chart(pdf).mark_bar().encode(
        x=alt.X("sample_id:N", title="Sample", sort=None,
                axis=alt.Axis(labelAngle=-45, labelLimit=100)),
        y=alt.Y("value:Q", title="Value"),
        color=alt.Color("bam_type:N", title="BAM Type", scale=_color_scale("bam_type")),
        column=alt.Column("metric:N", title="Metric"),
    ).properties(
        title="BAM Metrics per Sample by BAM Type",
        width=alt.Step(12),
    )
    chart = chart.resolve_scale(x="independent", y="independent")

    _save_chart(chart, "metrics_sample_wise", output_dir)
    return chart


def plot_bam_coverage_distribution(bam_stats_df, output_dir: str):
    """Grouped bar chart of mean coverage bin percentages by bam_type and set_number.

    Groups by bam_type and set_number, computes mean of coverage-related
    percentage columns (cov_1x_pct through cov_100x_pct). Creates a grouped
    bar chart with bars grouped by bam_type and colored by coverage bin.
    """
    if bam_stats_df is None or (hasattr(bam_stats_df, 'is_empty') and bam_stats_df.is_empty()):
        return

    # Find coverage percentage columns
    cov_pct_cols = [c for c in bam_stats_df.columns
                    if c.endswith("_pct") and "cov_" in c]
    if not cov_pct_cols:
        return

    has_set = "set_number" in bam_stats_df.columns
    group_cols = ["bam_type"]
    if has_set:
        group_cols.append("set_number")

    # Group by bam_type (and set_number) and compute mean of each cov_pct column
    agg_exprs = [pl.col(c).mean().alias(f"mean_{c}") for c in cov_pct_cols]
    grouped = bam_stats_df.group_by(group_cols).agg(agg_exprs)

    if has_set:
        grouped = grouped.with_columns(pl.col("set_number").cast(pl.Utf8))

    # Melt mean coverage columns for plotting
    mean_cols = [f"mean_{c}" for c in cov_pct_cols]
    melted = grouped.unpivot(
        index=group_cols, variable_name="coverage_bin", value_name="mean_pct"
    )
    # Clean up bin names: strip "mean_cov_" prefix and "_pct" suffix
    melted = melted.with_columns(
        pl.col("coverage_bin").str.replace("mean_cov_", "").str.replace("_pct", "")
    )

    pdf = melted.to_pandas()
    if pdf.empty:
        return

    base = alt.Chart(pdf)
    bars = base.mark_bar().encode(
        x=alt.X("bam_type:N", title="BAM Type",
                axis=alt.Axis(labelAngle=0)),
        y=alt.Y("mean_pct:Q", title="Mean % of Bases Covered"),
        xOffset="coverage_bin:N",
        color=alt.Color("coverage_bin:N", title="Coverage Threshold"),
        tooltip=["bam_type", "coverage_bin", "mean_pct"],
    )
    chart = bars.properties(
        title="Mean Coverage Distribution by BAM Type",
        width=250,
    )
    if has_set:
        chart = chart.facet(
            facet=alt.Facet("set_number:N", title="Set"),
            columns=2,
        ).resolve_scale(x="independent", y="shared")

    _save_chart(chart, "coverage_distribution", output_dir)
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
