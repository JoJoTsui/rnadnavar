## Context

See `fix-cross-sample-aggregation-oom` and `fix-visualizer-lazy-frame-and-oom` for prior work. This change builds on those fixes and addresses remaining data quality and visualization bugs discovered during deep code review of the 65-sample pipeline run.

## Goals / Non-Goals

**Goals:** Fix silent data failures from column pruning, parameterize all 14 chart functions for wise-based generation, add 5 new chart types, apply quality fixes (%, log scale, box+violin), implement CLI wise loop producing ~105 chart×wise combinations, comprehensive tests, shell script documentation.

**Non-Goals:** Changing rescue VCF semantics, modifying caller VCF parsing logic, adding new callers, changing the Nextflow pipeline itself.

---

## Part 1: Data Foundation — `_CROSS_SAMPLE_COLS` Fix

### Problem

`_ensure_eager()` filters every cross-sample operation to `_CROSS_SAMPLE_COLS` (40 columns). Per-caller VAF/DP/AD columns (24 columns) are absent, causing:

| Function | Symptom |
|----------|---------|
| `compute_caller_wise_summary()` | All 6 callers: `n_with_vaf=0, n_with_dp=0` |
| `compute_vaf_threshold_sweep()` | Empty DataFrame → `vaf_threshold_sweep.tsv` missing → Chart 31 missing |
| `plot_vaf_distribution()` | Per-caller VAF columns pruned → empty |
| `plot_vaf_boxplot_per_tier()` | Same |
| `plot_dp_boxplot_per_tier()` | Same |
| `plot_dna_vs_rna_per_caller()` | Same |

### Fix

Add 24 per-caller columns to `_CROSS_SAMPLE_COLS`:

```python
# Per-caller VAF columns (6)
"DNA_mutect2_VAF", "RNA_mutect2_VAF",
"DNA_deepsomatic_VAF", "RNA_deepsomatic_VAF",
"DNA_strelka_VAF", "RNA_strelka_VAF",

# Per-caller DP columns (6)
"DNA_mutect2_DP", "RNA_mutect2_DP",
"DNA_deepsomatic_DP", "RNA_deepsomatic_DP",
"DNA_strelka_DP", "RNA_strelka_DP",

# Per-caller AD columns (12)
"DNA_mutect2_AD_REF", "DNA_mutect2_AD_ALT",
"RNA_mutect2_AD_REF", "RNA_mutect2_AD_ALT",
"DNA_deepsomatic_AD_REF", "DNA_deepsomatic_AD_ALT",
"RNA_deepsomatic_AD_REF", "RNA_deepsomatic_AD_ALT",
"DNA_strelka_AD_REF", "DNA_strelka_AD_ALT",
"RNA_strelka_AD_REF", "RNA_strelka_AD_ALT",
```

Also fix the PerformanceWarning in `_ensure_eager` by using `df.collect_schema().names()` instead of `df.columns` for LazyFrame column checks.

---

## Part 2: Chart Function Parameterization

### 2.1 Add `group_col` Parameter (10 functions)

These functions currently hardcode `"set_number"` and need a `group_col` parameter with dynamic titles:

```python
# Before
def plot_caller_overlap(df, output_dir: str):
    counts = df.group_by(["set_number", "final_tier"]).agg(...)
    # title: "Variant Tier Distribution by Set"

# After
def plot_caller_overlap(df, output_dir: str, group_col: str = "set_number"):
    counts = df.group_by([group_col, "final_tier"]).agg(...)
    # title: f"Variant Tier Distribution by {group_col.replace('_',' ').title()}"
```

Functions: `plot_caller_overlap`, `plot_variant_type_distribution`, `plot_ti_tv_ratio`, `plot_cross_modality`, `plot_dna_vs_rna_dp`, `plot_ref_alt_dp_scatter`, `plot_filter_distribution`, `plot_redi_evidence`, `plot_gt_concordance`, `plot_cosmic_gnomad_annotation`

`plot_vc_distribution` already has `group_col` param — only needs title fix.

### 2.2 Add `color_col` Parameter (5 functions)

These functions need wise-based coloring instead of hardcoded color sources:

```python
# Before
def plot_vaf_distribution(df, output_dir: str):
    # color by caller name (from melt)

# After
def plot_vaf_distribution(df, output_dir: str, color_col: str = None):
    # if color_col: add color encoding by color_col
```

Functions: `plot_vaf_distribution`, `plot_dna_vs_rna_vaf`, `plot_dna_vs_rna_per_caller`, `plot_bam_coverage_violin`, `plot_caller_concordance_vs_vaf`

### 2.3 Add `facet_col` Parameter (6 functions)

These tier-based functions need additional wise faceting:

```python
# Before
def plot_vaf_boxplot_per_tier(df, output_dir: str):
    chart = ... column=alt.Column("caller_tier:N", ...)

# After
def plot_vaf_boxplot_per_tier(df, output_dir: str, facet_col: str = None):
    if facet_col:
        chart = ... column=alt.Column(f"{facet_col}:N", ...)
        # nested: outer=facet_col, inner=caller_tier
```

Functions: `plot_vaf_boxplot_per_tier`, `plot_dp_boxplot_per_tier`, `plot_gt_concordance_per_tier`, `plot_tiered_caller_overlap`, `plot_tiered_variant_types`, `plot_tier_quality_distribution`

---

## Part 3: Chart Quality Fixes

### 3.1 Percentage Text Marks (4 functions)

Apply `_plot_bar_wise(pct=True)` pattern: after creating bar chart, overlay text marks showing percentage values.

Functions: `plot_tiered_caller_overlap`, `plot_tiered_variant_types`, `plot_filter_distribution`, `plot_variant_type_distribution`

### 3.2 Log Scale

`plot_tiered_caller_overlap`: y-axis should use log scale (`scale=alt.Scale(type="log")`) because N_SUPPORT_CALLERS counts span orders of magnitude (C1D1 has thousands, C7D0 has millions).

### 3.3 Per-Sample Chart

`plot_per_sample_distribution`: remove `top_n=30`, enable horizontal scrolling for all 65 samples via `alt.Chart(...).properties(height=len(samples)*15)` with `alt.Scale(padding=...)`.

### 3.4 Box+Violin Overlays

Wire existing VAF/DP distribution charts through `_plot_box_violin_wise` factory instead of plain `mark_boxplot()`:

Functions: `plot_vaf_distribution` (chart 03), new `plot_dp_distribution` (chart 12), `plot_vaf_boxplot_per_tier` (chart 14), `plot_dp_boxplot_per_tier` (chart 15), `plot_bam_coverage_violin` (chart 21)

---

## Part 4: New Chart Functions

### Chart 12: `plot_dp_distribution(df, output_dir, color_col=None)`
Per-caller DP box+violin. Melts `{caller}_DP` columns, same structure as `plot_vaf_distribution` but for DP values. Generated for set, disease, tier, caller wises.

### Chart 35: `plot_mean_vaf_per_group(summary_df, output_dir, group_col, wise_name)`
Bar chart from pre-computed wise summary data. Uses `_plot_bar_wise` factory. One bar per group showing mean DNA/RNA VAF side-by-side.

### Chart 36: `plot_mean_dp_per_group(summary_df, output_dir, group_col, wise_name)`
Same structure as chart 35 but for mean DP values.

### Chart 37: `plot_n_support_callers_dist(df, output_dir, group_col)`
N_SUPPORT_CALLERS histogram (1-6) faceted or colored by `group_col`. Generated for set, disease, tier, chromosome wises.

### Chart 38: `plot_caller_tier_heatmap(df, output_dir, facet_col=None)`
Caller (6) × Tier (C1D1..C7D0) detection matrix. Rows = callers, columns = tiers, color = % of variants detected by that caller in that tier. Faceted by set or disease.

### Chart 39: `plot_sample_overview_scatter(sample_stats_df, output_dir)`
Global overview: mean VAF vs mean DP per sample, colored by set_number or disease_normalized. Point size proportional to total_variants.

---

## Part 5: CLI Wise Loop Architecture

### Chart Registry

```python
WISE_CHART_REGISTRY = {
    "set": [
        (plot_vc_distribution, {"group_col": "set_number"}),
        (plot_caller_overlap, {"group_col": "set_number"}),
        (plot_vaf_distribution, {"color_col": "set_number"}),
        (plot_dp_distribution, {"color_col": "set_number"}),
        (plot_dna_vs_rna_vaf, {"color_col": "set_number"}),
        (plot_dna_vs_rna_dp, {"group_col": "set_number"}),
        (plot_gt_concordance, {"group_col": "set_number"}),
        (plot_cosmic_gnomad_annotation, {"group_col": "set_number"}),
        (plot_variant_type_distribution, {"group_col": "set_number"}),
        (plot_ti_tv_ratio, {"group_col": "set_number"}),
        (plot_cross_modality, {"group_col": "set_number"}),
        (plot_ref_alt_dp_scatter, {"group_col": "set_number"}),
        (plot_filter_distribution, {"group_col": "set_number"}),
        (plot_redi_evidence, {"group_col": "set_number"}),
        (plot_caller_agreement_matrix, {}),  # per-set
        (plot_tier_quality_distribution, {"facet_col": "set_number"}),
        (plot_caller_concordance_vs_vaf, {"color_col": "set_number"}),
        (plot_n_support_callers_dist, {"group_col": "set_number"}),
        (plot_mean_vaf_per_group, {"group_col": "set_number"}),
        (plot_mean_dp_per_group, {"group_col": "set_number"}),
        # ~22 entries
    ],
    "disease": [
        # ~18 entries — similar to set but fewer (no caller_overlap, fewer facets)
    ],
    "sample": [
        (plot_per_sample_distribution, {}),
        (plot_per_sample_tier_distribution, {}),
        (plot_vc_distribution, {"group_col": "sample_id"}),
        (plot_ti_tv_ratio, {"group_col": "sample_id"}),
        (plot_cross_modality, {"group_col": "sample_id"}),
        (plot_cosmic_gnomad_annotation, {"group_col": "sample_id"}),
        # ~8 entries
    ],
    "tier": [
        (plot_vc_distribution, {"group_col": "final_tier"}),
        (plot_vaf_boxplot_per_tier, {}),
        (plot_dp_boxplot_per_tier, {}),
        (plot_gt_concordance_per_tier, {}),
        (plot_tiered_caller_overlap, {}),
        (plot_tiered_variant_types, {}),
        (plot_per_tier_vaf_boxplot, {}),
        (plot_database_enrichment_by_tier, {}),
        # ~15 entries
    ],
    "caller": [
        (plot_vaf_distribution, {}),
        (plot_dp_distribution, {}),
        (plot_caller_agreement_matrix, {}),
        (plot_dna_vs_rna_per_caller, {}),
        (plot_caller_tier_heatmap, {}),
        # ~8 entries
    ],
    "chromosome": [
        (plot_chromosome_density, {}),
        (plot_vc_distribution, {"group_col": "CHROM"}),
        (plot_variant_type_distribution, {"group_col": "CHROM"}),
        (plot_ti_tv_ratio, {"group_col": "CHROM"}),
        (plot_n_support_callers_dist, {"group_col": "CHROM"}),
        # ~8 entries
    ],
}
```

### Generation Loop

```python
for wise_name in active_wises:
    wise_dir = output_dir / "plots" / wise_name
    wise_dir.mkdir(parents=True, exist_ok=True)
    for chart_fn, kwargs in WISE_CHART_REGISTRY.get(wise_name, []):
        try:
            fig = chart_fn(combined_df, str(wise_dir), **kwargs)
            if fig is not None:
                figs.append(fig)
        except Exception as e:
            print(f"  WARNING: {chart_fn.__name__} ({wise_name}) failed: {e}")
```

### `_save_chart` Update

```python
def _save_chart(chart, name, output_dir):
    plots_dir = Path(output_dir) / "plots"
    # output_dir already includes /plots/{wise}/, so the chart goes to the right place
    plots_dir.mkdir(parents=True, exist_ok=True)
    ...
```

Wait — if output_dir is `stats/full/plots/set/`, then `_save_chart` would create `stats/full/plots/set/plots/`. Need to adjust: pass the wise plot directory directly, not the base output dir.

Fix: `output_dir` parameter of chart functions becomes the plot directory itself (e.g., `stats/full/plots/set/`), not the base output dir. `_save_chart` writes directly into it.

---

## Part 6: Shell Script & Documentation

Create `examples/seq2neo/run_stats.sh`:

```bash
#!/usr/bin/env bash
# Run variant statistics pipeline on seq2neo output.
#
# Usage:
#   bash run_stats.sh                          # all sets, all wises
#   bash run_stats.sh --set 1                  # set 1 only
#   bash run_stats.sh --wise set tier caller   # specific wises
#   bash run_stats.sh --bed /path/to/exome.bed # WES coverage
#   bash run_stats.sh --no-pileup              # skip BAM pileup

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${HERE}/data/processed/sample_manifest.parquet"
OUTDIR="${HERE}/stats/full"

PYTHONPATH="${HERE}/../../bin" .venv/bin/python -m vcf_stats.seq2neo.cli \
    --manifest "$MANIFEST" \
    --output-dir "$OUTDIR" \
    --sample-workers 4 \
    --bam-workers 8 \
    "$@"
```

---

## Part 7: Chart Numbering (Final)

```
01 vc_distribution           — set, disease, sample, tier, chromosome (5 wises)
02 caller_overlap             — set, disease (2 wises)
03 vaf_distribution           — set, disease, tier, caller (4 wises)
04 dna_vs_rna_vaf             — set, disease, tier (3 wises)
05 dna_vs_rna_dp              — set, disease, tier (3 wises)
06 gt_concordance             — set, disease (2 wises)
07 cosmic_gnomad              — set, disease, sample, tier, chromosome (5 wises)
08 per_sample_distribution    — sample (1 wise)
09 variant_type_distribution  — set, disease, sample, tier, chromosome (5 wises)
10 ti_tv_ratio                — set, disease, sample, tier, chromosome (5 wises)
11 cross_modality             — set, disease, sample, tier, chromosome (5 wises)
12 dp_distribution            — NEW: set, disease, tier, caller (4 wises)
13 validation_heatmap         — global (1 wise)
14 vaf_violin_per_tier        — tier, set, disease (3 wises)
15 dp_per_tier                — tier, set, disease (3 wises)
16 gt_concordance_per_tier    — tier, set, disease (3 wises)
17 ref_alt_dp_scatter         — set, disease, tier (3 wises)
18 caller_overlap_per_tier    — tier, set, disease (3 wises) + % + log
19 variant_types_per_tier     — tier, set, disease (3 wises) + %
20 bam_metrics                — sample (1 wise)
21 bam_coverage_violin        — set, disease, tier (3 wises)
22 sample_tier_dist           — sample (1 wise)
23 per_tier_vaf               — REMOVED (consolidated into 25)
24 filter_distribution        — set, disease, tier, chromosome (4 wises) + %
25 per_tier_vaf               — tier (1 wise) [was 23+25 consolidated]
26 dna_vs_rna_per_caller      — set, disease, tier, caller (4 wises)
27 tier_quality               — set, disease (2 wises)
28 caller_agreement           — set, disease (2 wises)
29 redi_evidence              — set, disease, tier (3 wises)
30 chromosome_density         — chromosome (1 wise)
31 vaf_threshold_sweep        — global, set (2 wises)
32 caller_concordance_vs_vaf  — set, disease, tier (3 wises)
33 filter_effectiveness       — global, set, disease (3 wises)
34 db_enrichment_by_tier      — tier (1 wise)
35 mean_vaf_per_group         — NEW: set, disease, sample, tier, caller, chromosome (6 wises)
36 mean_dp_per_group          — NEW: set, disease, sample, tier, caller, chromosome (6 wises)
37 n_support_callers_dist     — NEW: set, disease, tier, chromosome (4 wises)
38 caller_tier_heatmap        — NEW: set, disease (2 wises)
39 sample_overview_scatter    — NEW: global (1 wise)

Total chart×wise: ~100 combinations
```

---

## Risks / Trade-offs

- **Memory**: Adding 24 per-caller columns to `_CROSS_SAMPLE_COLS` increases cross-sample memory from ~40 cols to ~64 cols. At 58M variants, this is ~11 GB extra. Acceptable — currently uses ~20 GB peak, well within 200 GB limit.
- **Chart count**: ~100 charts is a lot. Mitigation: `--wise` flag generates only requested wises; all wises generated by default.
- **Chart generation time**: Each chart scans parquet columns. polars query optimizer shares scans across charts reading from the same lazy frame. Should be efficient.
- **BAM pileup is slow**: Indexed BAM queries per variant × 3 BAM types. Mitigation: `--no-pileup` flag.
- **Test data**: New chart tests use synthetic polars DataFrames, not real VCF files. Fast and self-contained.
