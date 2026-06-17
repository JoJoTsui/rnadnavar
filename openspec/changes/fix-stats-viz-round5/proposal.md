## Why

After the round 4 pipeline run, 23 requirements were identified across 6 review agents. These map to 10 systematic patterns — not individual fixes. The core issues are: (1) no unified color registry for 16 categorical entities across ~40+ chart callsites, (2) broken faceting architecture using `alt.Column()` in `.encode()` (doesn't wrap into grids) affecting 17 functions, (3) `transform_density` violins produce empty SVGs in static export, (4) per-FILTER variant sub-plots computed but discarded, (5) ~20 count axes lack symlog scale, (6) rescue boxplots drop all non-rescued data, (7) 4/7 heatmaps lack text labels. The fix requires building 6 centralized helper functions first, then applying them systematically across all chart functions — not one-by-one patching.

## What Changes

### Infrastructure — 6 New Helpers
- **`_color_scale(entity)`**: Returns unified `alt.Scale(domain=..., range=...)` for any of 10+ categorical entities (callers, classifications, rescue, variant type, set, disease, BAM type, tier, modality, etc.)
- **`_apply_faceting(chart, group_col)`**: Applies correct faceting pattern based on group_col: set→2×2, disease→4×3, sample→x-independent, FILTER→3-col
- **`_add_heatmap_text(base, x_enc, y_enc, text_col, pdf, fontSize=8)`**: Adds auto-colored text overlay to any heatmap (white-on-dark, black-on-light)
- **`_add_bar_labels(chart, pdf, label_col, fmt, dy, fontSize)`**: Adds count/percentage labels to stacked bar charts
- **`_clip_dp(pdf, col, cap=2000)`**: Standard DP clipping with outlier count return
- **`_count_scale()`**: Returns `alt.Scale(type="symlog", constant=1)` for count axes

### Color Unification (Pattern 1)
- Define 10 global constant pairs (DOMAIN + COLORS) at module top for all categorical entities
- Replace all ~40+ ad-hoc `alt.Color(...)` references with `_color_scale()` calls

### Faceting Migration (Pattern 2)
- Migrate 11 functions from `alt.Column()` in `.encode()` to `.facet()` method with proper `columns=` wrapping
- Add `resolve_scale(x="independent")` for per-sample charts so each facet shows only its own samples
- Standardize: set→2×2, disease→4×3 (shared y), sample→independent x, FILTER→3-col

### Violin → Boxplot (Pattern 3)
- Remove `transform_density` + `mark_area` from 5 functions (empty in static SVG export)
- Replace with `mark_boxplot()` only (natively supported in vl-convert)

### Per-FILTER Sub-Plots + Variant-Category-Wise (Pattern 4)
- Stop filtering out per-classification data in threshold sweep chart functions
- Add per-FILTER faceting to threshold sweeps, per-caller, and rescue charts
- Add new `variant-category` wise dimension with `group_col="FILTER"` — 11 charts using shared codebase

### Symlog Count Axes (Pattern 5)
- Apply `_count_scale()` to ~20 chart functions with variant count y-axes

### Rescue Data Bug Fix (Pattern 6)
- Fix `plot_rescue_vaf_boxplot` and `plot_rescue_dp_boxplot`: normalize RESCUED before `drop_nulls()`

### Text Labels on Heatmaps + Bars (Pattern 7)
- Add `_add_heatmap_text()` to 4 heatmaps currently missing text (#28, #13, #33, #38)
- Reduce fontSize from 9→7-8 in existing heatmaps (#45, #50)
- Include Somatic in FP cross-tab, rename to "Classification × Caller Support"
- Add `_add_bar_labels()` to rescue tier (#54) and caller support (#56) stacked bars

### BAM DP Clipping (Pattern 8)
- Add `_clip_dp()` to `plot_bam_dp_distribution` (missed in round 4)

### Resume Per-Sample Fallback (Pattern 9)
- Recompute `sample_summary` and `bam_stats` from parquet when TSVs missing on `--resume`

## Capabilities

### New Capabilities
- `variant-category-wise`: New wise dimension grouping by FILTER (Somatic/Germline/Reference/Artifact/RNAedit/NoConsensus) using shared chart codebase
- `visualization-helpers`: 6 centralized helper functions for colors, faceting, text labels, clipping, and count scales

### Modified Capabilities
- `variant-visualization`: Unified colors (16 entities), faceting migration (17 functions), violin→boxplot (5 functions), symlog counts (~20 functions), heatmap text (7 charts), bar labels, DP clipping
- `rescue-analytics`: Fix non-rescued data bug, add per-FILTER sub-plots, unified rescue coloring
- `depth-threshold-sweep`: Stop filtering per-classification data, add FILTER faceting

## Impact

- `bin/vcf_stats/seq2neo/visualizer.py` — 6 new helpers (~150 lines), ~50 chart functions modified (~500 lines changed)
- `bin/vcf_stats/seq2neo/cli.py` — variant-category-wise registry, resume fallback (~30 lines)
- `bin/vcf_stats/seq2neo/statistics.py` — Remove Somatic exclusion from FP cross-tab (~5 lines)
