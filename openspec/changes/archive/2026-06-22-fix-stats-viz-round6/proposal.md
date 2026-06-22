## Why

Six-agent adversarial review of 100 SVGs and 35 TSVs from the 61-sample pipeline run identified 7 systematic patterns affecting ~40 chart functions. The most critical: (1) `symlog` scale compresses 88% of chart height into the near-zero region, making bar comparisons impossible across decades; (2) chart 43 renders 3 overlapping text labels per cell due to un-aggregated partition data (a data correctness bug); (3) sample-wise charts lack `set_number` metadata entirely, preventing the requested per-set faceting. All issues trace to shared code patterns and can be fixed by pattern, not one-by-one.

## What Changes

### Pattern 1 — Log scale for count axes (1 line → 20 functions → 40+ charts)
- Change `_count_scale()` from `alt.Scale(type="symlog", constant=1)` to `alt.Scale(type="log")` in `visualizer.py:274-276`. All 20 callsites inherit the fix automatically.

### Pattern 2 — BAM chart redundant legend (2 lines → 2 functions)
- Add `legend=None` to color encoding in `plot_bam_coverage_violin` and `plot_bam_dp_distribution` where x-axis and color both encode `metric:N`.

### Pattern 3 — Chart 43 partition aggregation bug (5 lines → 1 function)
- In `plot_filter_vaf_dp_heatmap`, aggregate data by summing across the `partition` field (train/val/test) before charting. Currently 3 text labels render at each cell position.

### Pattern 4 — Heatmap cell sizing + number abbreviation (5 lines → 1-3 functions)
- Increase cell size in `plot_fp_cross_tab_heatmap` from 20×20px to ≥35×35px.
- Use abbreviated number format (`.2s` → "94M") for large counts in cramped heatmaps.

### Pattern 5 — Sample-wise charts missing set_number (10 lines → data prep + 5 charts)
- Join `set_number` from sample manifest to sample-level aggregation data before passing to sample-wise chart functions in the `sample/` dimension.
- Add `.facet(facet="set_number:N", columns=2)` to the 5 sample-wise charts.

### Pattern 6 — Chart 22 axis swap + log scale (5 lines → 1 function)
- In `plot_per_sample_tier_distribution`, swap x/y: put `sample_id` on x-axis (rotated labels) and `n_variants` on y-axis with `_count_scale()` (log).

### Pattern 7 — Threshold sweep vertical → 3-column grid (10 lines → 2 functions)
- In `plot_vaf_threshold_sweep`, replace `row=alt.Row("classification:N")` with `.facet(facet="classification:N", columns=3)`.
- Apply similar layout to `plot_dp_threshold_sweep`.

### Bonus fixes from prior unmet requirements
- Pie chart: filter `SOMATIC_MODALITY_DOMAIN` to data-present categories (remove phantom "Unknown" legend entry).
- GT concordance: add `.fill_null("Unknown")` before faceting (fix "undefined" facet title).

## Capabilities

### New Capabilities
_None — all changes modify existing capabilities._

### Modified Capabilities
- `variant-visualization`: Log scale for all count axes, chart 22 axis swap, threshold sweep grid layout, BAM legend cleanup, sample-wise set faceting, pie/GT concordance fixes.
- `visualization-helpers`: `_count_scale()` return type change from symlog to log.
- `rescue-analytics`: Log scale propagation to 5 rescue chart count axes.

## Impact

- `bin/vcf_stats/seq2neo/visualizer.py` — ~60 lines changed across ~15 functions (most via 1-line `_count_scale()` fix)
- `bin/vcf_stats/seq2neo/cli.py` — ~10 lines for sample-wise set_number data join
- `bin/vcf_stats/tests/test_seq2neo_stats.py` — Update `_count_scale` test, add partition aggregation test
- Total: ~90 lines across 3 files
