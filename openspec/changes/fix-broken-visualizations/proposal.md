## Why

The seq2neo statistics dashboard (`dashboard.html`) is broken — raw JavaScript fragments from the vega-embed library are rendered as visible text in every chart, the file is 100MB due to the library being inlined 106 times, and multiple per-wise plots (disease, set, tier, chromosome, BAM) have incorrect faceting, sorting, or rendering. These issues make the visualization output unusable for QC review.

## What Changes

- **Fix dashboard generation**: Stop inlining the full vega-embed JS library per chart; load it once via CDN `<script>` tag and embed only the vega-lite JSON spec per chart. Fix `_extract_body_content` to match the real `<body>` tag, not `<body>` strings inside the inlined JS library. Expected size reduction: ~100MB → <1MB.
- **Fix `coverage_distribution` plot**: Add a null-value guard so the chart returns early (or shows a "no data" message) when all `cov_*_pct` columns are null, instead of rendering an empty chart with null bars. Fix tooltip type from `nominal` to `quantitative`.
- **Fix `metrics_sample_wise` plot**: Replace `column=alt.Column("metric:N")` (single horizontal row) with `alt.Facet("metric:N", columns=3)` for a wrapped grid layout. Add `.facet(facet=alt.Facet("set_number:N"), columns=2)` so the chart is faceted by set. This aligns with all other sample-as-axis plots.
- **Fix `20_bam_metrics` bar width**: Replace fixed `width=350` with `width=alt.Step(20)` so bars scale with sample count instead of being squeezed into a fixed pixel width.
- **Fix `10_ti_tv_ratio` chromosome order**: Add `sort=chrom_order` to the text layer's x-encoding (the bar layer already has it). This is the only chromosome plot missing the sort on its text layer.
- **Fix disease-wise and set-wise plot faceting**: Pass `facet_col="disease_normalized"` and `facet_col="set_number"` in the wise chart registries so all chart functions facet by the wise dimension. Currently only 3 of ~17 chart functions internally call `_apply_faceting()`; the rest require explicit `facet_col`.
- **Fix tier-wise plot sorting**: Import `FINAL_TIER_ORDER` from `tiering_stats.py` into `visualizer.py` and apply `sort=FINAL_TIER_ORDER` to all tier x-axes. Currently all tier plots rely on alphabetical sort (works for C1-C9, breaks for C10+).
- **Fix non-reproducible sampling**: Add a fixed `seed=42` to all `df.sample()` calls in `_sample_if_large` so visualizations are deterministic across runs.
- **Fix orphaned plot functions**: Either wire up the 5 orphaned plot functions (charts 38-42: `plot_caller_tier_heatmap`, `plot_sample_overview_scatter`, `plot_mean_vaf_per_group`, `plot_mean_dp_per_group`, `plot_n_support_callers_dist`) into the chart registry, or remove them as dead code with stale docstrings.

## Capabilities

### New Capabilities

- `dashboard-generation`: Correct generation of the aggregate `dashboard.html` — single shared vega-embed runtime, correct body-content extraction, deterministic chart ordering.

### Modified Capabilities

- `variant-visualization`: Fix 8 broken/incorrect plot behaviors: dashboard body extraction, coverage null guard, metrics_sample_wise grid faceting, bam_metrics bar width, 10_ti_tv chromosome sort, disease/set wise faceting, tier-axis sorting, and sampling seed.
- `variant-tiering-stats`: Export `FINAL_TIER_ORDER` for consumption by the visualizer (currently defined but never imported by visualizer.py).

## Impact

- `bin/vcf_stats/seq2neo/visualizer.py` — ~120 lines: dashboard generation (~40 lines), `_extract_body_content` fix (~10 lines), `plot_bam_coverage_distribution` null guard (~10 lines), `plot_bam_metrics_sample_wise` faceting (~15 lines), `plot_bam_metrics_bars` width (~5 lines), `plot_ti_tv_ratio` text sort (~5 lines), `_sample_if_large` seed (~5 lines), tier sort imports + application (~20 lines), orphan cleanup (~10 lines)
- `bin/vcf_stats/seq2neo/cli.py` — ~20 lines: add `facet_col=` parameters to disease-wise and set-wise chart registry entries
- `bin/vcf_stats/seq2neo/tiering_stats.py` — ~0 lines (FINAL_TIER_ORDER already defined at line 67-68, just needs to be exported/used)
- Output: `dashboard.html` shrinks from ~100MB to <1MB; all per-wise plots render with correct faceting and sorting
