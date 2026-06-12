## Why

The visualizer module (`visualizer.py`) has 17 crash bugs from calling `.height`, `.iter_rows()`, and subscript access on `pl.LazyFrame` objects — methods that only exist on eager `DataFrame`. These bugs were masked while the pipeline OOM-killed before reaching visualizations. Now that the cross-sample aggregation is fixed, every pipeline run will crash at the first chart (`plot_cosmic_gnomad_annotation` line 206: `df.height` on a LazyFrame).

Beyond the crashes, the visualizer reads all 58M+ rows from 24 parquet files for each of 23 charts — 23 sequential full scans. While column projection keeps memory manageable, several functions load 50K-70M rows before sampling, and two functions (`plot_gt_concordance`, `plot_gt_concordance_per_tier`) use `.iter_rows()` to process all 58M rows in Python.

## What Changes

- **Fix 17 lazy-frame crash bugs**: Replace `.height` on LazyFrame with `.select(pl.len()).collect().item()`. Replace `.iter_rows()` on LazyFrame with `.collect().iter_rows()`. Replace subscript + `.to_list()` on LazyFrame with `.select().collect()["col"].to_list()`.
- **Reorganize chart data loading**: Extract a shared `_count_rows(df)` helper and a `_sample_if_large(df, max_rows)` helper used by all chart functions. Standardize the pattern: select needed columns → count → sample if needed → collect → pandas.
- **Fix `plot_bam_coverage_violin`**: It accesses `BAM_DP_*` columns from `combined_df` but those columns don't exist in the per-sample parquets. Either remove the chart or source data from BAM validation instead.
- **Optimize expensive charts**: `plot_cosmic_gnomad_annotation` (3 count queries) and `plot_caller_agreement_matrix` (36 count queries) — use lazy `.select(pl.len()).collect().item()` instead of full collects.
- **Remove unused `_maybe_collect`**: All charts will use explicit `.collect()` instead of the pipe pattern, making the data flow explicit.

## Capabilities

### Modified Capabilities
- `variant-visualization`: All 23 chart functions updated for LazyFrame compatibility, OOM safety, and polars 1.41.2 API conformance.

## Impact

- **Python**: `visualizer.py` — ~25 functions updated, ~15 new helper lines, no API changes.
- **Memory**: No chart loads more than 50K rows into memory (sampled charts) or uses efficient count queries (aggregate charts). Peak ~2 GB down from potential ~70M-row unpivot.
- **Speed**: 23 sequential parquet scans remain (I/O bound), but each scan now uses column projection efficiently. Count queries don't load data.
