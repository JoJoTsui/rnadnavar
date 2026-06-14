## Why

The current pipeline accumulates all 65 samples' variant DataFrames in memory simultaneously, peaking at ~80 GB. At 12+ samples the process is OOM-killed by cgroup limits. This blocks scaling to the full 65-sample dataset. Additionally, visualization functions copy millions of rows to pandas and melt 78M-row DataFrames, join_caller_columns builds 1.2 GB of intermediate Python objects per sample, and aggregation functions redundantly scan the combined DataFrame 20+ times.

## What Changes

- **Streaming CLI architecture**: Each sample's variant_details is written to per-sample parquet files immediately after processing, the DataFrame is freed, and cross-sample aggregation uses `pl.scan_parquet()` for lazy evaluation. Peak memory drops from ~80 GB to ~500 MB.
- **Lazy aggregation functions**: All cross-sample statistics functions (`dataset_summary`, `disease_summary`, `tier_summary`, `sample_tier_summary`, etc.) accept lazy frames and materialize once.
- **Pre-aggregated visualization**: Chart functions receive pre-computed summary DataFrames instead of the 13M-row combined dataset. `.to_pandas()` calls removed — altair uses polars directly. VAF boxplots receive pre-computed quartiles instead of melted 78M-row DataFrames.
- **Memory-efficient join_caller_columns**: Use `pl.Series` with explicit dtypes directly from column vectors instead of row dicts, reducing intermediate Python object overhead.
- **Correctness fixes**: Fix broken `filter_distribution` sort, add warning on silent sample drops.

## Capabilities

### New Capabilities
- `streaming-variant-pipeline`: Per-sample parquet streaming with lazy cross-sample aggregation, enabling 65-sample processing within reasonable memory limits
- `memory-regression-tests`: Tests that verify DataFrames are freed after processing, lazy scan parity, and peak memory under thresholds

### Modified Capabilities
- `multi-level-aggregation`: Aggregation functions now accept lazy frames; scalar aggregate computation uses lazy evaluation from parquet scans
- `variant-visualization`: Chart functions accept pre-aggregated DataFrames; `.to_pandas()` calls replaced with polars-native altair

## Impact

- **Python code**: `cli.py` (major refactor), `statistics.py` (lazy support), `visualizer.py` (pre-aggregated inputs), `caller_parser.py` (Series construction)
- **Tests**: New `TestMemoryEfficiency` class, updates to existing test fixtures for lazy frames
- **Output**: Same CSV files and dashboard, identical format and values
- **Dependencies**: No new dependencies needed
