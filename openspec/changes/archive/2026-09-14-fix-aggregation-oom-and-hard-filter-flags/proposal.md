## Why

The seq2neo stats pipeline OOM-kills at the cross-sample aggregation step when processing all 66 samples (~7.9M variants after unified filtering). The `combined_df.collect()` at line 1223 materializes all 165 columns into a single DataFrame, consuming 219 GB RSS before the OOM killer intervenes. Additionally, the hard filter silently drops ~3.1M variants (28%) during statistics, conflating observation with exclusion — statistics should characterize, not filter; exclusion belongs in downstream dataset preparation.

## What Changes

- **Replace eager `.collect()` (all 165 columns) with column-pruned collect (~15 columns)** at the confidence tier computation point, keeping `combined_df` as a true `LazyFrame` throughout the aggregation pipeline
- **Remove the wasteful re-scan block** (lines 1251–1275) that re-scanned from parquet, re-applied all filters, and wrapped the in-memory DataFrame as `.lazy()`
- **Convert hard filter from variant-dropping to flag-column generation**: add `flag_hard_*` boolean columns, `hard_filter_flags` (comma-joined string), and `n_hard_flags` (count) to the lazy scan — no variants are excluded during statistics
- **Extract hard filter conditions into a Python config module** (`hard_filter_config.py`) stored alongside the code, matching the `tier_config.py` pattern
- **Add hard filter breakdown statistics**: count variants matching each condition, write `hard_filter_breakdown.tsv`, and add a breakdown bar chart to the visualizer
- **Remove `--no-hard-filter` CLI flag** — no longer applicable since the hard filter no longer drops variants

## Capabilities

### New Capabilities
- `hard-filter-flag-config`: Hard filter conditions defined as a Python config module (`bin/vcf_stats/seq2neo/hard_filter_config.py`) with name, description, severity, required columns, and polars expression builders. Loaded at runtime, version-controlled with the code.
- `hard-filter-statistics`: Counting and visualization of hard filter condition breakdowns — per-condition variant counts written to `hard_filter_breakdown.tsv` and rendered as a severity-colored bar chart.

### Modified Capabilities
- `streaming-cross-sample-aggregation`: Fix the OOM-causing `.collect()` at the confidence tier computation to follow the column-pruned lazy pattern already specified. The confidence tier + soft flags computation collects only ~15 columns instead of all 165, then joins results back to the lazy scan. Remove the re-scan block that duplicated work.
- `variant-tiering-stats`: Hard filter no longer drops variants before confidence tiering. Confidence tiers are now assigned to ALL variants that pass the unified filter, with hard filter conditions recorded as flag columns for later use.

## Impact

- **`bin/vcf_stats/seq2neo/cli.py`**: Lines 1206–1275 replaced — flag columns instead of drop, pruned collect, join-back, no re-scan
- **`bin/vcf_stats/seq2neo/statistics.py`**: `build_hard_filter_expr()` removed; `build_hard_filter_flag_exprs()` and `hard_filter_breakdown()` added; `_CROSS_SAMPLE_COLS` updated with new flag columns
- **`bin/vcf_stats/seq2neo/hard_filter_config.py`**: New file — conditions as Python data
- **`bin/vcf_stats/seq2neo/visualizer.py`**: New hard filter breakdown bar chart
- **`memory-regression-tests`**: Peak memory spec validated against new flow (now safely under 4 GB for 66 samples)
- **Removed**: `--no-hard-filter` CLI argument (obsolete since filter no longer drops)
