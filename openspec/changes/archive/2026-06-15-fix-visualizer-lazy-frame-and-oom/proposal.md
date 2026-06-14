## Why

The visualizer module has issues discovered during the first successful full-pipeline run (65 samples, June 2026):

1. **`_save_chart` crash**: PNG/SVG export fails if `vl-convert` is not installed, killing all visualization generation including HTML (which doesn't need vl-convert).

2. **`plot_gt_concordance_per_tier` bug**: When `facet_col` is present, `row[:-1]` incorrectly includes `caller_tier` in the GT agreement list, corrupting concordance counts. This bug is silently producing wrong data in tier-wise GT concordance charts.

3. **Performance**: Several chart functions (`plot_cosmic_gnomad_annotation`, `plot_caller_agreement_matrix`, `plot_database_enrichment_by_tier`, `plot_caller_tier_heatmap`) do O(n_groups) separate `_count_rows` calls, each triggering a full lazy scan of all parquet files. With 65 samples, this makes visualization generation extremely slow.

4. **Defensive hardening**: `_sample_if_large` and `_maybe_collect` have no `ColumnNotFoundError` handling. While the missing-column issue in this run was caused by truncated BAM files for one sample (PRJNA298376_4264), defensive handling prevents future edge cases from crashing the pipeline.

5. **Missing CLI flag**: No `--exclude-sample-ids` to skip broken samples without editing the manifest.

## What Changes

- **Fix `_save_chart`**: Wrap PNG/SVG saves in try/except `ImportError` — HTML always succeeds
- **Fix `plot_gt_concordance_per_tier`**: Use `row[:len(existing_gt)]` instead of `row[:-1]` to correctly exclude `caller_tier` and `facet_col` from GT list
- **Optimize count-based charts**: Batch count queries instead of per-group lazy scans
- **Harden helpers**: Add `ColumnNotFoundError` catch in `_sample_if_large` and `_maybe_collect`
- **Add `--exclude-sample-ids`**: CLI flag to skip specific samples
- **Remove dead code**: `_chromosome_sort_key` (unused), `plot_per_sample_violin` alias
- **Update tests**: 7 new test methods covering vl-convert crash, GT concordance bug, exclusion flag, and helper hardening

## Capabilities

### Modified Capabilities
- `variant-visualization`: All chart functions updated

## Impact

- **Python**: `visualizer.py` (~5 functions modified), `cli.py` (1 new flag, 1 import fix)
- **Memory**: No change — charts already use efficient patterns
- **Speed**: Count-based charts significantly faster after batching (eliminates redundant scans)
- **Robustness**: Visualizations survive missing vl-convert, missing columns, and broken-sample edge cases
