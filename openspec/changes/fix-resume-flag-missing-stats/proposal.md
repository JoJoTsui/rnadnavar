## Why

The `--resume` flag added in `566f6fc` skips variant processing and BAM stats, but fails to reload `all_stats` from the existing `sample_summary.csv` and `bam_stats_df` from `bam_stats.csv`. This causes downstream cross-sample CSVs (`disease_summary.csv`, `tier_summary.csv`, `dataset_summary.csv`, `sample_tier_summary.csv`) and per-sample/BAM dashboard charts to be silently skipped. Additionally, the variant count from existing parquet files ignores `--max-samples`/`--set`/`--sample-ids` filters, producing mismatched totals.

## What Changes

- Reload `all_stats` from `sample_summary.csv` when resuming — enables cross-sample CSVs and per-sample charts
- Reload `bam_stats_df` from `bam_stats.csv` when resuming — enables BAM charts
- Filter variant count to match active `--max-samples`/`--set`/`--sample-ids` filters

## Capabilities

- None — bugfix only.

## Impact

- **Python**: `cli.py` — 4 lines added to the `--resume` block
