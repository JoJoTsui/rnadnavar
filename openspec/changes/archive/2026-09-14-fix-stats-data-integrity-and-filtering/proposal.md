## Why

The seq2neo statistics and visualization pipeline has critical data integrity issues that corrupt all per-disease and per-set outputs, cause cross-run sample contamination, and leave BAM coverage charts empty in whole-genome mode. The existing `fix-broken-visualizations` change resolved rendering-level bugs (dashboard HTML, axis sorting, faceting toggle) but did not address the underlying data quality problems. Fixing these data integrity and filtering issues is a prerequisite for producing correct statistics and visualizations for downstream model training.

## What Changes

- **Fix `combined_df` sample contamination (CRITICAL)**: Add a sample_id filter to the lazy parquet scan so cumulative parquet files from prior `--set` runs don't contaminate current-run statistics and charts
- **Fix disease normalization (CRITICAL)**: Add "Unknown" sentinel in `process_single_sample` when disease values are empty string; fix upstream `parse_projects_to_json.py` to not default disease to `""`
- **Fix set_number validation (CRITICAL)**: Add warnings when all `set_number` values are 0; validate `partition_set` is present in manifest building
- **Fix `--resume` path filtering (HIGH)**: Filter `sample_summary.tsv` and `bam_stats.tsv` reloads to only the current run's sample_ids when `--set`/`--sample-ids` is specified
- **Fix WG coverage bins (HIGH)**: Extend Rust `whole_genome_stats_impl` to compute `cov_*_pct` for whole-genome mode via streaming histogram, so `coverage_distribution` chart renders with data in both WG and WES modes
- **Fix disease-wise chart registry (MEDIUM)**: Add `plot_caller_overlap` to the disease-wise registry to match the set-wise registry (16→17 charts)
- **Fix chart faceting coverage (MEDIUM)**: Wire `_apply_faceting()` into `plot_filter_distribution`, `plot_redi_evidence`, `plot_ref_alt_dp_scatter`, and `plot_dna_vs_rna_dp` so these charts facet by group dimension in set/disease wise directories
- **Fix `caller_tier` schema verification (MEDIUM)**: Verify `caller_tier` column exists in `combined_df` schema after parquet scan; recompute tiers if missing or warn
- **Fix BAM chart layouts (MEDIUM)**: Constrain `metrics_sample_wise` to fixed-width 3-column grid; pass `top_n` parameter to `plot_bam_metrics_bars` from CLI
- **Add confidence analysis pipeline (FEATURE)**: Add `confidence` wise dimension, `stats/confidence/` directory with TSV summaries, confidence-wise chart registry, `--min-confidence-tier` filter, and `--export-high-confidence` flag for ML training data export of HIGH-confidence Somatic/Germline/Reference/DeepSomatic variants

## Capabilities

### New Capabilities

- `confidence-analysis`: Confidence tier statistics and visualizations, `--min-confidence-tier` unified filter, `--export-high-confidence` ML data export. Includes `stats/confidence/` directory with per-confidence summaries, cross-tabulations with FILTER and tier, and filtered parquet output of HIGH-confidence variants suitable for model training.

### Modified Capabilities

- `variant-visualization`: Per-wise chart functions (`plot_filter_distribution`, `plot_redi_evidence`, `plot_ref_alt_dp_scatter`, `plot_dna_vs_rna_dp`) now call `_apply_faceting()` for consistent per-group faceting. `plot_bam_metrics_sample_wise` uses fixed-width 3-column grid. `plot_bam_metrics_bars` receives `top_n` from the CLI call site. Disease-wise registry gains `plot_caller_overlap`. `coverage_distribution` renders in both WG and WES modes. `caller_tier` schema is verified before chart generation.
- `variant-tiering-stats`: `FINAL_TIER_ORDER` is exported and consumed for confidence-tier x tier cross-tabulation. Confidence tier is added as a wise dimension alongside existing tier/disease/set dimensions.
- `bam-statistics`: `cov_*_pct` columns (1x, 10x, 20x, 50x, 100x) are now populated for whole-genome mode via streaming coverage histogram in Rust, not just WES/BED mode.

## Impact

- `bin/vcf_stats/seq2neo/cli.py` — ~40 lines: combined_df sample_id filter, disease/set defaults with "Unknown"/warnings, --resume filtering, confidence wise dimension, --min-confidence-tier, --export-high-confidence
- `bin/vcf_stats/seq2neo/visualizer.py` — ~50 lines: metrics_sample_wise fixed grid, bam_metrics_bars top_n pass-through, faceting for 4 chart functions, caller_tier verification, confidence-wise chart registry
- `bin/vcf_stats/seq2neo/statistics.py` — ~15 lines: confidence summary function, confidence metrics in wise kernel
- `bin/vcf_stats/seq2neo/bam_stats.py` — ~10 lines: WG coverage bins path (delegates to updated Rust function)
- `bin/vcf_stats/seq2neo/stats_core/src/bam.rs` — ~60 lines: WG coverage bin histogram computation in `whole_genome_stats_impl`
- `examples/seq2neo/build_sample_manifest.py` — ~5 lines: disease/set validation warnings
- `examples/seq2neo/scripts/parse_projects_to_json.py` — ~5 lines: disease default to last seen or "Unknown" instead of `""`
