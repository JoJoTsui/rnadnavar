## 1. Hard filter config module

- [x] 1.1 Create `bin/vcf_stats/seq2neo/hard_filter_config.py` with `HARD_FILTER_CONDITIONS` list (8 conditions, each with `name`, `description`, `severity`, `flag_column`, `required_columns`, `build_expr`)
- [x] 1.2 Add utility functions: `build_hard_filter_flag_exprs(columns)`, `build_hard_filter_summary_exprs(columns)`, `load_hard_filter_config()`
- [x] 1.3 Verify the module is importable and each condition's `build_expr` returns a valid `pl.Expr` for its required columns

## 2. Statistics module changes

- [x] 2.1 Define `_CONFIDENCE_TIER_COLS` constant in `statistics.py` — the ~17 columns needed by `compute_confidence_tier` and `compute_soft_flags`
- [x] 2.2 Remove `build_hard_filter_expr()` function from `statistics.py`
- [x] 2.3 Add `build_hard_filter_flag_exprs()` wrapper that imports from `hard_filter_config` and returns list of `pl.Expr` for existing columns
- [x] 2.4 Add `hard_filter_breakdown(df)` function that takes a LazyFrame, counts variants per hard filter condition, and returns a DataFrame suitable for TSV writing and charting
- [x] 2.5 Update `_CROSS_SAMPLE_COLS` to include `flag_hard_*` (8 columns), `hard_filter_flags`, and `n_hard_flags`

## 3. CLI aggregation block rewrite (lines 1206–1275)

- [ ] 3.1 At line 1206: replace hard filter drop (`filter(~hard_expr)`) with `.with_columns()` adding all `flag_hard_*` columns + `hard_filter_flags` + `n_hard_flags` via `build_hard_filter_flag_exprs()` and `build_hard_filter_summary_exprs()`
- [ ] 3.2 Add hard filter breakdown stats call after flag columns: compute `hard_filter_breakdown()` and write `hard_filter_breakdown.tsv` to `stats/hard_filter/`
- [ ] 3.3 At line 1223: replace `combined_df.collect()` (all columns) with `combined_df.select(_CONFIDENCE_TIER_COLS).collect()` (pruned to ~17 columns)
- [ ] 3.4 After confidence tier computation: left-join `confidence_tier` and `soft_flags` back to the main `combined_df` lazy scan using `combined_df.join(ct_lookup, on=["sample_id", "CHROM", "POS"], how="left")`
- [ ] 3.5 Delete the re-scan block (lines 1251–1275): remove the second `pl.scan_parquet()`, the re-filtering, the re-hard-filter, and the `combined_eager.lazy()` wrap

## 4. Remove `--no-hard-filter` CLI flag

- [ ] 4.1 Remove `--no-hard-filter` argument definition from the argparse parser (~line 804)
- [ ] 4.2 Remove the `if not args.no_hard_filter:` guard at line 1207 (and the duplicate at line 1266, which is removed as part of the re-scan deletion)
- [ ] 4.3 Remove `--no-hard-filter` from any shell scripts or documentation that reference it

## 5. Hard filter breakdown visualization

- [ ] 5.1 Add `hard_filter_breakdown_chart()` function to `visualizer.py` that reads `hard_filter_breakdown.tsv` and generates a horizontal bar chart colored by severity
- [ ] 5.2 Wire the chart call into `cli.py` after the hard filter breakdown TSV is written

## 6. Validation and testing

- [ ] 6.1 Run the full pipeline with `--max-samples 5` and verify it completes without OOM, produces `hard_filter_breakdown.tsv` and the chart
- [ ] 6.2 Run with the full 66-sample manifest and verify peak RSS is under 5 GB
- [ ] 6.3 Verify output parity: `sample_summary.tsv`, `disease_summary.tsv`, `tier_summary.tsv`, `caller_overlap_matrix.tsv` match the pre-change baseline for a small test run
- [ ] 6.4 Verify `hard_filter_breakdown.tsv` contains all 8 conditions with non-negative counts
- [ ] 6.5 Verify `--resume` still works with the new flow (load existing parquets, skip variant extraction, compute hard filter flags lazily)
- [ ] 6.6 Run existing nf-test suite to catch regressions in the stats pipeline (`nf-test test --profile test,docker`)
