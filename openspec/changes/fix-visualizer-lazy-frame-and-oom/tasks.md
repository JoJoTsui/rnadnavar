## 1. Fix `_save_chart` — vl-convert crash

- [x] 1.1 Wrap `chart.save(format="png")` in try/except `(ImportError, ModuleNotFoundError)`
- [x] 1.2 Wrap `chart.save(format="svg")` in try/except `(ImportError, ModuleNotFoundError)`
- [x] 1.3 Keep `chart.save(html_path)` unconditional (always works)

## 2. Fix `plot_gt_concordance_per_tier` row slicing bug

- [x] 2.1 Replace `row[:-1]` with `row[:len(existing_gt)]` for GT extraction
- [x] 2.2 Use explicit `row[n_gt]` for `caller_tier` instead of `row[-1]`
- [x] 2.3 When `facet_col` present, use `row[n_gt + 1]` for facet value

## 3. Harden helpers against ColumnNotFoundError

- [x] 3.1 `_sample_if_large`: wrap `.collect()` in try/except `pl.exceptions.ColumnNotFoundError`, return empty DataFrame
- [x] 3.2 `_maybe_collect`: same wrapping

## 4. Add `--exclude-sample-ids` CLI flag

- [x] 4.1 Add `--exclude-sample-ids` argument to argparse parser
- [x] 4.2 Add filter: `manifest = manifest.filter(~pl.col("sample_id").is_in(args.exclude_sample_ids))`
- [x] 4.3 Place filter after `--set`, `--sample-ids`, `--max-samples` filters

## 5. Optimize count-based charts (batch count queries)

- [x] 5.1 `plot_cosmic_gnomad_annotation`: pre-aggregate with `group_by().agg()` instead of per-group `_count_rows`
- [x] 5.2 `plot_database_enrichment_by_tier`: same optimization
- [x] 5.3 Verify results match pre-optimization output (test_cosmic_gnomad_batch_counts_match_per_group passes)

## 6. Remove dead code

- [x] 6.1 Remove `_chromosome_sort_key` function (unused)
- [x] 6.2 Remove `plot_per_sample_violin` alias if it exists — does NOT exist, nothing to remove

## 7. Tests

- [x] 7.1 `test_save_chart_survives_missing_vl_convert` — HTML saved, PNG/SVG skipped gracefully
- [x] 7.2 `test_gt_concordance_per_tier_facet_col_correct_gts` — `caller_tier` not included in GT list when `facet_col` present
- [x] 7.3 `test_gt_concordance_per_tier_without_facet_correct` — correct behavior without `facet_col`
- [x] 7.4 `test_sample_if_large_column_not_found_returns_empty` — returns empty DataFrame instead of crashing
- [x] 7.5 `test_maybe_collect_column_not_found_returns_empty` — same
- [x] 7.6 `test_cli_exclude_sample_ids` — `--exclude-sample-ids` correctly removes samples from manifest
- [x] 7.7 `test_cosmic_gnomad_batch_counts_match_per_group` — optimized batch counts match per-group results

## 8. Verification

- [x] 8.1 Run existing test suite → 0 regressions (188 passed, 1 pre-existing flaky GIL test failed)
- [ ] 8.2 Run pipeline with `--exclude-sample-ids PRJNA298376_4264` → all 64 good samples complete, all charts generate
- [ ] 8.3 Run pipeline with vl-convert not installed → HTML charts generate, PNG/SVG skipped with message
