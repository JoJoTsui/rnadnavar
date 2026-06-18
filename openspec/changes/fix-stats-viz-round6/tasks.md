## 1. Scale Change — Raw Linear Counts (Approach C)

- [x] 1.1 Changed `_count_scale()` to return default linear `alt.Scale()` — raw counts, no log transform
- [x] 1.2 Changed `_count_axis()` to return `alt.Axis(format="~s")` — SI formatted labels (1K, 10K, 1M)
- [x] 1.3 Added `_make_bar_text()` helper for systematic count text labels on all bar charts

## 2. Text Labels on All 19 Bar Charts

### Simple bars (7 functions) — NEW count labels
- [x] 2.1 plot_chromosome_density — count label
- [x] 2.2 plot_gt_concordance (grouped branch) — count label
- [x] 2.3 plot_gt_concordance (ungrouped branch) — count label
- [x] 2.4 plot_gt_concordance_per_tier — count label with column faceting
- [x] 2.5 plot_per_sample_distribution — count label
- [x] 2.6 plot_low_vaf_rna_support — count label
- [x] 2.7 plot_somatic_modality_bars — count label

### Stacked bars (10 functions) — NEW or UPGRADED count+pct labels
- [x] 2.8 plot_vc_distribution — NEW count+pct label
- [x] 2.9 plot_variant_type_distribution — UPGRADED from pct-only to count+pct
- [x] 2.10 plot_cross_modality — UPGRADED from pct-only to count+pct
- [x] 2.11 plot_filter_distribution — UPGRADED from pct-only to count+pct
- [x] 2.12 plot_redi_evidence — NEW count+pct label
- [x] 2.13 plot_tiered_variant_types — NEW count label (simple, not pct — facet data compatibility)
- [x] 2.14 plot_per_sample_tier_distribution — NEW count+pct label
- [x] 2.15 plot_rescue_breakdown — UPGRADED from pct-only to count+pct
- [x] 2.16 plot_rescue_sample_distribution — NEW count+pct label
- [x] 2.17 plot_rescue_by_tier — UPGRADED from pct-only to count+pct

### Grouped bars (2 functions)
- [x] 2.18 plot_rescue_by_filter — NEW count+pct label
- [x] 2.19 plot_rescue_caller_support — UPGRADED to count+pct

## 3. Prior Round 6 Fixes (retained)

- [x] 3.1 Chart 43 partition aggregation bug fix
- [x] 3.2 Chart 22 axis swap (sample_id on x, n_variants on y)
- [x] 3.3 BAM chart legend=None
- [x] 3.4 Sample-wise set faceting (5 functions)
- [x] 3.5 Threshold sweep 3-column grid
- [x] 3.6 Heatmap cell sizing + SI abbreviation
- [x] 3.7 Pie chart zero-count domain filter
- [x] 3.8 GT concordance null facet fix
- [x] 3.9 Boxplot size=30 removed (auto-fit, 10 occurrences)
- [x] 3.10 Validation heatmap transposed + set-faceted

## 4. Tests

- [x] 4.1 test_count_scale_linear — passes (default linear scale)
- [x] 4.2 test_count_axis_si_format — passes (format="~s")
- [x] 4.3 test_make_bar_text_simple — passes
- [x] 4.4 test_make_bar_text_with_pct — passes
- [x] 4.5 All 39 relevant tests pass

## 5. Verification (requires pipeline run)

- [ ] 5.1 Run `--resume` and verify bar charts show raw counts with visible bars
- [ ] 5.2 Verify text labels on all count bar charts (count for simple, count+pct for stacked/grouped)
- [ ] 5.3 Verify DP boxplots no longer overflow
- [ ] 5.4 Verify validation heatmap is transposed with set faceting
