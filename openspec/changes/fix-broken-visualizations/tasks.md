## 1. Dashboard Generation Fix

- [x] 1.1 Replace `inline=True` with `inline=False` in `generate_dashboard()` — CDN script tags added to head
- [x] 1.2 Fix `_extract_body_content()` to strip `<script>` blocks before searching for `<body>` tag
- [x] 1.3 Add single vega-embed/vega-lite/vega CDN `<script>` tags to the dashboard `<head>` template
- [ ] 1.4 Verify dashboard.html is under 2MB and renders correctly in a browser (requires full pipeline run)
- [x] 1.5 Add `SAMPLE_SEED = 42` module-level constant and pass `seed=SAMPLE_SEED` to all `df.sample()` calls

## 2. Coverage Distribution Null Guard

- [x] 2.1 Add null-value check in `plot_bam_coverage_distribution`: check `is_not_null().any()` for cov_*_pct columns, return early if all null
- [x] 2.2 Fix tooltip type from `"nominal"` to `"quantitative"` for `mean_pct`

## 3. BAM Metrics Plots Fix

- [x] 3.1 Replace `column=alt.Column("metric:N")` with `alt.Facet("metric:N", columns=3)` in `plot_bam_metrics_sample_wise`
- [x] 3.2 Add set_number faceting when multiple sets exist
- [x] 3.3 Replace fixed `width=350` with `width=alt.Step(20)` in `plot_bam_metrics_bars`

## 4. Chromosome Sort Fix

- [x] 4.1 Add `sort=chrom_order` to the text layer x-encoding in `plot_ti_tv_ratio`

## 5. Disease-wise and Set-wise Faceting Fix

- [ ] 5.1 Add `facet_col="disease_normalized"` to all disease-wise registry entries in `cli.py` (deferred — requires CLI registry update)
- [ ] 5.2 Add `facet_col="set_number"` to all set-wise registry entries in `cli.py` (deferred — requires CLI registry update)
- [ ] 5.3 Verify chart functions correctly facet (requires full pipeline run with multiple sets)

## 6. Tier-axis Sort Fix

- [x] 6.1 Import `FINAL_TIER_ORDER` from `tiering_stats.py` at the top of `visualizer.py`
- [x] 6.2 Apply `sort=FINAL_TIER_ORDER` to all `alt.X("final_tier:N")` and `alt.X("tier:N")` encodings
- [x] 6.3 Apply `sort=FINAL_TIER_ORDER` to all `alt.Column("caller_tier:N")` encodings

## 7. Orphaned Plot Functions

- [ ] 7.1 Decide: wire up or remove 5 orphaned functions (deferred — low priority)
- [ ] 7.2 If wiring up: fix bugs, add to registries (deferred)
- [ ] 7.3 If removing: delete functions, remove imports (deferred)

## 8. Verification

- [ ] 8.1 Run `nf-test test tests/default.nf.test --profile test,docker` (requires full pipeline)
- [ ] 8.2 Verify dashboard.html renders in a browser (requires full pipeline run)
- [ ] 8.3 Verify all per-wise plots facet correctly (requires full pipeline run)
- [ ] 8.4 Verify tier-axis plots show tiers in FINAL_TIER_ORDER (requires full pipeline run)
- [x] 8.5 Verify two consecutive runs produce identical sampled charts (seed test passed)
