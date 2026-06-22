## 1. Strelka TAR/TIR Field Inversion Fix (P0)

- [x] 1.1 Fix `caller_parser.py:447`: Change `("AD_REF", "TOR"), ("AD_ALT", "TAR")` to `("AD_REF", "TAR"), ("AD_ALT", "TIR")`
- [x] 1.2 Fix `cli.py:181`: Same mapping change as 1.1
- [x] 1.3 Fix `statistics.py:245` docstring: Update comment to say AD_ALT comes from TIR, AD_REF from TAR
- [x] 1.4 Add defensive `.clip(0.0, 1.0)` in `compute_vaf_columns()` after VAF = AD_ALT / DP for all callers
- [x] 1.5 Add `_repair_strelka_columns(df)` function in cli.py for post-load repair during `--resume`:
  - Detect: check if `AD_ALT == TAR` for non-null rows → old parquet → repair
  - Remap: `AD_ALT ← TIR`, `AD_REF ← TAR`
  - Recompute: `strelka_VAF = TIR / DP`, clipped to [0, 1]
  - Recompute: `DNA_VAF_mean`, `RNA_VAF_mean` (mean_horizontal of per-caller VAFs)
  - Recompute: `DNA_ALT_DP_mean`, `RNA_ALT_DP_mean`, `DNA_REF_DP_mean`, `RNA_REF_DP_mean`
- [x] 1.6 Wire `_repair_strelka_columns()` into the `--resume` path after parquet load, before statistics
- [x] 1.7 Add test: verify `AD_ALT == TIR` and `AD_REF == TAR` after Strelka join

## 2. DP Data Clipping (P1)

- [x] 2.1 `plot_dna_vs_rna_dp`: Add `.clip(0, 2000)` on `DNA_DP_mean` and `RNA_DP_mean` columns before plotting. Add outlier count subtitle.
- [x] 2.2 `plot_dp_boxplot_per_tier`: Add `.clip(0, 2000)` on DP column. Add outlier count subtitle.
- [x] 2.3 `plot_ref_alt_dp_scatter`: Add `.clip(0, 2000)` on both REF_DP and ALT_DP columns. Add outlier count subtitle.
- [x] 2.4 `plot_per_tier_dp_boxplot`: Add `.clip(0, 2000)` on DP column. Add outlier count subtitle.
- [x] 2.5 `plot_dp_distribution`: Add `.clip(0, 2000)` on DP column (dead code — will be wired in task 10.1)

## 3. DP Threshold Sweep Filter (P1)

- [x] 3.1 In `plot_dp_threshold_sweep()`: Add `classification.is_null()` filter before `to_pandas()`, mirroring `plot_vaf_threshold_sweep:1333`

## 4. VAF Violin Restore in Faceted Mode (P1)

- [x] 4.1 In `plot_vaf_distribution` faceted branch (line 649+): Build `transform_density` violin layer alongside boxplot
- [x] 4.2 Combine violin + boxplot with `alt.layer()` before applying column faceting
- [x] 4.3 Add explicit `width=alt.Step(60)` for per-caller spacing in both faceted and non-faceted branches

## 5. Sample Charts 2×2 Set Grid (P1)

- [x] 5.1 `plot_per_sample_distribution` (chart 08): Replace `enc["row"] = alt.Row("set_number:N")` with `.facet(facet='set_number:N', columns=2).resolve_scale(x='shared')`. Set per-panel `height=300`. Remove dynamic `n_samples * 12` height. Set sample as x-axis, counts as y-axis.
- [x] 5.2 `plot_bam_metrics_bars` (chart 20): Replace `enc["row"]` with `.facet(facet='set_number:N', columns=2).resolve_scale(y='shared')`. Set per-panel `width=350`.
- [x] 5.3 `plot_per_sample_tier_distribution` (chart 22): Replace `enc["row"]` with `.facet(facet='set_number:N', columns=2).resolve_scale(x='shared')`. Set per-panel `height=300`.

## 6. Caller Concordance Column Wrap (P1)

- [x] 6.1 In `plot_caller_concordance_vs_vaf`: Replace `enc["column"] = alt.Column(f"{color_col}:N")` with `.facet(facet=f"{color_col}:N", columns=3)`. Add `width=200` per panel.

## 7. Heatmap Text Labels (P1)

- [x] 7.1 `plot_filter_vaf_dp_heatmap` (chart 43): Add `mark_text(baseline="middle", fontSize=9)` layer with `text=alt.Text("count:Q", format=",d")` and conditional white/black color. Layer rect+text BEFORE `.facet()`.
- [x] 7.2 `plot_fp_cross_tab_heatmap` (chart 45): Same text overlay pattern as 7.1.

## 8. Pie Chart Labels (P1)

- [x] 8.1 `plot_somatic_modality_pie` (chart 46): Pre-compute percentage and label string (`"count\n(pct%)"`) in pandas. Add `mark_text(size=10, radiusOffset=20)` layer with `theta` + `text` encoding.

## 9. Rescue Analytics — Statistics (P1)

- [x] 9.1 Add `compute_rescue_breakdown(df, group_col="set_number")` → per-group rescued/non-rescued counts + proportions
- [x] 9.2 Add `compute_rescue_by_filter(df)` → RESCUED × FILTER cross-tab with counts and proportions
- [x] 9.3 Add `compute_rescue_cross_tab(df)` → RESCUED × FILTER × set_number 3-way cross-tab
- [x] 9.4 Add `compute_rescue_vaf_dp(df)` → mean/median/Q1/Q3 of DNA_VAF_mean, RNA_VAF_mean, DNA_DP_mean, RNA_DP_mean grouped by RESCUED
- [x] 9.5 Add `sample_rescue_summary(df)` → per-sample rescued/non-rescued counts with FILTER breakdown
- [x] 9.6 Add `compute_rescue_by_tier(df)` → rescue count and rate per final_tier (CxDy)
- [x] 9.7 Add `compute_rescue_by_caller_support(df)` → N_SUPPORT_CALLERS distribution by RESCUED status

## 10. Rescue Analytics — Charts (P1)

- [x] 10.1 `plot_rescue_breakdown` (chart 48) — stacked bar per set: rescued/non-rescued counts with % labels
- [x] 10.2 `plot_rescue_by_filter` (chart 49) — grouped bar: rescued/non-rescued per FILTER category
- [x] 10.3 `plot_rescue_cross_tab_heatmap` (chart 50) — heatmap: FILTER × set_number × rescued, with text labels
- [x] 10.4 `plot_rescue_vaf_boxplot` (chart 51) — violin+box: DNA/RNA VAF distributions for rescued vs non-rescued
- [x] 10.5 `plot_rescue_dp_boxplot` (chart 52) — boxplot: DNA/RNA DP distributions for rescued vs non-rescued
- [x] 10.6 `plot_rescue_sample_distribution` (chart 53) — per-sample rescue counts, faceted by set using 2×2 grid
- [x] 10.7 `plot_rescue_by_tier` (chart 54) — stacked bar: rescue rate per CxDy tier
- [x] 10.8 `plot_rescue_rate_trend` (chart 55) — line chart: rescue proportion across sets 1-4
- [x] 10.9 `plot_rescue_caller_support` (chart 56) — grouped bar: N_SUPPORT_CALLERS distribution by rescue status

## 11. CLI Wiring (P1)

- [x] 11.1 Add rescue section in cli.py after threshold analysis: compute all rescue stats, write TSVs to `stats/rescue/`, generate charts to `plots/rescue/`
- [x] 11.2 Add "Rescue Analysis" section to `_chart_section()` in visualizer.py for dashboard routing
- [x] 11.3 Wire `plot_dp_distribution` into the chart registry (currently dead code) or remove import

## 12. Cleanup and Tests

- [x] 12.1 Add test: Strelka mapping correctness (`AD_ALT == TIR`, `AD_REF == TAR`)
- [x] 12.2 Add test: `_repair_strelka_columns` detects and fixes old parquet
- [x] 12.3 Add test: `_repair_strelka_columns` skips already-fixed parquet (idempotent)
- [x] 12.4 Add test: rescue statistics functions return expected schema
- [x] 12.5 Add test: heatmap chart functions include text marks in output

## 13. Verification

- [ ] 13.1 Run pipeline with `--resume` on existing parquet, verify Strelka repair activates and VAFs are corrected
- [ ] 13.2 Verify Strelka VAFs are now in [0, 1] range (no values > 1.0)
- [ ] 13.3 Verify all DP charts have data clipped — no overflow beyond 2000
- [ ] 13.4 Verify DP sweep chart shows clean monotonic curves (no anomalous increases)
- [ ] 13.5 Verify VAF distribution shows violin overlay in faceted mode
- [ ] 13.6 Verify sample charts 08, 20, 22 render as 2×2 grid with shared value axis
- [ ] 13.7 Verify caller concordance chart wraps into 3-column grid
- [ ] 13.8 Verify heatmaps 43, 45 show count numbers on cells
- [ ] 13.9 Verify pie chart 46 shows count + percentage labels
- [ ] 13.10 Verify rescue TSVs generated in stats/rescue/
- [ ] 13.11 Verify rescue charts 48-56 generated in plots/rescue/
- [ ] 13.12 Run test suite, confirm no regressions
