## 1. Data Integrity — combined_df Sample Filtering

- [x] 1.1 Add `current_sample_ids = [r["sample_id"] for r in rows]` after manifest filtering at cli.py:885
- [x] 1.2 Add `combined_df = combined_df.filter(pl.col("sample_id").is_in(current_sample_ids))` immediately after `pl.scan_parquet("*")` at cli.py:1103 (and again after re-scan at line 1174)
- [x] 1.3 Add log message: `f"combined_df filtered to {len(current_sample_ids)} sample(s) — {n_filtered} variants"` with count of variants pre/post filter
- [ ] 1.4 Verify: run `--set 1` twice, confirm second run does not include set 1 samples twice

## 2. Data Integrity — Disease Normalization

- [x] 2.1 In `parse_projects_to_json.py` line 98: change `disease = ""` to track last-seen disease from README lines; fall back to `"Unknown"` if no README line found in any block
- [x] 2.2 In `cli.py:process_single_sample` after line 557: add sentinel normalization — `if not disease_norm or disease_norm.strip() == "": disease_norm = "Unknown"` and `if not disease_val or disease_val.strip() == "": disease_val = "Unknown"`
- [x] 2.3 In `normalize_disease()` (common.py:18-20): add `if not result: return "unknown"` after strip/lower
- [ ] 2.4 Verify: check `variant_details/*.parquet` with `pl.scan_parquet` → assert no empty-string `disease_normalized` values

## 3. Data Integrity — set_number Validation

- [x] 3.1 After all samples processed, check if all `set_number` values are 0 or 1 unique value == 0: log warning `"WARNING: All samples have set_number=0 — manifest may be missing partition_set"`
- [x] 3.2 In `build_sample_manifest.py`: add validation that `partition_set` key exists in merged.json sample dicts; warn if any sample is missing it
- [x] 3.3 In `cli.py:process_single_sample` at line 553: keep default `0` but add `if set_num == 0: print(f"  WARNING: [{sample_id}] set_number=0 — manifest may be missing partition_set")`

## 4. Resume Path Filtering

- [x] 4.1 After reloading `sample_summary.tsv` at cli.py:955, filter to current sample_ids when `--set`/`--sample-ids` is active
- [x] 4.2 After reloading `bam_stats.tsv` at cli.py:975, filter to current sample_ids when `--set`/`--sample-ids` is active
- [x] 4.3 Add log message showing how many rows were filtered out of each TSV
- [ ] 4.4 Verify: run with `--set 1` on existing output dir → only set 1 samples in reloaded TSV data

## 5. BAM Coverage Bins — Whole-Genome Mode

- [x] 5.1 Extend `BamStats` struct in `bam.rs` to include `cov_1x_pct` through `cov_100x_pct` as `Option<f64>` fields
- [x] 5.2 In `whole_genome_stats_impl`: initialize 5 depth threshold counters (`bases_1x` through `bases_100x`) as `u64`
- [x] 5.3 During the existing CIGAR walk (after alignment match check): increment per-base depth; when depth crosses each threshold, increment that threshold's counter once per base position
- [x] 5.4 At end of scan: compute percentages as `100.0 * bases_Nx / total_bases` where `total_bases` = sum of reference sequence lengths from BAM header
- [x] 5.5 Update `whole_genome_stats()` and `whole_genome_stats_bed()` return dicts in Python `_compute_bam_stats_rust()` to include `cov_*_pct` fields
- [x] 5.6 In `bam_stats.py:compute_bam_stats`: for WG mode, use `cov_*_pct` from the Rust result instead of setting them to `None`
- [x] 5.7 Update `_BAM_STATS_COLUMNS` in `bam_stats.py` if needed (columns already exist, just now populated)
- [ ] 5.8 Test: run without `--bed` flag, verify `bam_stats.tsv` has non-null `cov_1x_pct` through `cov_100x_pct` values
- [ ] 5.9 Test: run with `--bed` flag (WES), verify coverage bins unchanged from current behavior

## 6. Chart Faceting Coverage

- [x] 6.1 Add `chart = _apply_faceting(chart, group_col)` to `plot_filter_distribution` before `_save_chart` call
- [x] 6.2 Add `chart = _apply_faceting(chart, group_col)` to `plot_redi_evidence` before `_save_chart` call
- [x] 6.3 Add `chart = _apply_faceting(chart, group_col)` to `plot_ref_alt_dp_scatter` before `_save_chart` call (note: this function generates two charts — facet both)
- [x] 6.4 Add `chart = _apply_faceting(chart, group_col)` to `plot_dna_vs_rna_dp` before `_save_chart` call
- [ ] 6.5 Verify: set-wise and disease-wise directories now have faceted versions of these 4 chart types

## 7. Disease Registry Completeness

- [x] 7.1 Add `(plot_caller_overlap, {"group_col": "disease_normalized"})` to `_WISE_CHART_REGISTRY["disease"]` between `plot_vc_distribution` and `plot_vaf_distribution` entries (matching position in set registry)

## 8. caller_tier Schema Verification

- [x] 8.1 Before the tier chart generation loop (cli.py ~1691), check `combined_df.collect_schema().names()` for `caller_tier`, `database_tier`, `final_tier`
- [x] 8.2 If any tier column is missing, compute `combined_eager = combined_df.collect()` and run `compute_tiers_for_dataframe(combined_eager)`, then replace `combined_df = combined_eager.lazy()`
- [x] 8.3 Log warning: `"WARNING: Tier columns missing from parquet schema — recomputed tiers for {n} variants"`
- [x] 8.4 If tier computation fails (e.g., missing FILTERS_NORMALIZED column), log error and skip tier-dependent charts with individual try/except guards

## 9. BAM Chart Layout Fixes

- [x] 9.1 In `plot_bam_metrics_sample_wise`: change `width=alt.Step(15)` to `width=200`, ensure `columns=3` is explicit in the facet call
- [x] 9.2 In `cli.py:1709`: change `plot_bam_metrics_bars(bam_stats_df, str(output_dir))` to `plot_bam_metrics_bars(bam_stats_df, str(output_dir), top_n=30)`
- [x] 9.3 In `plot_bam_metrics_bars`: change default `top_n=20` to `top_n=30`
- [ ] 9.4 Verify: with 12 samples, `20_bam_metrics.html` has readable bar widths; with >30 samples, only top 30 shown

## 10. Confidence Analysis Pipeline

- [x] 10.1 Add `"confidence"` to `all_wise_names` list at cli.py:1582
- [x] 10.2 Add `("confidence", ["confidence_tier"])` to `wise_configs` list at cli.py:1317
- [x] 10.3 Add `_WISE_CHART_REGISTRY["confidence"]` entry with: `plot_vc_distribution` (group_col="confidence_tier"), `plot_variant_type_distribution` (group_col="confidence_tier"), `plot_vaf_distribution` (color_col="confidence_tier"), `plot_caller_overlap` (group_col="confidence_tier"), `plot_filter_distribution` (group_col="confidence_tier"), `plot_cosmic_gnomad_annotation` (group_col="confidence_tier"), `plot_cross_modality` (group_col="confidence_tier"), `plot_ti_tv_ratio` (group_col="confidence_tier")

## 11. Verification

- [ ] 11.1 Run with `--set 1` — verify only set 1 samples in all TSVs and charts (no cross-contamination)
- [ ] 11.2 Run with `--set 2` after set 1 — verify no accumulation from set 1 parquet files
- [ ] 11.3 Verify `disease_normalized` values are actual disease names (not empty string)
- [ ] 11.4 Verify `set_number` values are 1, 2, 3, 4 (not all 0) — or warning is logged
- [ ] 11.5 Verify `coverage_distribution.html` renders with grouped bars in WG mode (no `--bed`)
- [ ] 11.6 Verify `18_caller_overlap_per_tier.html` has bar panels for each caller tier
- [ ] 11.7 Verify `08_per_sample_distribution.html` has only current run's samples with correct set faceting
- [ ] 11.8 Verify `20_bam_metrics.html` has readable bars with at most 30 samples
- [ ] 11.9 Verify `dashboard.html` is under 2MB and renders in browser
- [ ] 11.10 Verify `stats/confidence/` directory exists with TSVs and charts
- [ ] 11.11 Verify `high_confidence_variants.parquet` exports with correct FILTER values
