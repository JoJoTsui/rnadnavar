## 1. Pre-Decomposition Caller VCF Parsing

- [x] 1.1 Add pre-norm VCF path discovery to `manifest_loader.py` — add `pre_norm_vcf` key to each caller config in `CALLER_CONFIGS`, using `base_output_dir/{dir_name}/variant_calling/{caller}/.../*.vcf.gz` pattern
- [x] 1.2 Implement `parse_pre_norm_multiallelic()` in `caller_parser.py` for Mutect2 — parse FORMAT/AD (full array), FORMAT/AF (per-allele), FORMAT/F1R2, FORMAT/F2R1, FORMAT/GT from records where ALT contains comma, return polars DataFrame keyed on (CHROM, POS, REF, ALT)
- [x] 1.3 Implement pre-norm parsing for DeepSomatic — parse FORMAT/AD (full array), FORMAT/VAF (per-allele), FORMAT/GT
- [x] 1.4 Add `--no-pre-norm` CLI flag to skip pre-norm parsing (for environments where pre-norm VCFs are unavailable)
- [x] 1.5 In `process_single_sample()`, call `parse_pre_norm_multiallelic()` for Mutect2 and DeepSomatic callers, build allele registry DataFrame, left-join with normalized parse result on (CHROM, POS, REF, ALT)
- [x] 1.6 Tests: `test_parse_pre_norm_mutect2_multiallelic`, `test_parse_pre_norm_deepsomatic`, `test_allele_registry_join`, `test_pre_norm_biallelic_passthrough`, `test_pre_norm_vcf_missing_graceful`

## 2. Multi-Allelic Site Analysis

- [x] 2.1 Add `compute_multi_allelic_metrics()` to `statistics.py` — group by (CHROM, POS), compute n_alleles_at_site, vaf_sum, allele_balance_ratio, category_conflict, classify into normalization_artifact/noise/true_multi_allelic
- [x] 2.2 When pre-norm allele registry data is available, enrich multi-allelic metrics with per-allele AF from pre-norm AD arrays, per-allele strand balance from F1R2/F2R1, and GT co-occurrence
- [x] 2.3 Call `compute_multi_allelic_metrics()` in `process_single_sample()` after caller joining and VAF computation
- [x] 2.4 Add multi-allelic columns to `_CROSS_SAMPLE_COLS` — n_alleles_at_site, vaf_sum, allele_balance_ratio, category_conflict, multiallelic_class
- [x] 2.5 Tests: `test_multi_allelic_detection`, `test_normalization_artifact_classification`, `test_true_multi_allelic_classification`, `test_vaf_sum_computation`, `test_category_conflict_detection`, `test_pre_norm_enrichment`

## 3. modality_evidence Three-Way Classification

- [x] 3.1 Implement `compute_modality_evidence()` in `statistics.py` — compute `modality_evidence_caller` (based on N_DNA_CALLERS_SUPPORT/N_RNA_CALLERS_SUPPORT → C1-C4 mapping) and `modality_evidence_dp` (based on DNA_DP_mean/RNA_DP_mean > 30)
- [x] 3.2 Call `compute_modality_evidence()` in `process_single_sample()` after `compute_tiers_for_dataframe()`
- [x] 3.3 Compute deprecated RESCUED and CROSS_MODALITY columns from modality_evidence: RESCUED="YES" when modality_evidence_caller ∈ {cross_modality, rna_rescued}
- [x] 3.4 Add `modality_evidence_caller` and `modality_evidence_dp` to `_CROSS_SAMPLE_COLS`; remove RESCUED and CROSS_MODALITY from `_CROSS_SAMPLE_COLS` (retain in parquet only)
- [x] 3.5 Emit deprecation warning once per run when RESCUED or CROSS_MODALITY columns are accessed from cross-sample data
- [x] 3.6 Tests: `test_modality_evidence_caller_c1`, `test_modality_evidence_caller_c3_rna_rescued`, `test_modality_evidence_caller_c7_low_confidence`, `test_modality_evidence_dp_cross_modality`, `test_modality_evidence_dp_rna_rescued`, `test_rescued_column_maps_from_modality_evidence`, `test_cross_modality_column_maps_from_modality_evidence`

## 4. Biological Evidence Flags

- [x] 4.1 Implement `compute_biological_flags()` in `statistics.py` — compute flag_vaf_overflow, flag_multi_allelic_heterogeneity, flag_category_conflict, flag_germline_low_vaf, flag_somatic_high_vaf, flag_reference_with_signal, flag_rna_rescued
- [x] 4.2 Call `compute_biological_flags()` in `process_single_sample()` after multi-allelic metrics and modality_evidence
- [x] 4.3 Add all 8 flag columns to `_CROSS_SAMPLE_COLS`
- [x] 4.4 Tests: `test_vaf_overflow_flag`, `test_heterogeneity_flag`, `test_germline_low_vaf_flag`, `test_somatic_high_vaf_flag`, `test_reference_with_signal_flag`, `test_rna_rescued_flag`, `test_no_false_positive_flags_on_normal_variants`

## 5. RNA VAF Correction in Filtering

- [x] 5.1 Modify `--min-vaf` filter logic — remove `RNA_VAF_mean >= threshold` condition, replace with `(N_RNA_CALLERS_SUPPORT >= 2) & (RNA_DP_mean >= 10)`
- [x] 5.2 Update `--min-vaf` help text to document that RNA evidence uses caller support + depth, not RNA VAF (due to allele-specific expression)
- [x] 5.3 Tests: `test_min_vaf_uses_rna_caller_support_not_vaf`, `test_min_vaf_rna_leg_requires_two_callers`, `test_min_vaf_rna_leg_requires_min_depth`, `test_min_vaf_fallback_when_no_rna_columns`

## 6. Unified Filter Pipeline

- [x] 6.1 Implement `build_unified_filter(args, combined_df)` in `cli.py` — compose filter expression from: include/exclude FILTER categories, min DNA/RNA callers, min VAF (corrected), min DP, min evidence tier, max gnomAD AF, multi-allelic conflict exclusion, disease exclusion, VAF physics constraint
- [x] 6.2 Apply unified filter at combined_df creation point (after `pl.scan_parquet()`, before any statistics function)
- [x] 6.3 Add new CLI flags: `--include-filters`, `--exclude-filters`, `--min-dna-callers`, `--min-rna-callers`, `--min-evidence-tier`, `--max-gnomad-af`, `--exclude-multiallelic-conflict`
- [x] 6.4 Add `--no-filter` escape hatch (disables all filtering, including VAF physics constraint)
- [x] 6.5 Deprecate `--pileup-mode` — map "filtered" to `--exclude-filters NoConsensus`, emit deprecation warning
- [x] 6.6 Remove separate `variant_details_filtered/` output path — unified filter covers the use case; the existing `--min-vaf`, `--min-dp`, `--exclude-disease` flags now affect ALL output
- [x] 6.7 Tests: `test_unified_filter_include_filters`, `test_unified_filter_multiple_conditions`, `test_unified_filter_vaf_physics_constraint`, `test_no_filter_escape_hatch`, `test_pileup_mode_deprecation_warning`, `test_filter_propagates_to_statistics`, `test_filter_propagates_to_viz`

## 7. Statistics Migration

- [x] 7.1 Update `_WISE_METRICS` in `statistics.py` — add n_cross_modality, n_dna_confident, n_rna_rescued, n_low_confidence (counts per modality_evidence category); keep deprecated n_rescued and n_cross_modality as aliases
- [x] 7.2 Update rescue statistics functions (`compute_rescue_breakdown`, `compute_rescue_by_filter`, `compute_rescue_cross_tab`, `compute_rescue_vaf_dp`, `sample_rescue_summary`, `compute_rescue_by_tier`, `compute_rescue_by_caller_support`) to use `modality_evidence_caller` as primary grouping column
- [x] 7.3 Update `sample_summary()` and `dataset_summary()` to include new modality_evidence and biological flag counts
- [x] 7.4 Tests: `test_wise_metrics_modality_evidence`, `test_rescue_breakdown_with_modality_evidence`, `test_sample_summary_with_new_columns`, `test_backward_compat_rescued_count`

## 8. Visualization Updates (minimal — full viz in Proposal 2)

- [x] 8.1 Update `plot_cross_modality` signature to accept optional `modality_evidence` column; add fallback that maps RESCUED+CROSS_MODALITY to modality_evidence when new column isn't available
- [x] 8.2 Update rescue chart functions (`plot_rescue_*`) in `visualizer.py` to accept `evidence_col` parameter defaulting to `"modality_evidence_caller"` with fallback to `"RESCUED"`
- [x] 8.3 Register `plot_dna_vs_rna_per_caller` in variant-category wise chart list in `cli.py`
- [x] 8.4 Update `_RESCUE_COLOR_REGISTRY` to `_MODALITY_EVIDENCE_COLOR_REGISTRY` with 4 colors
- [x] 8.5 Tests: `test_cross_modality_migrates_to_modality_evidence`, `test_rescue_charts_accept_new_column`, `test_variant_category_wise_has_per_caller_chart`

## 9. Integration & Verification

- [x] 9.1 Run full pipeline on a single sample with `--no-filter` — verify parity with existing output (new columns present, old columns unchanged)
- [x] 9.2 Run full pipeline with `--include-filters Somatic Germline` — verify only those FILTER categories in output statistics
- [x] 9.3 Run full pipeline with `--min-vaf 0.05` — verify RNA caller support logic (not RNA VAF) is used, and that biologically rescued variants survive
- [x] 9.4 Run `--resume` — verify new columns are computed and added to existing parquet files without re-processing
- [x] 9.5 Run all 240 existing tests — verify no regressions; update any tests that assert exact column lists or deprecated metric names
- [x] 9.6 Run new tests (from sections 1-8) — verify all pass
