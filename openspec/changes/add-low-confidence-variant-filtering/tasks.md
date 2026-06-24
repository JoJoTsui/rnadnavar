## 1. New Metrics Computation

- [x] 1.1 Add `n_recurrent_samples` computation: `pl.len().over(['CHROM','POS','REF','ALT'])` on `combined_df` in cli.py
- [x] 1.2 Add `n_recurrent_samples` to `_CROSS_SAMPLE_COLS` in `statistics.py`
- [ ] 1.3 Add `max_strand_bias_fisher_p` computation (deferred — requires scipy Fisher exact; placeholder for now)
- [ ] 1.4 Add `max_strand_bias_fisher_p` to parquet schema and `_CROSS_SAMPLE_COLS` (deferred with 1.3)

## 2. New Biological Flags

- [x] 2.1 Add `flag_germline_high_vaf` in `compute_biological_flags`: FILTER=="Germline" AND DNA_VAF_mean > 0.85
- [x] 2.2 Add `flag_somatic_loh` in `compute_biological_flags`: FILTER=="Somatic" AND DNA_VAF_mean > 0.60 AND DNA_REF_DP_mean < 0.10 × DNA_DP_mean
- [x] 2.3 Add `flag_no_caller_support` in `compute_biological_flags`: N_SUPPORT_CALLERS == 0
- [x] 2.4 Add `flag_low_rna_mapq` in `compute_biological_flags`: BAM_RT_mean_MQ < 2 AND RNA_DP_mean >= 20
- [x] 2.5 Add all 4 new flags + category_conflict_resolution to `_CROSS_SAMPLE_COLS`

## 3. Stage 0 Hard Exclusions

- [x] 3.1 Implement `build_hard_filter_expr()` in `statistics.py` — returns boolean expression for all Stage 0 conditions
- [x] 3.2 Wire hard filter in cli.py with `--no-hard-filter` escape hatch
- [ ] 3.3 Write dropped variants to `failed.tsv` with drop_reason column (deferred — needs separate tracking)
- [x] 3.4 Add `--no-hard-filter` CLI flag

## 4. Stage 1 Soft Flags

- [x] 4.1 Implement `compute_soft_flags()` returning comma-separated string of flag names
- [x] 4.2 Add `soft_flags` column to `_CROSS_SAMPLE_COLS`

## 5. Stage 2 Confidence Tiers

- [x] 5.1 Implement `compute_confidence_tier()`: HIGH/MEDIUM/LOW based on final_tier + soft flags
- [x] 5.2 Add `confidence_tier` column to `_CROSS_SAMPLE_COLS`
- [x] 5.3 Generate `confidence_tier_summary.tsv` in cli.py
- [x] 5.4 Add `--confidence-tier` CLI flag

## 6. Both Filtered and Unfiltered Statistics

- [x] 6.1 sample_summary.tsv (unfiltered) is already computed from all_stats before filtering
- [ ] 6.2 Add sample_summary_filtered.tsv (computed from filtered combined_df) — deferred to integration
- [ ] 6.3 Same for set_summary_filtered — deferred
- [ ] 6.4 Same for dataset_summary_unfiltered — deferred
- [ ] 6.5 Same for rescue_validation_summary — deferred

## 7. Caller Overlap Outputs

- [x] 7.1 Implement `caller_overlap_matrix()` in `statistics.py`: N×N pairwise co-occurrence matrix
- [x] 7.2 Write `caller_overlap_matrix.tsv` in cli.py
- [x] 7.3 Wire `caller_support_distribution()` into CLI output → write `caller_support_distribution.tsv`
- [x] 7.4 Rename `caller_overlap.tsv` to `tier_distribution.tsv` in cli.py

## 8. Visualization

- [ ] 8.1 Register `confidence_tier` in `_COLOR_REGISTRY` (deferred to fix-broken-visualizations)
- [ ] 8.2 Add confidence_tier_distribution plot (deferred to fix-broken-visualizations)
- [ ] 8.3 Add cross_sample_recurrence plot (deferred to fix-broken-visualizations)

## 9. Tests

- [x] 9.1 Test: hard filter drops correct variants (N_SUPPORT=0, flag_vaf_overflow, noise, germline_low_vaf, ref_with_signal) — verified
- [x] 9.2 Test: soft flags mark correct variants (category_conflict, rna_rescued, low_confidence, high_recurrence) — verified
- [x] 9.3 Test: confidence tiers assign correctly (HIGH for C1D1 Somatic, LOW for soft-flagged) — verified
- [ ] 9.4 Test: caller_overlap_matrix is symmetric — needs caller VAF column test data
- [ ] 9.5 Run `nf-test test tests/default.nf.test --profile test,docker` (requires full pipeline)

## 10. Integration Verification

- [ ] 10.1 Full pipeline run verification (requires test dataset + pipeline execution)
- [ ] 10.2 Verify `--no-hard-filter` preserves all variants with soft flags only
- [ ] 10.3 Verify `--confidence-tier high` filters to HIGH-confidence variants only
