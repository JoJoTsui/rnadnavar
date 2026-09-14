## 1. Fix Multi-allelic Classification

- [x] 1.1 Replace the `when()`-chain classification in `compute_multi_allelic_metrics` with per-allele criteria: normalization_artifact (REF/ALT length difference OR VAF≤0.001 AND DP≤1), noise (exactly one allele with VAF≥0.01 AND DP≥5), true_multi_allelic (≥2 alleles with VAF≥0.01 AND DP≥5 AND different ALT bases)
- [x] 1.2 Handle null DNA_ALT_DP_mean in classification: treat null as VAF=0, DP=0 (fill_null before signal detection)
- [ ] 1.3 Verify classification distribution shifts from 0.03%/73%/27% to closer to 60%/35%/5% (requires full pipeline re-run)

## 2. Fix allele_balance_ratio Formula

- [x] 2.1 Change `allele_balance_ratio` computation from `alt_dp_2nd_max / alt_dp_max` to `max_alt_dp / total_alt_dp`
- [x] 2.2 Ensure the formula handles 3+ alleles (denominator includes all alleles)
- [ ] 2.3 When pre-norm per-allele AF is available, prefer `max(pre_norm_af) / sum(pre_norm_af)` (deferred — pre-norm wiring in task 4)

## 3. Fix Biological Flags

- [x] 3.1 Fix `flag_multi_allelic_heterogeneity`: change to `== "true_multi_allelic"` only (excludes normalization_artifact)
- [x] 3.2 Fix `flag_category_conflict`: add biological-significance guard — only flag when FILTER is in {Somatic, Germline, RNAedit}
- [x] 3.3 Add `category_conflict_resolution` column with values: keep_germline, keep_somatic, keep_rnaedit, drop_both
- [x] 3.4 Remove `DNA_DP_mean >= 10` constraint from `flag_germline_low_vaf`
- [ ] 3.5 Add `flag_extreme_strand_bias` flag (deferred — requires strand_balance from pre-norm wiring, task 4)

## 4. Wire Pre-decomposition Data

- [ ] 4.1 Add derived scalar columns (`gt_cooccurrence`, `strand_balance`) to `_CROSS_SAMPLE_COLS` (deferred to integration phase — requires pre-norm data flow verification)
- [ ] 4.2 Compute `gt_cooccurrence` from `pre_norm_*_gt_alleles` in `process_single_sample()`
- [ ] 4.3 Compute `strand_balance` as `F1R2_alt / (F1R2_alt + F2R1_alt)`
- [ ] 4.4 Use pre-norm per-allele AF in `vaf_sum` when available
- [ ] 4.5 Add `OLD_MULTIALLELIC` validation
- [ ] 4.6 Add `flag_extreme_strand_bias` and `category_conflict_resolution` to parquet schema and `_CROSS_SAMPLE_COLS`

## 5. Unify modality_evidence_caller

- [x] 5.1 Replace `compute_modality_evidence` with tier-based lookup: C1→cross_modality, C2→dna_confident, C3+C4→rna_rescued, C5-C7→low_confidence
- [x] 5.2 Keep `modality_evidence_dp` as depth-based (DNA_DP>30 AND RNA_DP>30)
- [ ] 5.3 Update `_WISE_METRICS` and rescue statistics functions to use the new modality_evidence_caller definition (deferred — metrics auto-adapt via column name)
- [x] 5.4 Verify that 1-DNA+1-RNA variants (C4) are now classified as `rna_rescued` not `cross_modality` (verified via Python test)

## 6. Fix --exclude-multiallelic-conflict

- [x] 6.1 Change `build_unified_filter` to check `flag_category_conflict` instead of `flag_vaf_overflow`
- [x] 6.2 Update help text to match the new behavior

## 7. Fix Multi-allelic Visualizations

- [ ] 7.1 Fix `plot_allele_balance_scatter` to group by (CHROM, POS) for per-allele VAF balance (deferred to fix-broken-visualizations)
- [ ] 7.2 Fix `plot_category_conflict_summary` to use per-sample flag data (deferred to fix-broken-visualizations)
- [ ] 7.3 Register `multiallelic_class` in `_COLOR_REGISTRY` (deferred to fix-broken-visualizations)

## 8. Remove Cyvcf2 from Seq2neo Stats

- [x] 8.1 Remove `from cyvcf2 import VCF` from `caller_parser.py`
- [x] 8.2 Remove cyvcf2 fallback path in `_parse_one_caller` — replaced with hard RuntimeError
- [x] 8.3 Remove top-level cyvcf2 import from `caller_parser.py` (pre-norm keeps local import as optional)
- [x] 8.4 Remove `from cyvcf2 import VCF` from `rescue_parser.py` — function now delegates to Rust parser
- [x] 8.5 Remove `--vcf-parser` CLI flag from `cli.py`
- [x] 8.6 Rust startup check exists in bam_stats.py; caller_parser raises ImportError if Rust unavailable

## 9. Tests

- [x] 9.1 Tests for corrected classification: verified via Python test (true_multi_allelic with 3 alleles)
- [x] 9.2 Test for allele_balance_ratio: two-allele (0.75), three-allele (0.444) — both verified
- [x] 9.3 Tests for corrected flags: heterogeneity (true_multi_allelic only), category_conflict (biological guard), germline_low_vaf (no DP constraint) — all verified
- [x] 9.4 Test for modality_evidence_caller tier mapping: C1→cross_modality, C4→rna_rescued — verified
- [ ] 9.5 Update existing tests that use nonexistent file paths for pre-norm parsing (deferred)
- [ ] 9.6 Run `nf-test test tests/default.nf.test --profile test,docker` (requires full pipeline)
