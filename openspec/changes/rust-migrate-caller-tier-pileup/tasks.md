## 1. Update CALLER_CONFIGS for normalized VCF paths

- [x] 1.1 Update CALLER_CONFIGS in `manifest_loader.py` — change `subdir` and `pattern` to point to normalized VCF paths (`normalized/<caller>/*.dec.norm.vcf.gz` for DNA, `vcf_realignment/normalized/<caller>/*.dec.norm.vcf.gz` for RNA)
- [x] 1.2 Verify normalized VCF files exist for all 66 samples using the manifest

## 2. Caller VCF Rust parser (stats_core/src/caller.rs)

- [x] 2.1 Create `stats_core/src/caller.rs` — implement `parse_caller_vcf()` using noodles-vcf + noodles-bgzf
- [x] 2.2 Parse caller VCF header, find sample index by suffix (DT, RT, TUMOR)
- [x] 2.3 Iterate records, match on (CHROM, POS, REF, ALT) 4-tuple target set, accumulate FORMAT fields per caller type
- [x] 2.4 Implement early termination when all target positions found
- [x] 2.5 Handle missing FORMAT fields gracefully (None values)
- [x] 2.6 Register `parse_caller_vcf` as `#[pyfunction]` in `lib.rs` with `py.detach()` for GIL release
- [x] 2.7 Build, install .so, verify basic parsing works on real normalized VCF

## 3. Wire caller Rust parser into Python

- [x] 3.1 Update `caller_parser.py` — `_parse_one_caller` tries Rust `stats_core.parse_caller_vcf()` first, falls back to cyvcf2
- [x] 3.2 Update `parse_all_callers` to use ThreadPoolExecutor with the new Rust path (already parallel, just switch backend)
- [x] 3.3 Update `build_caller_results_lookup` to use (CHROM, POS, REF, ALT) 4-tuple keys

## 4. join_caller_columns polars optimization

- [x] 4.1 Replace Python for-loops with polars `.join()` on ["CHROM", "POS", "REF", "ALT"], how="left"
- [x] 4.2 Handle missing callers (null-fill all columns for that caller)
- [x] 4.3 Handle Strelka AD_REF/AD_ALT from TOR/TAR correctly in the DataFrame join
- [x] 4.4 Verify output columns match previous implementation exactly

## 5. Tiering Rust implementation (stats_core/src/tier.rs)

- [x] 5.1 Create `stats_core/src/tier.rs` — implement `compute_tiers()` as `#[pyfunction]`
- [x] 5.2 Implement FILTERS_NORMALIZED string parsing (regex: `DNA_strelka:Somatic|RNA_mutect2:Germline|...`)
- [x] 5.3 Implement category-concordant caller counting (matches Python `category_matcher.count_concordant_callers`)
- [x] 5.4 Implement database evidence checking (gnomAD_AF > 0.001, COSMIC_CNT > 0, REDI_EVIDENCE truthy)
- [x] 5.5 Implement C1-C7 caller tier assignment from `tier_config.CALLER_TIER_RULES`
- [x] 5.6 Implement D0-D1 database tier assignment
- [x] 5.7 Compute final_tier, caller_tier, database_tier, tier_quality for each variant
- [x] 5.8 Register in `lib.rs` with `py.detach()` for GIL release
- [x] 5.9 Build, test parity: assert 1000 random variants produce identical tiers to Python TieringEngine

## 6. Wire tiering Rust into Python

- [x] 6.1 Update `tiering_stats.py` — `compute_tiers_for_dataframe` delegates to Rust `stats_core.compute_tiers()` when available
- [x] 6.2 Keep Python TieringEngine fallback for when Rust is unavailable
- [x] 6.3 Verify `tier_summary()` output unchanged (uses same tier columns regardless of backend)

## 7. BAM pileup Rust implementation (stats_core/src/pileup.rs)

- [x] 7.1 Verify noodles-bam `reader.query()` compiles and works on a real BAM with CSI/BAI index
- [x] 7.2 Create `stats_core/src/pileup.rs` — implement `pileup_variants()` as `#[pyfunction]`
- [x] 7.3 Read BAM header and CSI/BAI index
- [x] 7.4 For each (CHROM, POS) in target positions, query the BAM region, iterate reads, count bases
- [x] 7.5 Compute DP, REF_DP, ALT_DP, F1R2_ref, F2R1_ref, F1R2_alt, F2R1_alt, mean_BQ, mean_MQ per position
- [x] 7.6 Handle unmapped reads, duplicates, supplementary reads correctly
- [x] 7.7 Register in `lib.rs` with `py.detach()` for GIL release
- [x] 7.8 Build, test parity: assert 100 random positions produce identical metrics to pysam pileup

## 8. Wire BAM pileup Rust into Python

- [x] 8.1 Update `rust_bam.py` — `pileup_variants` tries Rust `stats_core.pileup_variants()` first, falls back to pysam
- [x] 8.2 Update `HAS_RUST_BAM` flag — check for `pileup_variants` on stats_core (was always False before)
- [x] 8.3 Verify BAM validation step in CLI works end-to-end with Rust pileup

## 9. Tests

- [x] 9.1 Add `TestRustCallerParser` — test normalized VCF parsing, 4-column matching, FORMAT extraction per caller type, early termination, GIL release
- [x] 9.2 Add `TestRustTiering` — test tier parity with Python TieringEngine for 100 variants, boundary cases (C1D1, C7D0), FILTERS_NORMALIZED edge cases
- [x] 9.3 Add `TestRustPileup` — test pileup output matches pysam, empty regions, missing BAM/index handling, GIL release
- [x] 9.4 Add `TestJoinOptimization` — test 4-column join correctness, multiallelic matching, missing caller handling
- [x] 9.5 Extend `TestGILRelease` with tests for new Rust functions (caller, tier, pileup)
- [x] 9.6 Run full test suite, ensure 0 regressions
- [x] 9.7 Run e2e tests on 4 samples, verify output parity with pre-Rust baseline
