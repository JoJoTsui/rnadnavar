## Phase 1: Rust Windowed Queries

- [ ] 1.1 Implement `pileup_variants_windowed()` — group positions into ~1Mb genomic windows
- [ ] 1.2 Build HashMap<(chrom, window_start), Vec<(pos, orig_idx, ref_byte, alt_byte)>> for window grouping
- [ ] 1.3 Within each window, build pos_map HashMap<i64, Vec<(orig_idx, ref_byte, alt_byte)>>
- [ ] 1.4 One `reader.query()` per window, iterate reads, match to positions via pos_map
- [ ] 1.5 Handle reads overlapping multiple positions (150bp read → up to 150 position checks)
- [ ] 1.6 Return results in original input order via orig_idx
- [ ] 1.7 Fall back to per-position query if BAI index missing or windowed query fails
- [ ] 1.8 Expose as `#[pyfunction]` in lib.rs

## Phase 2: 4-Column Join Fix

- [ ] 2.1 Add REF and ALT columns to pileup result in `rust_bam.py` (line 46-47)
- [ ] 2.2 Update pysam fallback in `_pileup_pysam()` to include REF and ALT columns
- [ ] 2.3 Change join key in `cli.py` from `["CHROM", "POS"]` to `["CHROM", "POS", "REF", "ALT"]`
- [ ] 2.4 Update column rename logic to exclude REF and ALT from prefix renaming

## Phase 3: Wire --pileup-mode Flag

- [ ] 3.1 Pass `pileup_mode` through `process_single_sample()` and `_process_worker()` signatures
- [ ] 3.2 When `pileup_mode == "filtered"`, exclude `FILTER == "NoConsensus"` positions before pileup
- [ ] 3.3 Add progress log: print position count and BAM type before each pileup call
- [ ] 3.4 Thread `pileup_mode` from CLI args through all call sites (spawn, thread, sequential)

## Phase 4: Tests

- [ ] 4.1 `test_pileup_output_includes_ref_alt` — REF and ALT columns present in result
- [ ] 4.2 `test_pileup_multiallelic_join_correct` — 4-column join matches correct alleles
- [ ] 4.3 `test_pileup_windowed_matches_per_position` — windowed query parity with per-position for 100 positions
- [ ] 4.4 `test_pileup_10k_positions_completes_quickly` — 10K positions < 60 seconds
- [ ] 4.5 `test_pileup_mode_filtered_excludes_noconsensus` — filtered positions < all positions
- [ ] 4.6 `test_process_single_sample_writes_pileup_columns` — BAM_DN_DP etc. in parquet
- [ ] 4.7 `test_pileup_join_uses_4_columns` — join key is (CHROM, POS, REF, ALT)

## Phase 5: Verification

- [ ] 5.1 Rebuild Rust module: `maturin develop --release`
- [ ] 5.2 Run existing pileup tests → 0 regressions
- [ ] 5.3 Run full unit test suite → 0 failures
- [ ] 5.4 Run 2-sample pipeline with pileup enabled → completes in < 30 minutes
- [ ] 5.5 Verify pileup columns in parquet, join correctness at multiallelic sites
- [ ] 5.6 Verify `--pileup-mode filtered` reduces position count
