## Phase 1: Rust Windowed Queries (DONE)
- [x] 1.1-1.8 Windowed query implementation with chromosome boundary clamping

## Phase 2: 4-Column Join Fix (DONE)
- [x] 2.1-2.4 REF/ALT columns in pileup output, 4-column join

## Phase 3: Wire --pileup-mode Flag (DONE)
- [x] 3.1-3.4 pileup_mode threaded through all call sites

## Phase 4: Tests (4/7 DONE)
- [x] 4.1-4.5, 4.7 Tests for REF/ALT, multiallelic join, filtered mode, 10K perf, join key
- [ ] 4.6 Integration test for pileup columns in parquet

## Phase 5: Verification (3/6 DONE)
- [x] 5.1-5.3 Rust rebuild, pileup tests, unit suite
- [ ] 5.4-5.6 Pipeline verification

## Phase 6: Shared BED Processing (NEW)

- [ ] 6.1 Replace `sum_bed_regions()` with `read_and_merge_bed(path, gap=100_000)`
- [ ] 6.2 Return both `bed_total` (int) and `bed_regions` (list of (chrom, start, end))
- [ ] 6.3 Merge adjacent intervals within `gap` bp to reduce query count
- [ ] 6.4 Update `compute_all_bam_stats()` signature to accept `bed_regions`
- [ ] 6.5 Pass `bed_regions` to Rust `whole_genome_stats()` and pysam fallback

## Phase 7: BAM Stats On-Target Coverage (NEW)

- [ ] 7.1 Add `bed_regions: Option<Vec<(String, i64, i64)>>` parameter to Rust `whole_genome_stats()`
- [ ] 7.2 Build per-chromosome BED interval list for O(reads + regions) lookup
- [ ] 7.3 Track `on_target_bases` separately from `total_query_length`
- [ ] 7.4 Coverage = on_target_bases / bed_total when BED provided
- [ ] 7.5 Update pysam fallback `_compute_bam_stats_pysam()` with same logic
- [ ] 7.6 Remove post-hoc coverage recalculation hack (lines 213-218 in bam_stats.py)

## Phase 8: Rust Combined Multi-BAM Pileup (NEW)

- [ ] 8.1 New function `pileup_variants_multi(bam_paths, positions, bed_regions?)` in pileup.rs
- [ ] 8.2 Group positions into regions once (shared across all BAMs)
- [ ] 8.3 Region source: BED intervals if provided, else 1Mb sliding windows
- [ ] 8.4 Per-region: query all BAMs, process reads with binary search inner loop
- [ ] 8.5 Binary search: sorted position array, per-read `binary_search_by` for overlap range
- [ ] 8.6 Clamp regions to chromosome lengths from BAM header
- [ ] 8.7 Return `Vec<(bam_label, Vec<PileupResult>)>` — all results in one FFI call
- [ ] 8.8 Expose as `#[pyfunction]` in lib.rs
- [ ] 8.9 Python wrapper in `rust_bam.py`: `pileup_variants_multi(bam_paths, positions, bed_regions?)`
- [ ] 8.10 Update `process_single_sample()` to use combined function instead of per-BAM loop
- [ ] 8.11 Pass `bed_regions` from CLI through to pileup when `--bed` is provided

## Phase 9: Parallelize BAM Stats + Variant Processing (NEW)

- [ ] 9.1 Launch BAM stats in background ThreadPoolExecutor before variant processing
- [ ] 9.2 Collect BAM stats result after variant processing completes
- [ ] 9.3 Handle BAM stats failure gracefully (log warning, continue without BAM data)

## Phase 10: Tests for New Features (NEW)

- [ ] 10.1 `test_read_and_merge_bed_total` — merged total matches sum of input intervals
- [ ] 10.2 `test_bed_merge_within_gap` — intervals within 100Kb gap are merged
- [ ] 10.3 `test_bam_stats_on_target_coverage` — BED-filtered coverage ≤ total coverage
- [ ] 10.4 `test_pileup_multi_bam_parity` — multi-BAM function matches individual calls
- [ ] 10.5 `test_pileup_binary_search_parity` — binary search matches HashMap approach
- [ ] 10.6 `test_pileup_bed_guided_region_count` — BED mode creates fewer regions
- [ ] 10.7 `test_pileup_multi_bam_bed_integration` — BED regions passed through correctly

## Phase 11: Verification (NEW)

- [ ] 11.1 Rebuild Rust module
- [ ] 11.2 Run full test suite → 0 failures
- [ ] 11.3 Run pipeline with `--bed` → verify on-target coverage accuracy
- [ ] 11.4 Run pipeline with `--bed --no-pileup` → verify BED-guided pileup regions
- [ ] 11.5 Run pipeline without `--bed` → verify 1Mb window fallback
- [ ] 11.6 Measure samples/hour with 4 workers → verify > 50% improvement
