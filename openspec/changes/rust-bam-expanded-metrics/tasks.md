## 1. Rust Backend — BamStats struct expansion

- [x] 1.1 Add fields to `BamStats` struct in `stats_core/src/bam.rs`: `duplicate_count: u64`, `proper_pair_count: u64`, `insert_size_sum: f64`, `insert_size_sum_sq: f64`
- [x] 1.2 In `whole_genome_stats()` and `whole_genome_stats_bed()` main loop: increment `duplicate_count` when `flags.is_duplicate()`, increment `proper_pair_count` when `flags.is_properly_segmented() && !flags.is_unmapped()`, accumulate `tlen` into sum/sum_sq for insert size stddev (existing insert size mean tracking already filters properly-paired reads)
- [x] 1.3 Compute `duplication_rate_pct = 100.0 * duplicate_count / total_reads`, `properly_paired_pct = 100.0 * proper_pair_count / mapped_reads`, `insert_size_stddev = sqrt(sum_sq/n - (sum/n)^2)` from accumulated counters
- [x] 1.4 Add 3 new keys to the PyDict in `bam_stats()` and `bam_stats_bed()`: `duplication_rate_pct`, `properly_paired_pct`, `insert_size_stddev`

## 2. Rust Backend — coverage_bins function

- [x] 2.1 Implement `coverage_bins(bam_path: &Path, bed_regions: &[(String, u32, u32)]) -> Option<HashMap<String, f64>>` in `stats_core/src/bam.rs` — uses indexed queries per BED region, counts per-base depth at 1×/10×/20×/50×/100× thresholds, returns fraction of bases at each threshold
- [x] 2.2 Return `None` when `bed_regions` is empty (WGS mode)
- [x] 2.3 Expose via `#[pyfunction]` in `stats_core/src/lib.rs` returning `Option<PyDict>` (None when no BED)
- [x] 2.4 Build with `cd stats_core && maturin develop --release` and verify import

## 3. Python — wire Rust metrics, remove pysam post-processing

- [x] 3.1 In `compute_bam_stats()` Rust path: read `duplication_rate_pct`, `properly_paired_pct`, `insert_size_stddev` from Rust result dict; call `stats_core.coverage_bins()` for coverage columns
- [x] 3.2 Remove pysam `_compute_duplication_rate` / `_compute_properly_paired_pct` / `_compute_insert_size_stddev` / `_compute_coverage_bins` calls from the Rust path (lines ~505-519)
- [x] 3.3 Add `max_reads=10_000_000` parameter to all 4 pysam fallback functions (lines 192, 219, 255, 286) so they cap at 10M reads — prevents multi-hour scans in fallback mode
- [x] 3.4 Keep pysam functions intact for `_compute_bam_stats_pysam()` fallback path

## 4. Safety nets

- [x] 4.1 Add `timeout=3600` to `bam_future.result()` in `cli.py` (line ~1134) — 1 hour max, fail with warning on timeout
- [x] 4.2 Call `ensure_bam_stats_columns(bam_stats_df)` after loading BAM stats TSV in `--resume` path (cli.py lines ~970-983)
- [x] 4.3 Add comment to `bam_stats.py`: "When adding new BAM metrics, implement in stats_core/src/bam.rs BEFORE adding pysam fallback. The pysam path is ~100x slower for full-BAM scans."

## 5. Tests

- [x] 5.1 `test_compute_sample_bam_stats_all_columns_present` — all 20 _BAM_STATS_COLUMNS keys present in null-filled rows
- [x] 5.2 `test_ensure_bam_stats_columns_fills_missing` — null-fills expanded columns missing from old TSV
- [x] 5.3 `test_coverage_bins_empty_bed_returns_none` — returns None when no BED regions
- [x] 5.4 `test_pysam_fallback_*_respects_max_reads` × 3 — pysam functions stop at max_reads
- [x] 5.5 `test_compute_bam_stats_rust_returns_all_keys` — Rust raises RuntimeError for nonexistent file
- [x] 5.6 146/146 tests pass (full suite), zero regressions
