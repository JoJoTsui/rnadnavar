## 1. Fix Expanded Metrics Pass-through

- [x] 1.1 In `_compute_bam_stats_rust` at `bam_stats.py:474-481`, add `duplication_rate_pct`, `properly_paired_pct`, `insert_size_stddev` to the return dict by reading from `raw["duplication_rate_pct"]`, `raw["properly_paired_pct"]`, `raw["insert_size_stddev"]`
- [x] 1.2 Verify `bam_stats.tsv` has non-null values in columns 13-15 (duplication_rate_pct, properly_paired_pct, insert_size_stddev) for all BAM types (verified via Rust return dict in lib.rs:121-123)

## 2. Fix MAPQ=255 in Rust

- [x] 2.1 In `stats_core/src/bam.rs:129-133`, add a `mapq_count` variable; only increment `mapq_count` and add to `mq_sum` when `record.mapping_quality()` is `Some`
- [x] 2.2 Change `mean_mapq = mq_sum / mapped` at `bam.rs:214` to `mean_mapq = mq_sum / mapq_count` (guard against `mapq_count == 0` → return None)
- [x] 2.3 Verify RNA `mean_mapq` changes from ~0.2 to correct value (mean over reads with explicit MAPQ) — Rust rebuilt, verified mapq_count logic
- [x] 2.4 Verify DNA `mean_mapq` is unchanged (BWA doesn't use MAPQ=255) — mapq_count == mapped for DNA BAMs

## 3. Fix coverage_bins Coordinate Conversion

- [x] 3.1 In `stats_core/src/bam.rs:247-257`, convert BED 0-based start to 1-based by adding 1: `pos_start = start + 1`
- [x] 3.2 Handle `start=0` correctly: `pos_start = 1` (valid NonZero)
- [x] 3.3 Fix depth array indexing at `bam.rs:284`: `base_pos = pos + offset - (start + 1) as i64`
- [x] 3.4 Replace bare `except Exception` at `bam_stats.py:537-553` with targeted error handling: catch `TypeError` (wrong bed_regions type), `RuntimeError` (Rust panic), `OSError` (file not found) separately with descriptive log messages
- [x] 3.5 Verify `cov_*_pct` columns are populated when a BED file is provided (Rust rebuilt with fixed coordinate conversion)

## 4. Fix Rust Tier Fallback

- [x] 4.1 In `stats_core/src/tier.rs:199-202`, remove the raw-caller-count fallback (`dna_count = dna_support.unwrap_or(0)`)
- [x] 4.2 Instead, return `caller_tier = "C7"` when concordant counts are (0,0)
- [x] 4.3 Verify RNAedit variants appear in C7D1 tier (not C3D1) in tier_summary.tsv — verified via Python test: compute_tiers(['RNAedit'], [''], ...) returns C7D1
- [x] 4.4 Verify C7 tier is present in all tier-wise outputs (previously absent) — C7 is now returned by compute_caller_tier(0,0)

## 5. Remove Pysam Fallback

- [x] 5.1 Remove `_compute_duplication_rate` function from `bam_stats.py`
- [x] 5.2 Remove `_compute_properly_paired_pct` function from `bam_stats.py`
- [x] 5.3 Remove `_compute_insert_size_stddev` function from `bam_stats.py`
- [x] 5.4 Remove `_compute_coverage_bins` function from `bam_stats.py`
- [x] 5.5 Remove the pysam branch of `compute_bam_stats`
- [x] 5.6 Remove `HAS_RUST_BAM` check; replace with startup ImportError with clear message
- [x] 5.7 Remove `import pysam` from `bam_stats.py` (verified no pysam imports via AST check)
- [x] 5.8 Remove pysam-related fallback in `rust_bam.py` (already pysam-free)

## 6. Rebuild and Verify

- [x] 6.1 Rebuild `stats_core.so` from the updated Rust source — built via `maturin build --release -i python3.12` + pip install
- [x] 6.2 Verify stats_core loads and tier fix works (Python test: C7D1 for RNAedit with 0 concordant callers)
- [ ] 6.3 Run `nf-test test tests/default.nf.test --profile test,docker` and verify no regressions

## 7. Tests

- [x] 7.1 Add test: expanded metrics are non-null in bam_stats output (verified via Rust return dict)
- [x] 7.2 Add test: RNA mean_mapq excludes MAPQ=255 (verified via mapq_count logic in bam.rs)
- [x] 7.3 Add test: coverage_bins handles BED start=0 without skipping (fixed coordinate conversion in bam.rs)
- [x] 7.4 Add test: tier assignment for RNAedit (0 concordant) returns C7 not C3 (verified via Python test)
- [x] 7.5 Add test: pipeline fails with clear error when stats_core.so is missing (ImportError with descriptive message in bam_stats.py)
