## Why

The `expand-stats-visualizations` proposal (P2) added 4 expanded BAM metrics: duplication rate, properly paired percentage, insert size standard deviation, and coverage bins (1×/10×/20×/50×/100× thresholds). These were implemented as pysam post-processing functions called AFTER the fast Rust `stats_core.bam_stats()` backend — each function opens the BAM independently and iterates every read. With 183 BAM files × 300M reads each × 3-4 full scans via slow pysam, the BAM stats phase takes ~23 hours instead of ~30 seconds. The `future.result()` call that collects BAM stats has no timeout, so the entire pipeline hangs indefinitely. Additionally, `ensure_bam_stats_columns()` is not called in the `--resume` path, leaving old TSV files with missing columns.

## What Changes

### 1. Rust Backend — 3 counter metrics in existing single-pass loop
- Add `duplicate_count`, `proper_pair_count`, `insert_size_sum`/`insert_size_sum_sq` to `BamStats` struct in `stats_core/src/bam.rs`
- Increment counters in the existing `for result in reader.records()` loop — zero additional I/O, 3 flag checks per record
- Compute `duplication_rate_pct`, `properly_paired_pct`, `insert_size_stddev` from accumulated counters
- Expose 3 new keys in the PyDict returned by `bam_stats()` and `bam_stats_bed()`

### 2. Rust Backend — coverage bins as a new function
- Add `coverage_bins(bam_path, bed_regions)` to `stats_core/src/bam.rs` using indexed BAM queries per BED region
- Count bases at 1×, 10×, 20×, 50×, 100× depth thresholds
- Expose via PyO3 in `stats_core/src/lib.rs`
- Returns `None` for WGS mode (no BED regions provided)

### 3. Python — wire Rust metrics, remove pysam post-processing
- `compute_bam_stats()`: read 3 new keys from Rust dict + call `coverage_bins()` for the 5 coverage columns
- Remove the pysam `_compute_duplication_rate()` / `_compute_properly_paired_pct()` / `_compute_insert_size_stddev()` / `_compute_coverage_bins()` calls from the Rust path
- Keep pysam functions for the fallback path (when `HAS_RUST_BAM` is false)
- Add `max_reads=10_000_000` parameter to pysam fallback functions so they never scan more than 10M reads

### 4. Safety nets
- `bam_future.result(timeout=3600)` — 1-hour timeout in `cli.py` prevents indefinite hang
- `ensure_bam_stats_columns()` called in `--resume` path so old TSVs get null-filled missing columns
- Prevention rule: BAM metrics must be implemented in Rust first, pysam fallback second

## Capabilities

### Modified Capabilities
- `bam-statistics`: Expanded metrics now computed in Rust during the single-pass BAM scan (duplication rate, properly paired pct, insert size stddev) or via a dedicated Rust function (coverage bins). Output columns unchanged — same names, same semantics, just populated with real data instead of None.
- `bam-alignment-visualization`: `plot_bam_coverage_distribution()` and `plot_bam_metrics_sample_wise()` now receive real data for all metrics instead of null-filled columns.

## Impact

- `bin/vcf_stats/seq2neo/stats_core/src/bam.rs` — ~60 lines: 4 new fields in BamStats, counter increments in existing loops, new `coverage_bins()` function
- `bin/vcf_stats/seq2neo/stats_core/src/lib.rs` — ~30 lines: expose 3 new keys + 1 new function via PyO3
- `bin/vcf_stats/seq2neo/bam_stats.py` — ~30 lines: read new Rust fields, remove pysam post-processing, add max_reads to pysam fallbacks
- `bin/vcf_stats/seq2neo/cli.py` — ~5 lines: timeout on future.result(), ensure_bam_stats_columns() in resume path
- Build: `cd stats_core && maturin develop --release` (~30 seconds)
