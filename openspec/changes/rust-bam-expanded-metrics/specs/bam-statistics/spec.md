# bam-statistics (delta)

## MODIFIED Requirements

### Requirement: Per-sample BAM statistics metrics
The system SHALL compute all 14 BAM metrics using the Rust `stats_core` backend when available. The 4 expanded metrics — `duplication_rate_pct`, `properly_paired_pct`, `insert_size_stddev`, and `cov_1x_pct` through `cov_100x_pct` — SHALL be computed in Rust and SHALL NOT use pysam post-processing in the Rust path. The pysam fallback path SHALL compute these metrics with a `max_reads` cap of 10,000,000 reads per scan.

#### Scenario: Rust path computes all 14 metrics
- **WHEN** `HAS_RUST_BAM` is true and a BED file is provided
- **THEN** `compute_bam_stats()` returns a dict with all 14 metric keys populated
- **AND** no pysam `for read in bam.fetch()` calls are made in the Rust path

#### Scenario: Pysam fallback with max_reads cap
- **WHEN** `HAS_RUST_BAM` is false
- **THEN** `_compute_duplication_rate()`, `_compute_properly_paired_pct()`, and `_compute_insert_size_stddev()` each stop after 10,000,000 reads
- **AND** the pipeline completes within 30 seconds per BAM file

#### Scenario: WGS mode (no BED) coverage bins are null
- **WHEN** no BED file is provided
- **THEN** `cov_1x_pct` through `cov_100x_pct` are null-filled
- **AND** the other 9 metrics (original 6 + 3 new Rust counters) are populated

## ADDED Requirements

### Requirement: BAM stats future has timeout
The `bam_future.result()` call in `cli.py` SHALL include a `timeout=3600` (1 hour) argument. If the timeout expires, a warning SHALL be printed and the pipeline SHALL continue with an empty BAM stats DataFrame.

#### Scenario: Timeout safety net
- **WHEN** BAM stats computation exceeds 1 hour
- **THEN** a warning is printed: "WARNING: BAM stats timed out after 3600s — BAM charts will be skipped"
- **AND** the pipeline continues with remaining statistics and visualizations

### Requirement: ensure_bam_stats_columns in resume path
The `--resume` path that reloads `bam_stats.tsv` SHALL call `ensure_bam_stats_columns()` on the loaded DataFrame to null-fill any columns missing from older pipeline versions.

#### Scenario: Resume with old TSV
- **WHEN** `--resume` loads a `bam_stats.tsv` lacking the 4 expanded metric columns
- **THEN** the missing columns are null-filled via `ensure_bam_stats_columns()`
- **AND** BAM charts gracefully skip panels with null data
