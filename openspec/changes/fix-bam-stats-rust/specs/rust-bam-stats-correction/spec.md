## ADDED Requirements

### Requirement: Expanded metrics are passed through from Rust backend
The Python wrapper `_compute_bam_stats_rust` SHALL read `duplication_rate_pct`, `properly_paired_pct`, and `insert_size_stddev` from the Rust backend's return dictionary and include them in the result. These metrics SHALL NOT be null-filled when the Rust backend successfully computes them.

#### Scenario: Expanded metrics populated
- **WHEN** the Rust backend computes bam_stats for a DNA BAM with 100M reads
- **THEN** bam_stats.tsv has non-null values for duplication_rate_pct, properly_paired_pct, and insert_size_stddev

#### Scenario: Expanded metrics null only on backend failure
- **WHEN** the Rust backend fails to compute expanded metrics for a specific BAM
- **THEN** the expanded metric columns are null for that row
- **AND** a warning is logged with the specific error

### Requirement: MAPQ=255 excluded from mean_mapq computation
The Rust backend SHALL exclude reads with MAPQ=255 (SAM-spec "not available" sentinel) from both the numerator and denominator of `mean_mapq`. Reads with MAPQ=255 SHALL still be counted in `mapped_reads` and `mapping_rate_pct`. The `mean_mapq` metric represents the mean mapping quality over reads that have an explicit MAPQ value.

#### Scenario: RNA BAM with STAR MAPQ=255
- **WHEN** an RNA BAM has 99% of reads with MAPQ=255 (STAR uniquely-mapped) and 1% with MAPQ=3 (multi-mappers)
- **THEN** mean_mapq is computed over only the 1% with explicit MAPQ (≈3.0)
- **AND** mapped_reads includes all reads (both MAPQ=255 and MAPQ=3)
- **AND** mapping_rate_pct is not affected

#### Scenario: DNA BAM with BWA MAPQ 0-60
- **WHEN** a DNA BAM has reads with MAPQ in range 0-60 (BWA, no 255 sentinel)
- **THEN** mean_mapq is computed over all mapped reads (same as before)

### Requirement: Coverage bins use correct 1-based coordinates
The Rust `coverage_bins` function SHALL convert BED 0-based half-open coordinates to noodles 1-based inclusive coordinates by adding 1 to the start position. BED intervals starting at position 0 SHALL be handled correctly (converted to 1-based start=1). The per-base depth array indexing SHALL account for the coordinate conversion.

#### Scenario: BED interval starting at 0
- **WHEN** a BED interval is (chr1, 0, 100)
- **THEN** the Rust backend processes it as 1-based [1, 100]
- **AND** coverage is computed for all 100 bases
- **AND** the interval is NOT skipped

#### Scenario: BED interval mid-chromosome
- **WHEN** a BED interval is (chr1, 1000, 2000)
- **THEN** the Rust backend processes it as 1-based [1001, 2000]
- **AND** per-base depth is correctly aligned to the interval

### Requirement: Coverage bins errors are logged not swallowed
The Python wrapper for `coverage_bins` SHALL replace the bare `except Exception` with targeted error handling. Specific error types (missing BAI, invalid BAM path, Rust panic) SHALL be caught and logged with descriptive messages. The `cov_*_pct` columns SHALL be null-filled only when a specific error is identified, not as a catch-all.

#### Scenario: Missing BAI index
- **WHEN** coverage_bins is called on a BAM without a .bai index
- **THEN** a warning is logged: "coverage_bins failed for {bam_path}: missing BAI index"
- **AND** cov_*_pct columns are null for that row

#### Scenario: Successful coverage computation
- **WHEN** coverage_bins is called with valid BAM, BAI, and BED regions
- **THEN** cov_1x_pct, cov_10x_pct, cov_20x_pct, cov_50x_pct, cov_100x_pct are populated with computed values

### Requirement: Pysam fallback removed
The bam_stats module SHALL use the Rust backend as the only BAM statistics computation path. All pysam-specific computation functions (`_compute_duplication_rate`, `_compute_properly_paired_pct`, `_compute_insert_size_stddev`, `_compute_coverage_bins`, and the pysam branch of `compute_bam_stats`) SHALL be removed. A startup check SHALL verify that `stats_core.bam_stats` is available; if not, the pipeline SHALL fail with a clear error message.

#### Scenario: Rust extension unavailable
- **WHEN** stats_core.so cannot be loaded or does not have bam_stats function
- **THEN** the pipeline exits with error: "Rust BAM stats backend (stats_core.bam_stats) is required but could not be loaded"

#### Scenario: No pysam imports in bam_stats
- **WHEN** bam_stats.py is inspected
- **THEN** no pysam-specific computation functions exist
- **AND** no `import pysam` statement exists in the module
