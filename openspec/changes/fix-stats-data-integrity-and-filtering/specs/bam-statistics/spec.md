## ADDED Requirements

### Requirement: Whole-genome coverage bin computation
The Rust BAM statistics backend SHALL compute coverage bin percentages (`cov_1x_pct`, `cov_10x_pct`, `cov_20x_pct`, `cov_50x_pct`, `cov_100x_pct`) for whole-genome mode in addition to the existing BED-guided mode. The computation SHALL use a streaming depth histogram during the existing whole-genome scan, accumulating per-threshold counters for every alignment match base. The denominator SHALL be the sum of reference sequence lengths from the BAM header.

#### Scenario: WG coverage bins from streaming scan
- **WHEN** a BAM file is processed in whole-genome mode (no BED regions)
- **THEN** `cov_1x_pct`, `cov_10x_pct`, `cov_20x_pct`, `cov_50x_pct`, and `cov_100x_pct` are populated with non-null float values

#### Scenario: WES coverage bins unchanged
- **WHEN** a BAM file is processed with BED regions
- **THEN** the existing `coverage_bins()` BAI-guided function is used
- **AND** coverage bin values are identical to current behavior

#### Scenario: BAM without index
- **WHEN** a BAM file has no BAI index in WG mode
- **THEN** coverage bins are still computed (streaming scan does not require index)
- **AND** mean_coverage, total_reads, and other core metrics are unchanged

### Requirement: Python-side WG coverage bin passthrough
The Python `compute_bam_stats()` function SHALL call the Rust coverage bin computation for both WG and WES modes. In WG mode, the coverage bin result SHALL be extracted from the same `_compute_bam_stats_rust()` call (no separate FFI call). The `coverage_bins()` Rust function SHALL continue to be called only for WES mode where BAI-indexed random access is appropriate.

#### Scenario: WG mode — coverage from bam_stats result
- **WHEN** `compute_bam_stats(bam_path)` is called without bed_regions
- **THEN** cov_*_pct values are extracted from the Rust result dict
- **AND** no separate `coverage_bins()` call is made

#### Scenario: WES mode — coverage from coverage_bins
- **WHEN** `compute_bam_stats(bam_path, bed_regions=regions)` is called with BED regions
- **THEN** `coverage_bins()` is called as before
- **AND** cov_*_pct values are populated from the coverage_bins result
