# bam-statistics (delta)

## MODIFIED Requirements

### Requirement: Per-sample BAM statistics metrics
The system SHALL compute 4 additional metrics per BAM file beyond the current 6: `duplication_rate_pct`, `properly_paired_pct`, `insert_size_stddev`, and `coverage_bins` (cov_1x_pct, cov_10x_pct, cov_20x_pct, cov_50x_pct, cov_100x_pct). Coverage bins SHALL be computed only when a BED file is provided (WES mode). The BAM stats TSV output SHALL include these new columns alongside the existing `total_reads`, `mapped_reads`, `mapping_rate_pct`, `mean_coverage`, `mean_insert_size`, `mean_mapq`.

#### Scenario: Expanded BAM stats TSV
- **WHEN** BAM stats are computed for a sample with a BED file
- **THEN** bam_stats.tsv contains 16 columns (existing 11 + 5 new: duplication_rate_pct, properly_paired_pct, insert_size_stddev, cov_1x_pct, cov_10x_pct, cov_20x_pct, cov_50x_pct, cov_100x_pct — note 11 + 5 = 16, but coverage_bins expand to 5 columns)

#### Scenario: No BED file provided
- **WHEN** BAM stats are computed without a BED file (WGS mode)
- **THEN** coverage bin columns are null-filled

## ADDED Requirements

### Requirement: Duplication rate computation
The system SHALL compute duplication rate as the fraction of reads with the duplicate flag (0x400) set in the BAM FLAG field: `n_duplicate / n_total`. Both the Rust backend (`stats_core.bam_stats`) and the pysam fallback SHALL compute this metric.

#### Scenario: Duplication rate
- **WHEN** a BAM file has 1M reads with 50K duplicates
- **THEN** duplication_rate_pct = 5.0

### Requirement: Properly paired fraction computation
The system SHALL compute properly paired fraction as the fraction of mapped reads where both mates are mapped in a proper pair (FLAG has both 0x1 and 0x2 set, and 0x4, 0x8 unset).

#### Scenario: Properly paired fraction
- **WHEN** a BAM has 900K mapped reads of which 850K are properly paired
- **THEN** properly_paired_pct = 94.4
