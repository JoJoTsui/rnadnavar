## MODIFIED Requirements

### Requirement: Per-sample per-modality BAM statistics
The system SHALL compute per-sample BAM statistics for DNA and RNA modalities including: total reads, mapped reads, mapping rate, mean coverage, mean insert size, mean mapq, duplication rate, properly paired rate, and insert size standard deviation. Statistics SHALL be computed from alignment BAM files using the Rust `stats_core.bam_stats_bed` backend. MAPQ=255 reads SHALL be excluded from mean_mapq computation but included in mapped_reads. Coverage metrics (cov_1x_pct through cov_100x_pct) SHALL be computed when a BED file is provided, using correct 1-based coordinate conversion.

#### Scenario: DNA BAM statistics for a sample
- **WHEN** a sample has a DNA alignment BAM file
- **THEN** a `bam_stats.tsv` file is generated with one row per sample per modality containing all 20 columns populated (when BED is provided)

#### Scenario: RNA BAM with STAR MAPQ=255
- **WHEN** a sample has an RNA alignment BAM file mapped by STAR
- **THEN** mean_mapq is computed only over reads with explicit MAPQ (excluding 255)
- **AND** mapping_rate_pct includes all mapped reads

#### Scenario: Sample without BAM file
- **WHEN** a sample's BAM file does not exist
- **THEN** that modality's BAM stats are null-filled and a warning is logged

#### Scenario: Coverage bins with BED starting at position 0
- **WHEN** coverage_bins is called with a BED interval starting at 0
- **THEN** the interval is processed (not skipped) with correct 1-based conversion
- **AND** cov_*_pct values are computed for that interval
