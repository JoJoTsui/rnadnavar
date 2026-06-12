## MODIFIED Requirements

### Requirement: Per-sample per-modality BAM statistics
The system SHALL compute per-sample BAM statistics for THREE modalities (DNA normal DN, DNA tumor DT, RNA tumor RT) instead of the current two (DNA, RNA). Statistics SHALL be computed from BAM files located in `preprocessing/mapped/{prefix}{suffix}/`.

#### Scenario: Three BAM types per sample
- **WHEN** a sample has DN, DT, and RT BAM files
- **THEN** `bam_stats.csv` contains three rows for that sample, one per BAM type

### Requirement: Per-variant BAM-level statistics in parquet
The system SHALL append per-variant BAM pileup columns to `variant_details.parquet` including: BAM_DP_{type}, BAM_REF_{type}, BAM_ALT_{type}, BAM_F1R2_ref_{type}, BAM_F2R1_ref_{type}, BAM_F1R2_alt_{type}, BAM_F2R1_alt_{type}, BAM_MEAN_BQ_{type}, BAM_MEAN_MQ_{type} for each of DN, DT, RT.

#### Scenario: Per-variant BAM columns in parquet
- **WHEN** BAM pileup is complete
- **THEN** variant_details.parquet contains 27 new BAM pileup columns (9 metrics × 3 BAM types)

## ADDED Requirements

### Requirement: BAM file auto-discovery with index validation
The system SHALL automatically discover BAM files for DN, DT, RT types and validate that .bai index files exist and are current. Stale indices SHALL trigger automatic reindexing.

#### Scenario: Stale BAM index
- **WHEN** a .bai index is older than the .bam file
- **THEN** a warning is logged and samtools index is run to rebuild the index
