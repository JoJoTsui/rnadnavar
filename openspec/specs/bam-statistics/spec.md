# bam-statistics

## Purpose

Per-sample per-modality BAM statistics and per-variant BAM-level metrics for the seq2neo variant statistics pipeline.

## Requirements

### Requirement: Per-sample per-modality BAM statistics
The system SHALL compute per-sample BAM statistics for DNA and RNA modalities including: total reads, mapped reads, mapping rate, mean coverage, and mean insert size. Statistics SHALL be computed from alignment BAM files located in the sample's output directory.

#### Scenario: DNA BAM statistics for a sample
- **WHEN** a sample has a DNA alignment BAM file
- **THEN** a `bam_stats.csv` file is generated with one row per sample per modality containing read counts, mapping rate, and coverage metrics

#### Scenario: Sample without BAM file
- **WHEN** a sample's BAM file does not exist
- **THEN** that modality's BAM stats are null-filled and a warning is logged

### Requirement: Per-variant BAM-level statistics in parquet
The system SHALL append per-variant BAM-level metrics to `variant_details.parquet` from available caller FORMAT fields including: strand bias (SB from Mutect2), fragment allele depth (FAD from Mutect2), and per-allele counts for Strelka (AU/CU/GU/TU tier1).

#### Scenario: Mutect2 strand bias extraction
- **WHEN** a variant has a Mutect2 caller with SB in FORMAT
- **THEN** SB values are stored as `DNA_mutect2_SB` and `RNA_mutect2_SB` columns in the variant_details parquet

#### Scenario: Strelka allele count extraction
- **WHEN** a variant has a Strelka caller with AU/CU/GU/TU in FORMAT
- **THEN** AU[0] is stored as `DNA_strelka_AU` (A-allele count tier1) in the variant_details parquet
