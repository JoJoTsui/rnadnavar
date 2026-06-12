## ADDED Requirements

### Requirement: Dataset-wide statistics (Level 1)
The system SHALL compute dataset-wide aggregate statistics across all samples and all variants: total variant count, pass/filter distribution, VC classification distribution (Somatic/Germline/Reference/Artifact/RNAedit), variant type distribution (SNV/INS/DEL/MNV), Ti/Tv ratio, caller support histogram (N_SUPPORT_CALLERS 1-6), tier distribution (C1D1..C7D0), overall mean VAF/DP/REF_DP/ALT_DP (DNA and RNA), cross-modality counts, RESCUED counts, COSMIC/gnomAD coverage, REDIportal evidence distribution, and GT concordance summary.

#### Scenario: Dataset summary output
- **WHEN** all 65 samples are processed
- **THEN** `dataset_summary.csv` contains one row with all Level 1 metrics

### Requirement: Per-tier statistics (Level 2)
The system SHALL compute statistics aggregated across all samples, grouped by CxDy final tier: for each tier, compute variant count, mean VAF/DP/REF_DP/ALT_DP (DNA and RNA), variant type distribution, Ti/Tv ratio, N_SUPPORT_CALLERS distribution, cross-modality counts, RESCUED counts, and GT concordance.

#### Scenario: Tier summary output
- **WHEN** all samples are processed with tier assignments
- **THEN** `tier_summary.csv` contains one row per CxDy tier (up to 14 rows: C1D1..C7D0) with all Level 2 metrics

#### Scenario: Tier with no variants
- **WHEN** a CxDy tier has zero variants across all samples
- **THEN** that tier is absent from tier_summary.csv (not a zero row)

### Requirement: Per-sample statistics (Level 3)
The system SHALL compute per-sample aggregate statistics: total variants, pass/filter distribution, VC classification counts, variant type distribution, Ti/Tv ratio, mean VAF/DP/REF_DP/ALT_DP (DNA and RNA), N_SUPPORT_CALLERS distribution, cross-modality counts, RESCUED count, COSMIC/gnomAD/REDIportal counts, per-tier variant counts, and BAM stats (reads, mapping rate, coverage, insert size, MAPQ for DN/DT/RT).

#### Scenario: Sample summary output
- **WHEN** a sample is processed
- **THEN** `sample_summary.csv` contains one row per sample with all Level 3 metrics including BAM stats

### Requirement: Per-sample per-tier statistics (Level 4)
The system SHALL compute statistics for each sample broken down by CxDy tier: for each (sample_id, final_tier) combination, compute variant count, mean VAF/DP/REF_DP/ALT_DP, variant type distribution, Ti/Tv ratio, and N_SUPPORT_CALLERS distribution.

#### Scenario: Sample-tier summary output
- **WHEN** all samples are processed
- **THEN** `sample_tier_summary.csv` contains one row per (sample_id, final_tier) pair with Level 4 metrics
