# ref-alt-dp-stats

## Purpose

REF_DP and ALT_DP (reference and alternate allele depth) statistics for the seq2neo variant statistics pipeline.

## Requirements

### Requirement: REF_DP and ALT_DP per-variant columns
The system SHALL compute per-variant REF_DP and ALT_DP columns at the per-caller and aggregate level. REF_DP represents the depth of reference-supporting reads and ALT_DP represents depth of alternate-supporting reads.

#### Scenario: Per-caller REF_DP and ALT_DP
- **WHEN** caller VCF data contains AD_REF and AD_ALT for Mutect2/DeepSomatic, and TOR/TAR for Strelka
- **THEN** columns `{caller}_REF_DP` (from AD_REF or TOR[0]) and `{caller}_ALT_DP` (from AD_ALT or TAR[0]) are present in variant_details.parquet

#### Scenario: Aggregated REF_DP and ALT_DP means
- **WHEN** per-caller REF_DP and ALT_DP columns exist for at least one caller
- **THEN** `DNA_REF_DP_mean`, `DNA_ALT_DP_mean`, `RNA_REF_DP_mean`, `RNA_ALT_DP_mean` columns are computed as row-wise means across DNA and RNA callers respectively

### Requirement: REF_DP and ALT_DP statistics in sample summary
The system SHALL include mean REF_DP and mean ALT_DP (DNA and RNA) in the per-sample summary statistics.

#### Scenario: Sample summary with REF/ALT DP
- **WHEN** a sample's variants have REF_DP and ALT_DP columns
- **THEN** `sample_summary.csv` contains columns `mean_dna_ref_dp`, `mean_dna_alt_dp`, `mean_rna_ref_dp`, `mean_rna_alt_dp`
