## ADDED Requirements

### Requirement: Variants are tiered using CxDy hybrid system
The system SHALL compute CxDy tiers for each variant using the existing `TieringEngine` from `bin/vcf_stats/tiering_engine.py`. Tiers SHALL be based on category-concordant DNA/RNA caller counts (C1-C7) and database evidence (D0-D1), derived from FILTERS_NORMALIZED / FILTER_NORMALIZED_* INFO fields in the rescue VCF.

#### Scenario: Variant with ≥2 DNA + ≥2 RNA callers and database support
- **WHEN** a variant has 2 DNA callers and 2 RNA callers concordant with its FILTER category, AND has gnomAD or COSMIC annotation
- **THEN** the variant is assigned tier C1D1 with quality score 160

#### Scenario: Variant with no caller support and no database
- **WHEN** a variant has 0 DNA callers and 0 RNA callers concordant with its FILTER category, AND no database annotations
- **THEN** the variant is assigned tier C7D0 with quality score 20

### Requirement: Per-tier variant statistics
The system SHALL compute per-tier aggregate statistics including: variant count, mean VAF (DNA and RNA), mean DP (DNA and RNA), mean REF_DP (DNA and RNA), mean ALT_DP (DNA and RNA), variant type distribution (SNV/INS/DEL/MNV), and Ti/Tv ratio for each CxDy tier.

#### Scenario: Tier distribution summary
- **WHEN** variant statistics are computed for a sample
- **THEN** a `tier_summary.csv` file is generated with one row per CxDy tier containing counts and mean metrics

### Requirement: Per-tier VAF violin plots
The system SHALL generate violin plots showing VAF distribution per caller, split by CxDy caller tier (C1-C7). DNA and RNA callers SHALL be plotted in separate facets.

#### Scenario: VAF violin plot per tier
- **WHEN** variant data contains tier assignments and per-caller VAF columns
- **THEN** a multi-panel violin plot is generated with one panel per caller tier, showing the VAF distribution of all 6 callers

### Requirement: Per-tier DP and GT analysis
The system SHALL compute and visualize DP and GT statistics broken down by tier. GT concordance among the 4 GT-bearing callers SHALL be computed per tier.

#### Scenario: GT concordance per tier
- **WHEN** variants have GT fields from Mutect2 and DeepSomatic callers
- **THEN** GT agreement counts (2/3/4 callers agree) are computed separately for each CxDy tier
