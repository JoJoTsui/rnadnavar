## MODIFIED Requirements

### Requirement: Variants are tiered using CxDy hybrid system
The system SHALL compute CxDy tiers for each variant using the existing `TieringEngine` from `bin/vcf_stats/tiering_engine.py`. Tiers SHALL be based on category-concordant DNA/RNA caller counts (C1-C7) and database evidence (D0-D1), derived from FILTERS_NORMALIZED / FILTER_NORMALIZED_* INFO fields in the rescue VCF. Tier columns (caller_tier, database_tier, final_tier, tier_quality, dna_caller_count, rna_caller_count) SHALL be verified to exist in `combined_df`'s schema before tier-wise statistics and charts are generated. If tier columns are absent from the schema, the system SHALL recompute tiers for all variants and issue a warning.

#### Scenario: Variant with ≥2 DNA + ≥2 RNA callers and database support
- **WHEN** a variant has 2 DNA callers and 2 RNA callers concordant with its FILTER category, AND has gnomAD or COSMIC annotation
- **THEN** the variant is assigned tier C1D1 with quality score 160

#### Scenario: Variant with no caller support and no database
- **WHEN** a variant has 0 DNA callers and 0 RNA callers concordant with its FILTER category, AND no database annotations
- **THEN** the variant is assigned tier C7D0 with quality score 20

#### Scenario: Tier columns missing from combined_df schema
- **WHEN** combined_df's schema lacks caller_tier, database_tier, or final_tier columns
- **THEN** tiers are recomputed via compute_tiers_for_dataframe()
- **AND** a warning is logged indicating tier recomputation

## ADDED Requirements

### Requirement: Confidence x tier cross-tabulation
The system SHALL compute a cross-tabulation of confidence_tier by final_tier as part of tier-wise statistics, writing `confidence_tier_cross_tab.tsv` to `stats/tier/`. This SHALL allow assessment of which CxDy tiers contribute to each confidence level.

#### Scenario: Confidence x tier breakdown
- **WHEN** variants have both confidence_tier and final_tier columns
- **THEN** a TSV is generated with rows like (HIGH, C1D1, N), (HIGH, C2D1, M), (MEDIUM, C2D0, P), etc.
