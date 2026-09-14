# variant-tiering-stats

## Purpose

Variant statistics organized by the CxDy hybrid tiering system (C1-C7 caller support × D0-D1 database evidence) integrated from the existing `TieringEngine` in `bin/vcf_stats/tiering_engine.py`.

## ADDED Requirements

### Requirement: Hard filter flags included in tier statistics
The system SHALL include hard filter flag columns (`flag_hard_*`, `hard_filter_flags`, `n_hard_flags`) in the confidence tier computation and all downstream tier-level statistics, enabling analysis of how hard filter conditions correlate with confidence tiers.

#### Scenario: Confidence tier breakdown includes hard filter context
- **WHEN** `confidence_tier_summary.tsv` is generated
- **THEN** variants with no hard filter flags SHALL be distinguishable from those with active flags
- **AND** the confidence tier computation SHALL consider hard filter flags alongside other soft flags when assigning LOW confidence

#### Scenario: Hard filter flag counts per tier
- **WHEN** tier-level summary statistics are computed
- **THEN** each tier's output SHALL include the distribution of `n_hard_flags` values for variants in that tier

## MODIFIED Requirements

### Requirement: Variants are tiered using CxDy hybrid system
The system SHALL compute CxDy tiers for each variant using the existing `TieringEngine` from `bin/vcf_stats/tiering_engine.py`. Tiers SHALL be based on category-concordant DNA/RNA caller counts (C1-C7) and database evidence (D0-D1), derived from FILTERS_NORMALIZED / FILTER_NORMALIZED_* INFO fields in the rescue VCF. Variants that match hard filter conditions SHALL still receive confidence tier assignments; hard filter flags SHALL be recorded alongside tier assignments rather than causing exclusion.

#### Scenario: Variant with ≥2 DNA + ≥2 RNA callers and database support
- **WHEN** a variant has 2 DNA callers and 2 RNA callers concordant with its FILTER category, AND has gnomAD or COSMIC annotation
- **THEN** the variant is assigned tier C1D1 with quality score 160
- **AND** the variant's hard filter flags are preserved in the output

#### Scenario: Variant with no caller support and no database
- **WHEN** a variant has 0 DNA callers and 0 RNA callers concordant with its FILTER category, AND no database annotations
- **THEN** the variant is assigned tier C7D0 with quality score 20
- **AND** the variant's hard filter flags are preserved in the output

#### Scenario: Variant with hard filter flags still tiered
- **WHEN** a variant matches the `no_caller_support` hard filter condition
- **THEN** the variant SHALL still receive a CxDy tier assignment
- **AND** its `n_hard_flags` value SHALL be at least 1
- **AND** its confidence tier MAY reflect the hard filter condition via soft flag logic

## REMOVED Requirements

### Requirement: Hard filter drops variants before tiering
**Reason**: Hard filter conditions are now recorded as flag columns rather than causing variant exclusion. This separates observation (statistics stage) from filtering (dataset preparation stage). The `--no-hard-filter` CLI argument is also removed since the hard filter no longer performs any exclusion.
**Migration**: Replace `build_hard_filter_expr()` with `build_hard_filter_flag_exprs()`. Remove the `filter(~hard_expr)` call. Remove the `--no-hard-filter` argument from the CLI parser.
