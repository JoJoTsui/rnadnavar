## MODIFIED Requirements

### Requirement: Variants are tiered using CxDy hybrid system
The system SHALL compute CxDy tiers for each variant using a Rust implementation that replicates the existing `TieringEngine` logic from `bin/vcf_stats/tiering_engine.py`. Tiers SHALL be based on category-concordant DNA/RNA caller counts (C1-C7) and database evidence (D0-D1), derived from FILTERS_NORMALIZED INFO fields in the rescue VCF. The Rust implementation SHALL produce identical output to the Python `TieringEngine` for all inputs. The computation SHALL release the GIL via `py.detach()`.

#### Scenario: Variant with ≥2 DNA + ≥2 RNA callers and database support
- **WHEN** a variant has 2 DNA callers and 2 RNA callers concordant with its FILTER category, AND has gnomAD or COSMIC annotation
- **THEN** the variant is assigned tier C1D1 with quality score 160

#### Scenario: Variant with no caller support and no database
- **WHEN** a variant has 0 DNA callers and 0 RNA callers concordant with its FILTER category, AND no database annotations
- **THEN** the variant is assigned tier C7D0 with quality score 20

#### Scenario: Rust tiering matches Python TieringEngine
- **WHEN** 1000 random variants are tiered by both the Rust function and the Python TieringEngine
- **THEN** final_tier, caller_tier, database_tier, dna_caller_count, rna_caller_count, and tier_quality are identical for every variant
