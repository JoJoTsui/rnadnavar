## MODIFIED Requirements

### Requirement: Variants are tiered using CxDy hybrid system
The system SHALL compute CxDy tiers for each variant using the existing `TieringEngine` from `bin/vcf_stats/tiering_engine.py`. Tiers SHALL be based on category-concordant DNA/RNA caller counts (C1-C7) and database evidence (D0-D1), derived from FILTERS_NORMALIZED / FILTER_NORMALIZED_* INFO fields in the rescue VCF. The system SHALL export `FINAL_TIER_ORDER` (the full ordered list of 14 tier identifiers from C1D1 through C7D0) for consumption by downstream modules including the visualizer.

#### Scenario: Variant with ≥2 DNA + ≥2 RNA callers and database support
- **WHEN** a variant has 2 DNA callers and 2 RNA callers concordant with its FILTER category, AND has gnomAD or COSMIC annotation
- **THEN** the variant is assigned tier C1D1 with quality score 160

#### Scenario: Variant with no caller support and no database
- **WHEN** a variant has 0 DNA callers and 0 RNA callers concordant with its FILTER category, AND no database annotations
- **THEN** the variant is assigned tier C7D0 with quality score 20

#### Scenario: FINAL_TIER_ORDER exported
- **WHEN** a downstream module imports from tiering_stats
- **THEN** FINAL_TIER_ORDER is available as ["C1D1","C1D0","C2D1","C2D0","C3D1","C3D0","C4D1","C4D0","C5D1","C5D0","C6D1","C6D0","C7D1","C7D0"]

## ADDED Requirements

### Requirement: modality_evidence_caller maps from final_tier
The system SHALL compute `modality_evidence_caller` by direct lookup from `final_tier`: C1D0/C1D1 → `cross_modality`, C2D0/C2D1 → `dna_confident`, C3D0/C3D1/C4D0/C4D1 → `rna_rescued`, C5-C7 → `low_confidence`. This replaces the previous raw-caller-count condition (`N_DNA ≥ 1 AND N_RNA ≥ 1` for cross_modality) which over-classified weak single-caller-per-modality variants as cross-modality.

#### Scenario: C1 tier maps to cross_modality
- **WHEN** final_tier=C1D1
- **THEN** modality_evidence_caller="cross_modality"

#### Scenario: C4 tier maps to rna_rescued (not cross_modality)
- **WHEN** final_tier=C4D0 (1 DNA + 1 RNA concordant)
- **THEN** modality_evidence_caller="rna_rescued"
