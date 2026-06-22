# variant-tiering-stats (delta)

## MODIFIED Requirements

### Requirement: Per-tier variant statistics
The system SHALL compute per-tier aggregate statistics including: variant count, mean VAF (DNA and RNA), mean DP (DNA and RNA), mean REF_DP (DNA and RNA), mean ALT_DP (DNA and RNA), variant type distribution (SNV/INS/DEL/MNV), Ti/Tv ratio, and modality evidence distribution (cross_modality, dna_confident, rna_rescued, low_confidence counts) for each CxDy tier.

#### Scenario: Tier distribution summary with modality evidence
- **WHEN** variant statistics are computed for a sample
- **THEN** a `tier_summary.csv` file is generated with one row per CxDy tier containing counts, mean metrics, AND breakdown by modality_evidence_caller category

## ADDED Requirements

### Requirement: WISE_METRICS includes modality evidence counts
The `_WISE_METRICS` aggregation list SHALL include per-category counts for the `modality_evidence_caller` column: `n_cross_modality` (count where modality_evidence_caller="cross_modality"), `n_dna_confident` (count where modality_evidence_caller="dna_confident"), `n_rna_rescued` (count where modality_evidence_caller="rna_rescued"), `n_low_confidence` (count where modality_evidence_caller="low_confidence"). The deprecated `n_rescued` and `n_cross_modality` metrics SHALL continue to be computed by mapping from modality_evidence for backward compatibility.

#### Scenario: Wise metrics with modality evidence
- **WHEN** a wise summary is computed (e.g., tier-wise, set-wise)
- **THEN** the output includes columns n_cross_modality, n_dna_confident, n_rna_rescued, n_low_confidence
- **AND** the legacy columns n_rescued and n_cross_modality are still present (deprecated)

### Requirement: Rescue statistics use modality_evidence
The `compute_rescue_breakdown`, `compute_rescue_by_filter`, `compute_rescue_cross_tab`, `compute_rescue_vaf_dp`, `sample_rescue_summary`, `compute_rescue_by_tier`, and `compute_rescue_by_caller_support` functions SHALL use `modality_evidence_caller` as their primary grouping column instead of `RESCUED`. The rescue definition SHALL be `modality_evidence_caller IN ("cross_modality", "rna_rescued")`.

#### Scenario: Rescue breakdown with four evidence categories
- **WHEN** rescue breakdown statistics are computed
- **THEN** the result includes counts for each of the 4 modality_evidence categories
- **AND** the "rescued" count is derived as cross_modality + rna_rescued

## REMOVED Requirements

None. All existing requirements remain with updated behavior.
