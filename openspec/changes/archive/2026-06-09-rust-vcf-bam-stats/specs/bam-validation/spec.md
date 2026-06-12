## ADDED Requirements

### Requirement: BAM DP vs caller DP cross-validation
The system SHALL cross-validate BAM pileup DP against caller VCF DP for all 6 callers. For each variant, the absolute difference between BAM_DP (from pileup) and {caller}_DP (from caller VCF) SHALL be computed. Per-sample mismatch rates SHALL be reported.

#### Scenario: DP validation for a sample
- **WHEN** BAM pileup and caller DP exist for a sample's variants
- **THEN** `bam_validation.csv` contains per-sample mismatch statistics (mean absolute difference, mismatch %, correlation) for each caller

### Requirement: BAM VAF vs caller VAF cross-validation
The system SHALL compare BAM-derived VAF (ALT_depth / DP) against caller-computed VAF for Mutect2 and DeepSomatic callers. Strelka callers SHALL be excluded from VAF validation (no pre-computed VAF in FORMAT).

#### Scenario: VAF validation scatter
- **WHEN** BAM pileup VAF and caller VAF exist
- **THEN** a validation report shows correlation and mean absolute difference per caller

### Requirement: Strand bias validation
The system SHALL compare BAM-measured strand bias (F1R2_alt / F2R1_alt ratio) against caller-reported strand bias (SB field for Mutect2). Significant discrepancies (>2× difference) SHALL be flagged.

#### Scenario: Strand bias comparison
- **WHEN** a variant has both BAM strand counts and caller SB field
- **THEN** strand bias metrics are stored in variant_details.parquet for downstream filtering
