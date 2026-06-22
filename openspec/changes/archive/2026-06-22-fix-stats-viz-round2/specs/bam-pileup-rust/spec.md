## MODIFIED Requirements

### Requirement: Pileup output matches pysam reference
Rust pileup output SHALL match pysam pileup output for the same positions within a tolerance of 1 count for depth metrics and 0.1 for quality means.

#### Scenario: Output parity with pysam
- **WHEN** 100 random variant positions are pileup'd by both Rust and pysam
- **THEN** DP, REF_DP, and ALT_DP match exactly, and mean_BQ/mean_MQ differ by no more than 0.1

#### Scenario: Pysam not available
- **WHEN** stats_core is built and available
- **THEN** Rust pileup SHALL be used exclusively without pysam fallback
- **AND** if the Rust pileup fails, a RuntimeError SHALL be raised instead of silently falling back

## REMOVED Requirements

### Requirement: Per-position BAM pileup with Rust noodles-bam
**Reason**: Pysam fallback removed. Rust is the only path. The "Missing BAM index" scenario that described pysam fallback is no longer applicable — missing BAI raises an error.
**Migration**: Ensure BAM indices are generated before running the pipeline. Use `samtools index` if indices are missing.

## ADDED Requirements

### Requirement: BAM validation column naming matches pileup output
The BAM validation module SHALL reference pileup columns using the naming convention `BAM_{bam_type}_{metric}` (e.g., `BAM_DT_DP`, `BAM_RT_ALT_DP`, `BAM_DN_REF_DP`). Column prefix checks SHALL use this format, not the inverted `BAM_{metric}_{bam_type}` format.

#### Scenario: BAM validation detects pileup data
- **WHEN** per-sample parquet contains `BAM_DT_DP` and `BAM_RT_DP` columns from pileup
- **THEN** `has_bam_data` SHALL be `True` in bam_validation.tsv
- **AND** per-caller DP/VAF validation metrics SHALL be computed
