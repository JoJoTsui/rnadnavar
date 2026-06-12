## MODIFIED Requirements

### Requirement: Per-sample level statistics
The system SHALL compute per-sample aggregate statistics including: total variants, pass/filter distribution, VC classification counts, variant type distribution, Ti/Tv ratio, mean VAF (DNA/RNA), mean DP (DNA/RNA), mean REF_DP/ALT_DP (DNA/RNA), caller support distribution, cross-modality counts, COSMIC/gnomAD coverage, REDIportal levels, and per-tier distribution. Caller FORMAT fields SHALL be matched to rescue variants on (CHROM, POS, REF, ALT) — all 4 coordinate columns — to ensure correct matching at multiallelic sites.

#### Scenario: Sample summary output
- **WHEN** a sample is processed
- **THEN** a row in `sample_summary.csv` contains all per-sample metrics

#### Scenario: Multiallelic site matching
- **WHEN** a caller VCF has two records at the same (CHROM, POS) with different REF/ALT
- **THEN** each record is joined to the rescue variant with the matching REF/ALT, and mismatched alleles receive null-filled caller columns
