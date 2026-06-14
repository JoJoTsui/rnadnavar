## ADDED Requirements

### Requirement: Per-sample level statistics
The system SHALL compute per-sample aggregate statistics including: total variants, pass/filter distribution, VC classification counts, variant type distribution, Ti/Tv ratio, mean VAF (DNA/RNA), mean DP (DNA/RNA), mean REF_DP/ALT_DP (DNA/RNA), caller support distribution, cross-modality counts, COSMIC/gnomAD coverage, REDIportal levels, and per-tier distribution.

#### Scenario: Sample summary output
- **WHEN** a sample is processed
- **THEN** a row in `sample_summary.csv` contains all per-sample metrics

### Requirement: Per-tier level statistics
The system SHALL compute statistics aggregated across all samples at the tier level: for each CxDy tier, compute total variant count, mean VAF/DP/REF_DP/ALT_DP, variant type distribution, GT concordance, and caller overlap distribution.

#### Scenario: Tier-level output
- **WHEN** all samples are processed
- **THEN** `tier_summary.csv` contains one row per CxDy tier with tier-level aggregates across all samples

### Requirement: Whole-dataset level statistics
The system SHALL compute overall statistics across the entire dataset: total variant count, overall variant type distribution, overall Ti/Tv ratio, overall VC classification distribution, overall caller support histogram, and overall tier distribution.

#### Scenario: Dataset-level output
- **WHEN** all samples are processed
- **THEN** `dataset_summary.csv` contains one row with whole-dataset aggregate metrics

### Requirement: Per-sample variant distribution plot
The system SHALL generate a per-sample variant count bar chart with sample ID on the x-axis, ordered by variant count descending, colored by disease or set.

#### Scenario: Per-sample variant counts
- **WHEN** variant statistics are computed for all samples
- **THEN** a horizontal bar chart is generated showing each sample's total variant count, with sample ID labels on the y-axis
