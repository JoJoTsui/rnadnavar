## ADDED Requirements

### Requirement: Caller overlap matrix and support distribution outputs
The system SHALL compute and write `caller_overlap_matrix.tsv` (pairwise caller co-occurrence matrix) and `caller_support_distribution.tsv` (histogram of N_SUPPORT_CALLERS). The existing `caller_overlap.tsv` (which contains a tier distribution, not a caller overlap) SHALL be renamed to `tier_distribution.tsv`.

#### Scenario: Caller overlap matrix generated
- **WHEN** variant statistics are computed for a dataset with 6 callers
- **THEN** caller_overlap_matrix.tsv contains a 6×6 symmetric matrix of pairwise co-occurrence counts

#### Scenario: Support distribution generated
- **WHEN** variant statistics are computed
- **THEN** caller_support_distribution.tsv contains one row per N_SUPPORT_CALLERS value (0-6) with variant counts

### Requirement: Both filtered and unfiltered sample and set summaries
The system SHALL emit both filtered and unfiltered versions of `sample_summary`, `set_summary`, and `dataset_summary` when the unified filter is active. Unfiltered versions SHALL be computed from `combined_df` before filter application. Filtered versions SHALL be computed from the filtered `combined_df`.

#### Scenario: Sample summary both versions
- **WHEN** the pipeline runs with a filter applied
- **THEN** sample_summary.tsv (unfiltered, all variants) and sample_summary_filtered.tsv (filtered) are both generated
