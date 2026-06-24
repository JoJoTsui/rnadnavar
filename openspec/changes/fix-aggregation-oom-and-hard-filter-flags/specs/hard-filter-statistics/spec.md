# hard-filter-statistics

## Purpose

Statistics and visualization of hard filter condition breakdowns — counting how many variants match each hard filter condition, without dropping any variants from the dataset.

## Requirements

### Requirement: Per-condition variant counts
The system SHALL compute the number of variants matching each hard filter condition by evaluating each flag column against the combined lazy DataFrame. Counts SHALL be written to `hard_filter_breakdown.tsv` with columns: `condition`, `severity`, `description`, `n_variants`, `pct_of_filtered`.

#### Scenario: Breakdown TSV is written
- **WHEN** hard filter flags have been computed for all variants passing the unified filter
- **THEN** `hard_filter_breakdown.tsv` SHALL contain one row per hard filter condition
- **AND** each row SHALL include the condition name, severity level, description, absolute count, and percentage of variants passing the unified filter

#### Scenario: No conditions match
- **WHEN** no variants match any hard filter condition
- **THEN** `hard_filter_breakdown.tsv` SHALL still be written with all 8 rows and `n_variants` = 0 for each

### Requirement: Hard filter breakdown bar chart
The system SHALL generate a horizontal bar chart showing the count of variants per hard filter condition, colored by severity level (`high` = red, `medium` = orange, `low` = yellow). The chart SHALL be written to `stats/hard_filter/hard_filter_breakdown.png`.

#### Scenario: Bar chart generated
- **WHEN** `hard_filter_breakdown.tsv` has been written
- **THEN** a horizontal bar chart SHALL be generated with conditions on the y-axis sorted by count descending
- **AND** bars SHALL be colored by severity level
- **AND** the chart SHALL include count labels on each bar

#### Scenario: Chart output location
- **WHEN** the output directory is `<outdir>`
- **THEN** the bar chart SHALL be written to `<outdir>/stats/hard_filter/hard_filter_breakdown.png`

### Requirement: Hard filter statistics are observational only
The hard filter breakdown statistics SHALL NOT exclude any variants from downstream analysis. All variants passing the unified filter SHALL remain in the `combined_df` LazyFrame for confidence tiering and all other statistics, regardless of which hard filter conditions they match.

#### Scenario: Variants preserved after hard filter stats
- **WHEN** hard filter breakdown statistics are computed
- **THEN** the `combined_df` LazyFrame SHALL contain the same number of variants before and after the hard filter stats computation
- **AND** downstream confidence tiering SHALL operate on all variants, not just those with `n_hard_flags == 0`

### Requirement: Hard filter flags persisted for downstream use
The `flag_hard_*`, `hard_filter_flags`, and `n_hard_flags` columns SHALL be present in all downstream outputs (confidence tier summary, tier distribution, etc.) so that dataset preparation stages can filter on them.

#### Scenario: Flag columns flow to downstream outputs
- **WHEN** `disease_summary` or `tier_summary` is computed
- **THEN** the output SHALL include counts broken down by `hard_filter_flags` where applicable
- **AND** the `n_hard_flags` column SHALL be available for filtering in the dataset preparation stage
