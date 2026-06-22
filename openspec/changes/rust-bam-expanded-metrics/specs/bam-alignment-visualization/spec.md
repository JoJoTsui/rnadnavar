# bam-alignment-visualization (delta)

## MODIFIED Requirements

### Requirement: Sample-wise BAM metrics chart
The `plot_bam_metrics_sample_wise()` chart SHALL render all available BAM metrics, including `duplication_rate_pct`, `properly_paired_pct`, and `insert_size_stddev` when they are populated in `bam_stats_df`. When these columns contain only null values, the corresponding chart panels SHALL be empty with no bars rendered.

#### Scenario: Populated metrics render bars
- **WHEN** `bam_stats_df` has non-null `duplication_rate_pct` values
- **THEN** the duplication rate panel shows bars per sample colored by bam_type

#### Scenario: Null metrics render empty panels
- **WHEN** `bam_stats_df` has null `insert_size_stddev` for all rows
- **THEN** the insert size stddev panel renders with axes but no bars

### Requirement: Coverage distribution chart
The `plot_bam_coverage_distribution()` chart SHALL render coverage bin percentages when `cov_1x_pct` through `cov_100x_pct` columns are populated. When all coverage columns are null, the chart SHALL be silently skipped.

#### Scenario: Coverage bins populated
- **WHEN** `bam_stats_df` has non-null coverage bin values (BED mode with Rust coverage_bins)
- **THEN** a grouped bar chart of mean coverage percentages per bam_type and set_number is rendered

#### Scenario: Coverage bins null
- **WHEN** all `cov_*_pct` columns contain only null values
- **THEN** `plot_bam_coverage_distribution()` returns early with no output file
