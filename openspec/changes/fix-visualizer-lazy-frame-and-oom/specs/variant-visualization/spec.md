## MODIFIED Requirements

### Requirement: Dashboard HTML renders all charts
The system SHALL generate a dashboard.html that renders all generated charts. All chart functions SHALL accept `pl.LazyFrame` input and SHALL NOT call methods exclusive to eager `DataFrame` (`.height`, `.iter_rows()`, subscript access) without first calling `.collect()`. Charts SHALL NOT load more than 50,000 rows into memory without explicit sampling.

#### Scenario: Dashboard with lazy frame input
- **WHEN** 23 charts are generated from a lazy scan of 24 parquet files with 58M+ total variants
- **THEN** dashboard.html contains all charts rendered correctly in a browser
- **AND** no `AttributeError` is raised for LazyFrame-incompatible methods

#### Scenario: Memory safety during chart generation
- **WHEN** any chart function processes a lazy frame with 58M+ rows
- **THEN** peak RSS SHALL NOT exceed 5 GB during chart generation
- **AND** the chart SHALL either sample to ≤50K rows or use efficient count queries (`.select(pl.len()).collect().item()`)

### Requirement: VAF distribution chart uses sampling
The VAF distribution per caller chart SHALL use lazy sampling to limit collected rows to 50,000. The row count check SHALL use `.select(pl.len()).collect().item()` instead of `.height`.

#### Scenario: VAF distribution with sampling
- **WHEN** 58M+ variants have per-caller VAF columns
- **THEN** the chart samples at most 50,000 rows before collecting
- **AND** no full VAF column set is materialized without sampling

### Requirement: GT concordance charts use eager iteration after collect
GT concordance charts (overall and per-tier) SHALL call `.collect()` before `.iter_rows()`. The collected DataFrame SHALL be limited to the needed GT columns only.

#### Scenario: GT concordance with lazy input
- **WHEN** the lazy frame contains 4 GT columns for 58M+ variants
- **THEN** `.collect()` is called before `.iter_rows()`
- **AND** only GT columns are materialized

### Requirement: Count-based charts use efficient queries
Charts that only need row counts (COSMIC/gnomAD annotation coverage, caller agreement matrix) SHALL use `.select(pl.len()).collect().item()` instead of loading data columns. No data SHALL be materialized for count queries.

#### Scenario: COSMIC annotation count
- **WHEN** computing "In COSMIC" count
- **THEN** `df.filter(pl.col("COSMIC_ID").is_not_null()).select(pl.len()).collect().item()` is used
- **AND** no COSMIC_ID column data is loaded into memory

### Requirement: Reorganized chart layout
The visualizer SHALL be organized into clear sections: helper functions (`_count_rows`, `_sample_if_large`), aggregate charts (group_by → small result), sampled charts (use `_sample_if_large`), count-based charts (use `_count_rows`), and per-sample/small-data charts.

#### Scenario: Code organization
- **WHEN** a developer reads visualizer.py
- **THEN** helper functions appear first, followed by chart functions grouped by data access pattern
- **AND** every function's data loading pattern is explicit (no hidden `.pipe(_maybe_collect)`)

### Requirement: BAM coverage chart uses correct data source
The BAM coverage violin chart (`plot_bam_coverage_violin`) SHALL be removed from the `combined_df`-based chart list. Its data source (`BAM_DP_*` columns) exists only in `bam_validation.csv`, not in the per-sample variant parquets. The chart SHALL be moved to use `bam_validation` data if re-added.

#### Scenario: BAM coverage chart removed
- **WHEN** generating charts from combined_df (lazy parquet scan)
- **THEN** no chart attempts to access `BAM_DP_*` columns that don't exist in the parquet files
