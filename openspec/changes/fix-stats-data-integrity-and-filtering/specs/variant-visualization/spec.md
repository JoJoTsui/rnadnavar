## ADDED Requirements

### Requirement: combined_df filtering by current run sample IDs
The system SHALL filter `combined_df` to only include variants belonging to the current run's sample IDs before any statistics, charts, or validation are computed. The filter SHALL be applied immediately after the lazy parquet scan at `cli.py:1103`. The filter SHALL use `pl.col("sample_id").is_in(current_sample_ids)` where `current_sample_ids` is derived from the filtered manifest rows.

#### Scenario: Set-specific run excludes other sets
- **WHEN** `--set 1` is specified and 12 samples are in the manifest
- **THEN** `combined_df` contains only variants from those 12 sample IDs
- **AND** variants from previously-run sets 2, 3, 4 (whose parquet files exist in variant_details/) are excluded

#### Scenario: Full run includes all manifest samples
- **WHEN** no `--set` or `--sample-ids` filter is specified
- **THEN** `combined_df` contains variants from all samples in the manifest
- **AND** orphaned parquet files from samples not in the manifest are excluded

#### Scenario: Missing sample filter logs warning
- **WHEN** some sample IDs in the manifest have no matching parquet file
- **THEN** a warning is logged with the count of missing samples
- **AND** processing continues with available parquet files

### Requirement: Disease value normalization to prevent empty string grouping
The system SHALL normalize empty or whitespace-only disease values to "Unknown" in `process_single_sample` before writing per-sample parquet files. Both `disease` and `disease_normalized` columns SHALL be checked. The normalization SHALL apply `str.strip()` and replace empty results with "Unknown".

#### Scenario: Empty disease normalized to Unknown
- **WHEN** a manifest row has disease="" and disease_normalized=""
- **THEN** the per-sample parquet has disease="Unknown" and disease_normalized="unknown"

#### Scenario: Valid disease preserved
- **WHEN** a manifest row has disease="Colorectal Cancer" and disease_normalized="colorectal cancer"
- **THEN** the per-sample parquet preserves both values unchanged

### Requirement: set_number validation warning
The system SHALL log a warning when all processed samples have `set_number=0`, indicating a possible manifest data quality issue. The warning SHALL be displayed once after all samples are processed, not per-sample.

#### Scenario: All sets are 0
- **WHEN** 12 samples are processed and all have set_number=0
- **THEN** a warning is logged: "All samples have set_number=0 — manifest may be missing partition_set"

#### Scenario: Mixed sets
- **WHEN** samples have set_number values 1, 2, 3, 4
- **THEN** no warning is logged

### Requirement: Chart faceting for filter_distribution, redi_evidence, ref_alt_dp_scatter, and dna_vs_rna_dp
The chart functions `plot_filter_distribution`, `plot_redi_evidence`, `plot_ref_alt_dp_scatter`, and `plot_dna_vs_rna_dp` SHALL call `_apply_faceting(chart, group_col)` to produce per-group faceted panels consistent with other chart functions in the same wise registries. The `_apply_faceting` function SHALL use existing facet rules for set_number (2 columns), disease_normalized (4 columns), FILTER (3 columns), and other group columns (3 columns default).

#### Scenario: filter_distribution faceted by set
- **WHEN** `plot_filter_distribution` is called with group_col="set_number" from the set-wise registry
- **THEN** the chart is faceted by set_number with 2 columns and independent x-axes

#### Scenario: redi_evidence faceted by disease
- **WHEN** `plot_redi_evidence` is called with group_col="disease_normalized" from the disease-wise registry
- **THEN** the chart is faceted by disease_normalized with 4 columns

### Requirement: Disease-wise registry includes plot_caller_overlap
The disease-wise chart registry (`_WISE_CHART_REGISTRY["disease"]`) SHALL include `(plot_caller_overlap, {"group_col": "disease_normalized"})` matching the set-wise registry which already includes this chart.

#### Scenario: Disease-wise caller overlap chart
- **WHEN** disease-wise charts are generated
- **THEN** `02_caller_overlap.html` is written to `plots/disease/` showing tier distribution per disease

### Requirement: caller_tier column existence verification
Before generating tier-wise charts, the system SHALL verify that `caller_tier`, `database_tier`, and `final_tier` columns exist in `combined_df`'s schema. If any tier column is missing, the system SHALL compute tiers for the entire `combined_df` using `compute_tiers_for_dataframe()` and warn the user. If tier computation fails, tier-wise charts that depend on the missing column SHALL be skipped with a warning.

#### Scenario: caller_tier present
- **WHEN** all per-sample parquets contain caller_tier column
- **THEN** tier-wise chart generation proceeds without warnings

#### Scenario: caller_tier missing — recompute
- **WHEN** one or more per-sample parquets lack caller_tier (detected via schema check)
- **THEN** tiers are recomputed for combined_df
- **AND** a warning is logged: "caller_tier missing from parquet schema — recomputing tiers"

### Requirement: BAM metrics_sample_wise fixed grid layout
The `plot_bam_metrics_sample_wise` chart SHALL use a fixed panel width (`width=200`) and a 3-column grid layout (`columns=3`) regardless of sample count. The chart SHALL resolve scales as `x="independent", y="independent"` per panel.

#### Scenario: 12 samples in 3-column grid
- **WHEN** 12 samples have BAM metrics
- **THEN** the chart renders as a 3-column wrapped grid with 200px-wide panels

#### Scenario: Multi-set faceting
- **WHEN** samples span multiple set_number values
- **THEN** the chart facets by row=set_number and column=metric with fixed-width panels

### Requirement: BAM metrics bars top_n parameter from CLI
The `plot_bam_metrics_bars` chart SHALL receive `top_n` from the CLI call site. When not specified, `top_n` SHALL default to 30 (changed from 20). The chart SHALL limit to `top_n` most variant-rich samples before rendering.

#### Scenario: top_n limits sample count
- **WHEN** 100 samples have BAM metrics and top_n=30
- **THEN** the chart shows bars for the 30 samples with the most total reads

## MODIFIED Requirements

### Requirement: BAM coverage chart uses correct data source
The BAM coverage violin chart (`plot_bam_coverage_violin`) SHALL be moved to use `bam_validation` data if re-added. Its data source (`BAM_DP_*` columns) exists only in `bam_validation.csv`, not in the per-sample variant parquets.

The BAM `coverage_distribution` chart (`plot_bam_coverage_distribution`) SHALL render with coverage bin data in both whole-genome and WES modes. In whole-genome mode, `cov_*_pct` columns (1x, 10x, 20x, 50x, 100x) SHALL be populated by the Rust BAM statistics backend via a streaming coverage histogram. In WES mode, the existing BED-guided `coverage_bins()` function SHALL continue to be used. The chart SHALL render grouped bars faceted by set_number when multiple sets are present.

#### Scenario: BAM coverage chart removed
- **WHEN** generating charts from combined_df (lazy parquet scan)
- **THEN** no chart attempts to access `BAM_DP_*` columns that don't exist in the parquet files

#### Scenario: WG coverage distribution with data
- **WHEN** BAM statistics are computed in whole-genome mode (no --bed flag)
- **THEN** `cov_1x_pct`, `cov_10x_pct`, `cov_20x_pct`, `cov_50x_pct`, and `cov_100x_pct` columns contain non-null values
- **AND** `coverage_distribution.html` renders with grouped bars for each BAM type and coverage threshold

#### Scenario: WES coverage distribution with BED regions
- **WHEN** BAM statistics are computed with --bed flag and BED regions
- **THEN** the existing BED-guided coverage_bins() function is used
- **AND** `coverage_distribution.html` renders identically to current behavior
