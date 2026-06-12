## MODIFIED Requirements

### Requirement: Chart functions accept pre-aggregated data
Chart functions SHALL accept pre-computed aggregation DataFrames instead of raw variant-level DataFrames. Functions that previously received the full combined dataset and performed their own aggregation SHALL now receive the aggregated result directly.

#### Scenario: VC distribution chart with pre-grouped data
- **WHEN** `plot_vc_distribution` receives a DataFrame with columns (set_number, VC, count) that was pre-computed via `group_by(["set_number", "VC"]).agg(pl.len())`
- **THEN** it SHALL produce the same stacked bar chart as when receiving raw variant data

#### Scenario: VAF boxplot with pre-computed quartiles
- **WHEN** `plot_vaf_boxplot_per_tier` receives a DataFrame with pre-computed quartile values per caller per tier
- **THEN** it SHALL produce box-and-whisker plots without loading per-variant VAF data into memory

### Requirement: Remove pandas dependency from visualization
Chart functions SHALL NOT call `.to_pandas()` on polars DataFrames. All altair chart specifications SHALL use polars DataFrames directly.

#### Scenario: No pandas conversion
- **WHEN** any chart function is called
- **THEN** no `.to_pandas()` call SHALL be made; altair SHALL receive polars DataFrames
