## MODIFIED Requirements

### Requirement: Count-based charts use efficient queries
Charts that display variant counts on an axis SHALL use a logarithmic scale (`alt.Scale(type="log")`) with equally-spaced decades (1, 10, 100, 1K, 10K, 100K, 1M, 10M). The `_count_scale()` helper SHALL return `alt.Scale(type="log")` instead of `alt.Scale(type="symlog", constant=1)`. Zero-count rows SHALL be omitted from the data before charting (they are naturally absent from `group_by` aggregations).

#### Scenario: Count axis uses log scale
- **WHEN** any chart with a variant count axis is rendered
- **THEN** the Y-axis ticks SHALL be at logarithmic intervals (1, 10, 100, 1K, 10K, 100K, 1M)
- **AND** the spacing between consecutive decades SHALL be visually equal

#### Scenario: Count axis with data spanning 6 orders of magnitude
- **WHEN** a VC distribution chart has counts from 10K (Somatic) to 90M (NoConsensus)
- **THEN** both the 10K bar and the 90M bar SHALL be visually distinguishable
- **AND** the 10K bar SHALL NOT be compressed to less than 5% of chart height

### Requirement: Per-sample variant distribution plot
The system SHALL generate a per-sample variant count bar chart showing each sample's total variant count. The chart SHALL be faceted by `set_number` with `columns=2` to create a 2×2 grid (4 sets). Each facet SHALL show only the samples belonging to that set, with sample_id on the x-axis (rotated labels) and total_variants on the y-axis (log scale).

#### Scenario: Per-sample distribution with set faceting
- **WHEN** 61 samples across 4 sets are processed
- **THEN** the chart shows a 2×2 grid with one panel per set
- **AND** each panel shows only its set's samples on the x-axis

### Requirement: Per-sample tier distribution uses correct axis orientation
The per-sample tier distribution chart (chart 22) SHALL display sample_id on the x-axis (with rotated labels at -45°) and variant counts on the y-axis using log scale. The chart SHALL be faceted by set_number.

#### Scenario: Chart 22 axis orientation
- **WHEN** the per-sample tier distribution is rendered
- **THEN** sample IDs appear on the x-axis (vertical bars)
- **AND** the y-axis uses log scale for variant counts
- **AND** the chart is faceted by set_number

### Requirement: BAM coverage chart uses correct data source
The BAM coverage boxplot chart (`plot_bam_coverage_violin`) SHALL NOT display a color legend when the color encoding encodes the same field as the x-axis. When `metric:N` is used for both x-axis and color, the legend SHALL be suppressed via `legend=None`.

#### Scenario: BAM chart without redundant legend
- **WHEN** the BAM coverage boxplot is rendered with metric on x-axis and color
- **THEN** no color legend is displayed
- **AND** the x-axis labels provide the only categorical identification

### Requirement: Threshold sweep charts use grid layout
The VAF threshold sweep (chart 31) and DP threshold sweep (chart 35) SHALL display classification panels in a grid layout with `columns=3` instead of a vertical stack. The layout SHALL use `.facet(facet="classification:N", columns=3)` producing a 3×N grid.

#### Scenario: VAF threshold sweep grid
- **WHEN** the VAF threshold sweep is rendered with 7 classifications
- **THEN** the chart displays panels in a 3-column grid (3×3 with 2 empty cells)
- **AND** the total chart height SHALL NOT exceed 1200px

### Requirement: Heatmap cells accommodate text labels
Heatmap charts SHALL size cells large enough to display their text labels without overflow. For `plot_fp_cross_tab_heatmap`, cells SHALL be at least 35×35px. Large count values (>999,999) SHALL use SI-abbreviated format (`.2s`, e.g., "94M" instead of "94,003,313").

#### Scenario: FP cross-tab heatmap cell sizing
- **WHEN** the FP cross-tab heatmap is rendered with counts up to 94M
- **THEN** cells are at least 35px wide
- **AND** text labels use abbreviated format ("94M")
- **AND** no text overflows the cell boundary

### Requirement: Filter VAF/DP heatmap aggregates across partitions
The filter VAF/DP cross-tabulation heatmap (chart 43) SHALL aggregate counts by summing across the `partition` column before charting. Each heatmap cell SHALL display exactly one text label showing the total count.

#### Scenario: Chart 43 no text duplication
- **WHEN** the filter VAF/DP heatmap is rendered from data with partition column
- **THEN** each cell displays exactly one text label
- **AND** the count is the sum across all partitions

## ADDED Requirements

### Requirement: Sample-wise charts include set_number metadata
All charts in the `sample/` wise dimension SHALL have access to `set_number` metadata joined from the sample manifest. Sample-wise charts SHALL be faceted by `set_number` with `columns=2`.

#### Scenario: Sample-wise VC distribution faceted by set
- **WHEN** the sample-wise VC distribution chart is rendered
- **THEN** the chart shows a 2×2 grid (4 sets)
- **AND** each panel shows only samples belonging to that set

### Requirement: Pie chart legend excludes zero-count categories
The somatic modality pie chart SHALL filter the color domain to include only categories present in the data. Categories with zero variants SHALL NOT appear in the legend.

#### Scenario: Pie chart with missing Unknown category
- **WHEN** the data has no "Unknown" somatic modality variants
- **THEN** the legend does NOT show an "Unknown" entry
- **AND** no arc slice is drawn for "Unknown"

### Requirement: GT concordance facet handles null group values
The GT concordance chart SHALL handle null values in the facet group column by replacing them with "Unknown" before faceting. No facet panel SHALL be titled "undefined" or "null".

#### Scenario: GT concordance with null set_number
- **WHEN** some variants have null set_number
- **THEN** the facet panel is titled "Unknown" (not "undefined")
