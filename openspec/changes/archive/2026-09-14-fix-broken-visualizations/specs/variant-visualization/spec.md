## MODIFIED Requirements

### Requirement: Dashboard HTML renders all charts
The system SHALL generate a dashboard.html that renders all generated charts. Charts SHALL be embedded as body-content fragments extracted from altair's `to_html(output_div=..., inline=False)` output, assembled into a single valid HTML document with shared vega-embed CDN loaded once in the `<head>`. The `_extract_body_content` function SHALL strip `<script>` blocks before searching for the `<body>` tag to avoid matching `<body>` strings inside JavaScript template literals. All chart functions SHALL accept `pl.LazyFrame` input and SHALL NOT call methods exclusive to eager `DataFrame` (`.height`, `.iter_rows()`, subscript access) without first calling `.collect()`. Charts SHALL NOT load more than 50,000 rows into memory without explicit sampling with a fixed seed.

#### Scenario: Dashboard with multiple charts
- **WHEN** 12 charts are generated and passed to generate_dashboard()
- **THEN** dashboard.html contains all 12 charts rendered correctly in a browser
- **AND** no raw JavaScript library fragments appear as visible text

#### Scenario: Dashboard with lazy frame input
- **WHEN** 23 charts are generated from a lazy scan of 24 parquet files with 58M+ total variants
- **THEN** dashboard.html contains all charts rendered correctly in a browser
- **AND** no `AttributeError` is raised for LazyFrame-incompatible methods

#### Scenario: Memory safety during chart generation
- **WHEN** any chart function processes a lazy frame with 58M+ rows
- **THEN** peak RSS SHALL NOT exceed 5 GB during chart generation
- **AND** the chart SHALL either sample to ≤50K rows with a fixed seed or use efficient count queries

## ADDED Requirements

### Requirement: Coverage distribution plot guards against null data
The `plot_bam_coverage_distribution` function SHALL check whether all `cov_*_pct` columns contain only null values before rendering. If all values are null, the function SHALL return early without producing a chart file. The tooltip for `mean_pct` SHALL use `"quantitative"` type, not `"nominal"`.

#### Scenario: All cov columns null
- **WHEN** bam_stats_df has cov_1x_pct, cov_10x_pct, etc. columns but all values are null
- **THEN** plot_bam_coverage_distribution returns None and no HTML/PNG/SVG files are written

#### Scenario: Some cov columns have data
- **WHEN** bam_stats_df has at least one non-null value in cov_*_pct columns
- **THEN** a bar chart is generated with mean_pct on y-axis as quantitative type

### Requirement: BAM metrics sample-wise plot uses grid faceting by metric and set
The `plot_bam_metrics_sample_wise` function SHALL use `alt.Facet("metric:N", columns=3)` for a wrapped grid layout instead of `column=alt.Column("metric:N")` (single horizontal row). The chart SHALL also facet by `set_number` using `.facet(facet=alt.Facet("set_number:N"), columns=2)` when multiple set numbers exist.

#### Scenario: Metrics sample-wise with multiple sets
- **WHEN** bam_stats_df contains samples from 4 different set_numbers and 7 metric columns
- **THEN** the chart produces a grid of metric panels (3 columns) nested within set_number facets (2 columns)

#### Scenario: Metrics sample-wise with single set
- **WHEN** all samples have set_number=0
- **THEN** the chart produces a grid of metric panels (3 columns) without set faceting

### Requirement: BAM metrics bar chart uses step-based width
The `plot_bam_metrics_bars` function (chart 20) SHALL use `width=alt.Step(N)` (where N ≥ 15) instead of a fixed pixel width, so bar width scales with the number of samples and BAM types rather than being squeezed into a fixed area.

#### Scenario: Many samples with 3 BAM types
- **WHEN** 20+ samples each have DN, DT, RT BAM types with xOffset encoding
- **THEN** bars are at least 15px wide and clearly visible

### Requirement: Ti/Tv ratio chart text layer uses chromosome sort
The `plot_ti_tv_ratio` function SHALL pass `sort=chrom_order` to the x-encoding of BOTH the bar layer and the text layer. The text layer SHALL NOT use an unsorted x-encoding.

#### Scenario: Chromosome-ordered Ti/Tv ratio
- **WHEN** the Ti/Tv ratio chart is generated for chromosome-wise grouping
- **THEN** both bar and text layers have chromosomes in natural order (chr1, chr2, ..., chr22, chrX, chrY)
- **AND** text labels are aligned with their corresponding bars

### Requirement: Disease-wise and set-wise plots facet by their wise dimension
All chart functions called from the disease-wise registry SHALL facet by `disease_normalized`. All chart functions called from the set-wise registry SHALL facet by `set_number`. The CLI registry SHALL pass `facet_col="disease_normalized"` or `facet_col="set_number"` respectively. Chart functions that use the wise dimension as their primary x-axis SHALL instead use a secondary dimension (e.g., FILTER, variant_type) on x-axis and the wise dimension as the facet.

#### Scenario: Disease-wise VC distribution with multiple diseases
- **WHEN** 3 disease categories exist and VC distribution is plotted disease-wise
- **THEN** the chart has 3 facet panels (one per disease), each showing FILTER on x-axis

#### Scenario: Set-wise filter distribution with multiple sets
- **WHEN** 4 set_numbers exist and filter distribution is plotted set-wise
- **THEN** the chart has 4 facet panels (one per set), each showing FILTER on x-axis

### Requirement: Tier-axis plots use explicit tier order
All charts with `final_tier` or `caller_tier` on an x-axis, column, or facet SHALL use `sort=FINAL_TIER_ORDER` imported from `tiering_stats.py`. The sort order SHALL be `["C1D1","C1D0","C2D1","C2D0","C3D1","C3D0","C4D1","C4D0","C5D1","C5D0","C6D1","C6D0","C7D1","C7D0"]`.

#### Scenario: Per-tier VAF boxplot with all tiers
- **WHEN** a per-tier VAF boxplot is generated with tiers C1D0 through C6D1
- **THEN** the x-axis shows tiers in FINAL_TIER_ORDER, not alphabetical order

#### Scenario: Tier axis with C10+ tier
- **WHEN** a tier C10D0 exists in the data
- **THEN** it appears after C9D0 in the sort order, not between C1D1 and C2D0
