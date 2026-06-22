# multi-allelic-visualization

## Purpose

Charts for multi-allelic site analysis consuming columns from the refactor-variant-evidence-foundation change: `multiallelic_class`, `vaf_sum`, `allele_balance_ratio`, `category_conflict`, `n_alleles_at_site`.

## ADDED Requirements

### Requirement: Multi-allelic classification distribution chart
The system SHALL generate `plot_multiallelic_classification`, a grouped bar chart showing the count of multi-allelic sites per `multiallelic_class` category (normalization_artifact, noise, true_multi_allelic) grouped by set_number or disease. The chart SHALL include count text labels on each bar using `_make_bar_text`.

#### Scenario: Multi-allelic classification per set
- **WHEN** combined_df contains multiallelic_class and set_number columns
- **THEN** a bar chart is generated showing 3 category bars per set with count labels

### Requirement: Allele balance scatter plot
The system SHALL generate `plot_allele_balance_scatter`, a scatter plot of major allele VAF vs minor allele VAF at multi-allelic sites where `n_alleles_at_site >= 2`. Points SHALL be colored by `multiallelic_class`. The chart SHALL include a diagonal reference line (y=x) showing the line of equal VAF.

#### Scenario: Allele balance scatter
- **WHEN** combined_df contains multi-allelic sites with 2+ alleles
- **THEN** a scatter plot is generated with major VAF on x-axis, minor VAF on y-axis, colored by classification

### Requirement: VAF sum histogram
The system SHALL generate `plot_vaf_sum_histogram`, a histogram of `vaf_sum` values across all multi-allelic sites. A vertical reference line SHALL be drawn at vaf_sum=1.0. Bins above 1.0 SHALL be colored distinctly to highlight physics violations.

#### Scenario: VAF sum histogram
- **WHEN** multi-allelic sites have vaf_sum values
- **THEN** a histogram is generated with reference line at 1.0 and distinct color for >1.0 bins

### Requirement: Category conflict summary chart
The system SHALL generate `plot_category_conflict_summary`, a horizontal bar chart showing the count of multi-allelic sites per conflict type (unique combination of FILTER categories). The top 10 conflict types SHALL be shown, ordered by count descending.

#### Scenario: Category conflict summary
- **WHEN** category_conflict=True for some multi-allelic sites
- **THEN** a bar chart shows the count per conflict type ordered by frequency
