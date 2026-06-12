## ADDED Requirements

### Requirement: Per-sample per-tier visualization
The system SHALL generate charts for per-sample per-tier (Level 4) statistics: a multi-panel bar chart showing per-tier variant counts per sample, with samples on y-axis and tier on color.

#### Scenario: Per-sample tier distribution
- **WHEN** Level 4 statistics are computed
- **THEN** a stacked bar chart is generated showing per-tier variant distribution for each sample

### Requirement: Per-tier comparison across samples
The system SHALL generate faceted charts comparing per-tier metrics across all samples, enabling identification of outlier samples within each tier.

#### Scenario: Per-tier VAF comparison
- **WHEN** Level 4 statistics exist
- **THEN** a boxplot is generated showing per-tier DNA_VAF_mean distribution across samples

### Requirement: BAM validation dashboard section
The dashboard SHALL include a dedicated BAM validation section with DP scatter plots, VAF scatter plots, coverage distributions, and the validation heatmap. This section SHALL appear after the tier-aware charts and before the per-sample charts.

#### Scenario: Dashboard with BAM section
- **WHEN** BAM pileup is enabled and validation data exists
- **THEN** dashboard.html includes a "BAM Statistics & Validation" section header separating BAM charts from variant charts

### Requirement: Consolidated dashboard with 25+ charts
The dashboard SHALL contain all charts organized in sections: Overview, Tier Analysis, BAM Statistics & Validation, Per-Sample. Each section SHALL have a header. Charts within sections SHALL be ordered consistently.

#### Scenario: Complete dashboard
- **WHEN** all processing completes
- **THEN** dashboard.html contains 25+ charts in 4 sections with section headers
