# variant-visualization

## Purpose

Interactive altair-based variant visualization dashboard with static PNG/SVG exports for the seq2neo variant statistics pipeline.

## Requirements

### Requirement: Dashboard HTML renders all charts
The system SHALL generate a dashboard.html that renders all generated charts. Charts SHALL be embedded as body-content fragments extracted from altair's `to_html()` output, assembled into a single valid HTML document with shared vega-embed CDN loaded once in the `<head>`.

#### Scenario: Dashboard with multiple charts
- **WHEN** 12 charts are generated and passed to generate_dashboard()
- **THEN** dashboard.html contains all 12 charts rendered correctly in a browser

### Requirement: VC distribution chart shows color legend
The VC distribution chart SHALL display a visible color legend mapping VC categories (Somatic, Germline, Reference, Artifact) to their respective colors.

#### Scenario: VC distribution with legend
- **WHEN** the VC distribution chart is rendered
- **THEN** a color legend with labels Somatic, Germline, Reference, and Artifact is visible

### Requirement: VAF distribution uses violin plots per tier
The VAF distribution chart SHALL use violin-style plots (layered density area charts) instead of box plots, with facets for each CxDy caller tier. DNA and RNA callers SHALL be shown in separate facets.

#### Scenario: VAF violin plot
- **WHEN** variant data contains per-caller VAF and tier assignments
- **THEN** a multi-panel violin plot is generated showing VAF distributions grouped by caller and faceted by tier

### Requirement: Per-sample variant chart uses sample ID on x-axis
The per-sample variant count chart SHALL use sample ID on the x-axis (or y-axis for horizontal bars), ordered by variant count descending, replacing the current disease-based grouping.

#### Scenario: Per-sample variant chart with sample IDs
- **WHEN** 65 samples are processed
- **THEN** the per-sample variant chart shows 65 bars with sample ID labels

### Requirement: Cross-modality chart includes percentages
The cross-modality and rescue analysis chart SHALL include percentage annotations showing the proportion of cross-modality and rescued variants relative to total variants per set.

#### Scenario: Cross-modality with percentages
- **WHEN** the cross-modality chart is rendered
- **THEN** each bar includes a text label showing the percentage relative to the set's total variant count

### Requirement: Per-sample variant distribution plot
The system SHALL generate a per-sample variant count bar chart showing each sample's total variant count, ordered by descending count.

#### Scenario: Per-sample distribution
- **WHEN** all samples are processed
- **THEN** a bar chart is generated with sample_id on y-axis and total_variants on x-axis, sorted by count

### Requirement: Per-tier VAF violin plot
The system SHALL generate violin plots showing VAF distribution per caller, with facets for each CxDy caller tier. DNA caller VAFs and RNA caller VAFs SHALL be plotted in separate rows.

#### Scenario: Per-tier VAF violin
- **WHEN** variant data has tier assignments and per-caller VAF columns
- **THEN** a faceted violin plot is generated with rows for DNA/RNA and columns for C1-C7 tiers, showing VAF distribution per caller

### Requirement: Per-tier DP violin plot
The system SHALL generate violin plots showing DP distribution per caller, faceted by CxDy caller tier.

#### Scenario: Per-tier DP violin
- **WHEN** variant data has tier assignments and per-caller DP columns
- **THEN** a faceted violin plot is generated showing DP distribution per caller per tier

### Requirement: Per-tier GT concordance chart
The system SHALL generate GT concordance breakdown charts faceted by CxDy caller tier, showing the number of variants where 2, 3, or 4 GT-bearing callers agree.

#### Scenario: Per-tier GT concordance
- **WHEN** variants have GT fields from Mutect2 and DeepSomatic callers and tier assignments
- **THEN** a faceted bar chart is generated with one panel per tier, showing GT agreement counts

### Requirement: Per-tier caller overlap chart
The system SHALL generate caller support overlap histograms (N_SUPPORT_CALLERS distribution) faceted by CxDy tier.

#### Scenario: Per-tier caller overlap
- **WHEN** variants have N_SUPPORT_CALLERS and tier assignments
- **THEN** a faceted histogram is generated showing caller support distribution per tier

### Requirement: REF_DP / ALT_DP scatter plots
The system SHALL generate DNA vs RNA REF_DP and ALT_DP scatter plots, colored by VC or tier.

#### Scenario: REF_DP scatter
- **WHEN** variant data has DNA_REF_DP_mean and RNA_REF_DP_mean columns
- **THEN** a scatter plot of DNA vs RNA mean REF_DP is generated

### Requirement: Per-tier variant type chart
The system SHALL generate variant type distribution charts (SNV/INS/DEL/MNV) faceted by CxDy tier.

#### Scenario: Per-tier variant types
- **WHEN** variants have variant_type and tier assignments
- **THEN** a stacked bar chart is generated showing variant type distribution per tier
