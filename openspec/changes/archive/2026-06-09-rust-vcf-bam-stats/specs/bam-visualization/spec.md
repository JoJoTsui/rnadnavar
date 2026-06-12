## ADDED Requirements

### Requirement: Per-sample BAM metrics bar chart
The system SHALL generate a grouped bar chart showing per-sample BAM metrics (total reads, mapped reads, mapping rate) for all three BAM types (DN, DT, RT) side by side. The chart SHALL order samples by total read count descending.

#### Scenario: BAM metrics for all samples
- **WHEN** BAM stats exist for 65 samples
- **THEN** a grouped bar chart with 65 sample groups, each with 3 bars (DN/DT/RT) is generated

### Requirement: BAM DP vs Caller DP scatter plot
The system SHALL generate scatter plots comparing BAM pileup DP against caller VCF DP for each of the 6 callers. Points SHALL be colored by CxDy caller tier.

#### Scenario: DNA Mutect2 DP vs BAM DP scatter
- **WHEN** BAM pileup DP and DNA_mutect2_DP columns exist
- **THEN** a scatter plot is generated with correlation line, colored by caller tier

### Requirement: BAM coverage distribution violin plot
The system SHALL generate violin plots showing the distribution of BAM coverage (DP at variant positions) for each BAM type (DN, DT, RT), split by CxDy caller tier.

#### Scenario: Coverage violin per tier
- **WHEN** BAM pileup data exists with tier assignments
- **THEN** violin plots are generated showing DP distribution per BAM type, faceted by caller tier

### Requirement: BAM validation heatmap
The system SHALL generate a heatmap showing per-sample BAM-vs-caller validation mismatch rates for each metric (DP, VAF) and each caller, enabling rapid identification of problematic samples.

#### Scenario: Validation heatmap
- **WHEN** cross-validation results exist for all samples
- **THEN** a heatmap is generated with samples on y-axis, (caller, metric) pairs on x-axis, colored by mismatch percentage
