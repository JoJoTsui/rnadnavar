# variant-visualization (delta)

## MODIFIED Requirements

### Requirement: Cross-modality chart includes percentages
The cross-modality and rescue analysis chart SHALL be replaced by `plot_modality_evidence`, a single stacked bar chart showing the four `modality_evidence` categories (cross_modality, dna_confident, rna_rescued, low_confidence) instead of the current two-panel hconcat of RESCUED and CROSS_MODALITY bar charts. The chart SHALL include count and percentage text labels on each bar segment.

#### Scenario: Modality evidence with percentages
- **WHEN** the modality evidence chart is rendered
- **THEN** a single stacked bar chart is generated showing 4 categories per group with count+percentage annotations
- **AND** no separate CROSS_MODALITY panel exists (redundant with RESCUED)

## REMOVED Requirements

### Requirement: Cross-modality chart includes percentages
**Reason**: Replaced by `plot_modality_evidence` using the four-category `modality_evidence` column. The two-panel RESCUED + CROSS_MODALITY chart was redundant because both columns always carried identical values.
**Migration**: Chart 11 (`plot_cross_modality`) is renamed to `plot_modality_evidence`. Rescue charts 48-56 accept `modality_evidence` as their grouping column instead of `RESCUED`. The `CROSS_MODALITY` column is no longer charted.

## ADDED Requirements

### Requirement: Rescue charts use modality_evidence column
All seven rescue chart functions (Charts 48-56: `plot_rescue_breakdown`, `plot_rescue_by_filter`, `plot_rescue_cross_tab_heatmap`, `plot_rescue_vaf_boxplot`, `plot_rescue_dp_boxplot`, `plot_rescue_sample_distribution`, `plot_rescue_by_tier`, `plot_rescue_rate_trend`, `plot_rescue_caller_support`) SHALL use the `modality_evidence` column as their primary grouping column instead of `RESCUED`. Charts SHALL display all four modality_evidence categories with distinct colors from the MODALITY_EVIDENCE_COLORS registry.

#### Scenario: Rescue breakdown with four categories
- **WHEN** the rescue breakdown chart is rendered
- **THEN** bars show four categories (cross_modality, dna_confident, rna_rescued, low_confidence) per group
- **AND** each category has a distinct color

### Requirement: DNA vs RNA per-caller chart in variant-category wise
The `plot_dna_vs_rna_per_caller` function SHALL be registered in the variant-category wise chart list. When called without explicit `color_col`, the function SHALL use `FILTER` as the color encoding (existing fallback behavior), providing per-caller DNA vs RNA scatter plots colored by variant category within each caller panel.

#### Scenario: Per-caller VAF scatter in variant-category wise
- **WHEN** variant-category wise charts are generated
- **THEN** `plot_dna_vs_rna_per_caller` is called with color_col="FILTER"
- **AND** the chart file `26_dna_vs_rna_per_caller` is saved in `plots/variant-category/`
