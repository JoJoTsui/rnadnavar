# per-caller-vaf-dp-scatters

## Purpose

DNA vs RNA per-caller scatter plots with both VAF and DP metrics, generating a 6-panel layout (3 callers × 2 metrics). Consumes per-caller columns from the normalized parse.

## ADDED Requirements

### Requirement: Per-caller VAF+DP dual scatter generation
The system SHALL extend `plot_dna_vs_rna_per_caller` to accept a `metrics` parameter (default `("VAF", "DP")`) and generate a subchart for each metric per caller pair. For each of the 3 caller pairs (Mutect2, DeepSomatic, Strelka), two scatter plots SHALL be generated: one for `{caller}_VAF` and one for `{caller}_DP`. DP scatters SHALL clip values to [0, 2000] consistent with `plot_dna_vs_rna_dp`.

#### Scenario: VAF+DP scatters generated
- **WHEN** plot_dna_vs_rna_per_caller is called with metrics=("VAF", "DP")
- **THEN** 6 scatter plots are generated (3 callers × 2 metrics)
- **AND** DP scatters have axes limited to [0, 2000]

#### Scenario: VAF-only backward compatibility
- **WHEN** plot_dna_vs_rna_per_caller is called with metrics=("VAF",)
- **THEN** only 3 VAF scatter plots are generated (existing behavior preserved)

### Requirement: Six-panel layout
The system SHALL arrange subcharts in a 2-row × 3-column grid: top row = VAF scatters (Mutect2, DeepSomatic, Strelka), bottom row = DP scatters. Each column SHALL share the caller name label. Color encoding SHALL be shared across all 6 panels via `.resolve_scale(color="shared")`.

#### Scenario: Layout structure
- **WHEN** the per-caller chart is rendered
- **THEN** the layout is 2 rows (VAF top, DP bottom) × 3 columns (callers)
- **AND** color legend is shared across all panels
