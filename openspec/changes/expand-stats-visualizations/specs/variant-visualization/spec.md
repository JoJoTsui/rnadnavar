# variant-visualization (delta)

## MODIFIED Requirements

### Requirement: Rescue charts use modality_evidence column
All rescue chart functions (Charts 48-56) SHALL use the `modality_evidence_caller` column as their primary grouping column via an `evidence_col` parameter. When `modality_evidence_caller` is not available in the DataFrame, the functions SHALL fall back to `RESCUED` for backward compatibility with older parquet files.

#### Scenario: Rescue breakdown with modality_evidence
- **WHEN** the rescue breakdown chart is rendered with modality_evidence_caller available
- **THEN** bars show four categories (cross_modality, dna_confident, rna_rescued, low_confidence) per group with distinct colors from MODALITY_EVIDENCE_COLORS (#2ca02c, #1f77b4, #ff7f0e, #d62728)

#### Scenario: Fallback to RESCUED column
- **WHEN** the rescue breakdown chart is rendered with old parquet lacking modality_evidence_caller
- **THEN** bars show two categories (YES, NO) from the RESCUED column with legacy colors

### Requirement: Cross-modality chart redesigned
The `plot_cross_modality` function SHALL be redesigned as `plot_modality_evidence`, producing a single 4-category stacked bar chart per group instead of the previous two-panel hconcat of RESCUED and CROSS_MODALITY. The chart SHALL be saved as `11_modality_evidence` (replacing `11_cross_modality`).

#### Scenario: Modality evidence chart
- **WHEN** the chart is rendered
- **THEN** a single stacked bar chart with 4 categories is generated
- **AND** the file is saved as 11_modality_evidence.{svg,png,html}

## REMOVED Requirements

### Requirement: Cross-modality chart includes percentages
**Reason**: Replaced by `plot_modality_evidence` in Proposal 1 (refactor-variant-evidence-foundation). The two-panel hconcat was redundant because both panels showed identical data.
**Migration**: Chart 11 now shows modality_evidence 4-category breakdown. The old CROSS_MODALITY column is no longer charted. The RESCUED column is retained in parquet for backward compatibility but not used by new charts.

## ADDED Requirements

### Requirement: Multi-allelic charts registered in chart pipeline
The system SHALL register `plot_multiallelic_classification`, `plot_allele_balance_scatter`, `plot_vaf_sum_histogram`, and `plot_category_conflict_summary` in the global chart list (not per-wise), generating output in `plots/multi_allelic/`.

#### Scenario: Multi-allelic charts generated
- **WHEN** the full chart pipeline runs
- **THEN** four multi-allelic chart files are created in plots/multi_allelic/

### Requirement: BAM charts registered in chart pipeline
The system SHALL register `plot_bam_metrics_sample_wise` and `plot_bam_coverage_distribution` in the global chart list, generating output in `plots/bam/`.

#### Scenario: BAM charts generated
- **WHEN** bam_stats_df is available
- **THEN** two BAM chart files are created in plots/bam/
