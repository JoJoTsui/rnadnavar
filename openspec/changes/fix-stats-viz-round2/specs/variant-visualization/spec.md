## MODIFIED Requirements

### Requirement: BAM coverage chart uses correct data source
The BAM coverage violin chart (`plot_bam_coverage_violin`) SHALL source its data from BAM pileup columns (`BAM_DT_DP`, `BAM_DT_REF_DP`, `BAM_DT_ALT_DP`, `BAM_RT_DP`, etc.) present in the per-sample variant parquet files. The chart SHALL handle partial BAM type availability (e.g., DT + RT present but DN absent) by skipping missing BAM types rather than returning early. The chart SHALL validate that sampled data is non-empty before calling `transform_density`.

#### Scenario: BAM coverage violin with partial BAM types
- **WHEN** only DT and RT BAM pileup columns exist (no DN)
- **THEN** the chart SHALL generate a violin plot for DT and RT metrics only
- **AND** SHALL NOT silently return due to `len(dp_cols) < 2` check

#### Scenario: BAM coverage violin with empty sampled data
- **WHEN** sampled data after `_sample_if_large` is empty
- **THEN** the chart SHALL return gracefully without calling `transform_density`

## ADDED Requirements

### Requirement: Chromosome ordering in all CHROM-grouped charts
All charts that group or sort by CHROM SHALL use `_sort_chromosomes()` to ensure natural chromosome ordering (chr1..22, chrX, chrY, chrM). This SHALL apply to: `plot_vc_distribution` (when group_col=CHROM), `plot_dna_vs_rna_per_caller` (when group_col=CHROM), and any chart using `_CHROMOSOME_ORDER` for sorting.

#### Scenario: VC distribution by chromosome
- **WHEN** `plot_vc_distribution` groups by CHROM with group_col="CHROM"
- **THEN** chromosomes SHALL appear in natural order (chr1, chr2, ..., chrX, chrY)
- **AND** chr10 SHALL appear after chr9, not after chr1

### Requirement: VAF distribution plots use consistent styling
VAF distribution plots SHALL use a consistent color palette across all VAF charts. The VAF axis SHALL be clamped to `[0, 1]` for visualization purposes (Strelka VAF > 1 SHALL be annotated separately). Vertical reference lines SHALL be drawn at VAF thresholds 0.005 and 0.01.

#### Scenario: VAF distribution with clamped axis
- **WHEN** Strelka variants have VAF > 1.0
- **THEN** the plot x-axis SHALL extend to at most 1.0
- **AND** an annotation SHALL note that Strelka VAF uses tier-1 depth denominator

### Requirement: Sample-wise charts faceted by set_number
Sample-wise bar charts (`plot_bam_metrics_bars`, `plot_per_sample_tier_distribution`) SHALL facet or group by `set_number` when the data includes multiple sets. All samples SHALL be shown (no top-N filtering). Charts SHALL use shared axis scales across facets for cross-set comparison.

#### Scenario: BAM metrics bars with 4 sets
- **WHEN** BAM stats contain samples from 4 sets
- **THEN** the chart SHALL show all samples in a faceted grid by set_number
- **AND** y-axis scales SHALL be shared across all facets

### Requirement: Database enrichment chart excludes D=0 tiers
The `plot_database_enrichment_by_tier` chart SHALL exclude CxD0 tiers (tiers ending in D0) from the display, or SHALL annotate them with a note explaining that D=0 means "no database evidence" by design.

#### Scenario: Database enrichment without D0 tiers
- **WHEN** variant data has COSMIC_ID and GNOMAD_AF columns
- **THEN** the enrichment chart SHALL either skip D=0 tiers or include them with a "no DB evidence by design" annotation
