# unified-variant-filtering

## Purpose

A single LazyFrame filter expression applied at `combined_df` creation point that propagates consistently to ALL downstream consumers — statistics, visualizations, BAM pileup, rescue analytics, and threshold sweeps. Replaces three separate disconnected filtering mechanisms (`--pileup-mode`, `--min-vaf/--min-dp` for filtered parquet, `--exclude-disease`).

## ADDED Requirements

### Requirement: Single filter injection point
The system SHALL apply all variant-level filters as a single polars `.filter()` expression on the `combined_df` LazyFrame immediately after it is created via `pl.scan_parquet()`. Because `combined_df` is lazy, the filter SHALL be pushed into the parquet scan by the polars optimizer with zero materialization overhead. Every downstream function that receives `combined_df` SHALL receive the filtered view.

#### Scenario: Filter propagated to all consumers
- **WHEN** `--include-filters Somatic Germline` is specified
- **THEN** the filter is applied once at combined_df creation
- **AND** statistics, visualizations, threshold sweeps, rescue analytics, and BAM pileup ALL receive the filtered data

### Requirement: Filter expression composition
The system SHALL compose multiple filter criteria into a single boolean expression using `&` (AND) logic. Filters SHALL include: FILTER category inclusion/exclusion (`--include-filters`, `--exclude-filters`), minimum caller support (`--min-dna-callers`, `--min-rna-callers`), VAF/DP thresholds (`--min-vaf`, `--min-dp`), evidence tier minimum (`--min-evidence-tier`), gnomAD population frequency maximum (`--max-gnomad-af`), multi-allelic conflict exclusion (`--exclude-multiallelic-conflict`), and disease exclusion (`--exclude-disease`).

#### Scenario: Multiple filters combined
- **WHEN** `--include-filters Somatic Germline --min-vaf 0.05 --max-gnomad-af 0.01` is specified
- **THEN** a variant passes only if FILTER ∈ {Somatic, Germline} AND (DNA_VAF_mean ≥ 0.05 OR (N_RNA_CALLERS_SUPPORT ≥ 2 AND RNA_DP_mean ≥ 10)) AND (GNOMAD_AF IS NULL OR GNOMAD_AF < 0.01)

#### Scenario: No filters specified (default)
- **WHEN** no filter flags are specified
- **THEN** the filter expression is `pl.lit(True)` (identity — keeps all variants)

### Requirement: RNA VAF is NOT used for filtering
The system SHALL NOT use `RNA_VAF_mean` in the `--min-vaf` filter threshold. RNA VAF is confounded by allele-specific expression (ASE). Instead, the RNA leg of the filter SHALL use `N_RNA_CALLERS_SUPPORT >= 2 AND RNA_DP_mean >= 10` — caller agreement and expression depth, which are not confounded by ASE.

#### Scenario: RNA VAF excluded from VAF filter
- **WHEN** `--min-vaf 0.10` is specified
- **THEN** a variant passes the VAF filter if DNA_VAF_mean ≥ 0.10 OR (N_RNA_CALLERS_SUPPORT ≥ 2 AND RNA_DP_mean ≥ 10)
- **AND** RNA_VAF_mean is NOT checked against the --min-vaf threshold

#### Scenario: No RNA data available
- **WHEN** a variant has no RNA caller data (all RNA columns are null)
- **THEN** only the DNA leg of the filter is evaluated
- **AND** the variant passes if DNA_VAF_mean ≥ threshold

### Requirement: VAF physics constraint
The system SHALL exclude variants where `vaf_sum > 1.1` at a multi-variant position, as this violates the physical constraint that total allele fraction at a diploid site cannot exceed 1.0. This filter SHALL be applied automatically (no CLI flag needed) because it represents a data integrity check, not an analytical choice.

#### Scenario: VAF overflow detected and filtered
- **WHEN** a (CHROM, POS) has alleles with VAF_sum=1.25
- **THEN** all alleles at that position are excluded from downstream analysis
- **AND** the exclusion is logged with the position and VAF_sum value

### Requirement: `--no-filter` escape hatch
The system SHALL provide a `--no-filter` flag that disables ALL filtering (including the VAF physics constraint), producing unfiltered output identical to the current default behavior. This flag SHALL be used for debugging and for comparing filtered vs unfiltered results.

#### Scenario: No-filter mode
- **WHEN** `--no-filter` is specified
- **THEN** no filters are applied, including VAF physics constraint
- **AND** all variants pass through to downstream processing

### Requirement: Deprecated flags emit warnings
The system SHALL accept the legacy `--pileup-mode` flag and map it to the unified filter: `--pileup-mode filtered` SHALL be equivalent to `--exclude-filters NoConsensus`. The system SHALL emit a deprecation warning directing users to use `--exclude-filters` instead.

#### Scenario: Legacy --pileup-mode flag
- **WHEN** `--pileup-mode filtered` is specified
- **THEN** a deprecation warning is printed: "`--pileup-mode` is deprecated. Use `--exclude-filters NoConsensus` instead."
- **AND** NoConsensus variants are excluded from ALL downstream processing, not just BAM pileup
