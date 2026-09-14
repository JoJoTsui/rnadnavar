## ADDED Requirements

### Requirement: Confidence tier assignment for all variants
The system SHALL assign a confidence tier (HIGH, MEDIUM, or LOW) to every variant using the `compute_confidence_tier()` function. HIGH confidence SHALL require: no soft flags (category_conflict, multi_allelic_heterogeneity, rna_rescued, germline_high_vaf, low_rna_mapq, low_confidence modality, high recurrence, single-caller Somatic/Germline), final_tier in {C1D0, C1D1, C2D1, C3D1, C4D1}, and at least 1 DNA caller or 2 RNA callers supporting. MEDIUM confidence SHALL require: no soft flags AND final_tier in {C2D0, C3D0, C4D0}. LOW confidence SHALL be assigned to all variants with any soft flag.

#### Scenario: HIGH confidence variant
- **WHEN** a variant has final_tier=C1D1, no soft flags, and 2 DNA callers supporting
- **THEN** the variant is assigned confidence_tier=HIGH

#### Scenario: LOW confidence variant with category conflict
- **WHEN** a variant has flag_category_conflict=True
- **THEN** the variant is assigned confidence_tier=LOW regardless of tier or caller support

### Requirement: Confidence wise dimension
The system SHALL include `confidence` as a wise dimension alongside set, disease, sample, tier, caller, chromosome, and variant-category. The `compute_wise_summary()` function SHALL accept `["confidence_tier"]` as group_cols and produce a `confidence_summary.tsv` in `stats/confidence/` with all 27 shared wise metrics.

#### Scenario: Confidence summary generation
- **WHEN** `--wise` includes confidence (or defaults to all wises)
- **THEN** `stats/confidence/confidence_summary.tsv` is written with one row per confidence tier (HIGH/MEDIUM/LOW) containing n_variants, n_somatic, mean_dna_vaf, mean_rna_vaf, and other shared metrics

### Requirement: Confidence x FILTER cross-tabulation
The system SHALL compute a cross-tabulation of confidence_tier by FILTER, writing `confidence_filter_breakdown.tsv` to `stats/confidence/`. This TSV SHALL contain columns for confidence_tier, FILTER, and count.

#### Scenario: Confidence x FILTER breakdown
- **WHEN** variants have both confidence_tier and FILTER columns
- **THEN** a TSV is generated with rows like (HIGH, Somatic, N), (HIGH, Germline, M), (MEDIUM, Somatic, P), etc.

### Requirement: Confidence-wise chart registry
The system SHALL include a `_WISE_CHART_REGISTRY["confidence"]` entry with at minimum: VC distribution per confidence, variant type distribution per confidence, VAF distribution per confidence, tier distribution per confidence, and FILTER distribution per confidence. Charts SHALL be saved to `plots/confidence/`.

#### Scenario: Confidence charts generated
- **WHEN** the confidence wise dimension is active
- **THEN** at least 5 chart HTML files are written to `plots/confidence/`
- **AND** each chart uses `group_col="confidence_tier"` or `color_col="confidence_tier"`

### Requirement: --min-confidence-tier unified filter
The system SHALL support a `--min-confidence-tier` CLI argument accepting values HIGH or MEDIUM. When specified, the unified filter pipeline SHALL exclude variants with confidence_tier below the specified threshold. The filtering SHALL use the ordered semantics: HIGH > MEDIUM > LOW.

#### Scenario: Filter to HIGH confidence only
- **WHEN** `--min-confidence-tier HIGH` is specified
- **THEN** only variants with confidence_tier=HIGH are included in all statistics, charts, and outputs
- **AND** MEDIUM and LOW confidence variants are excluded

#### Scenario: No confidence filter
- **WHEN** `--min-confidence-tier` is not specified
- **THEN** all confidence tiers are included in outputs

### Requirement: --export-high-confidence ML data export
The system SHALL support a `--export-high-confidence` CLI flag that writes `high_confidence_variants.parquet` to `stats/confidence/` containing only variants with confidence_tier=HIGH and FILTER in {Somatic, Germline, Reference}. The output SHALL include all variant columns from `combined_df` suitable for downstream model training.

#### Scenario: Export high-confidence variants
- **WHEN** `--export-high-confidence` is specified
- **THEN** `stats/confidence/high_confidence_variants.parquet` is written
- **AND** every variant in the output has confidence_tier=HIGH
- **AND** every variant has FILTER in {Somatic, Germline, Reference}

#### Scenario: Export with filtering combined
- **WHEN** both `--export-high-confidence` and `--min-confidence-tier HIGH` are specified
- **THEN** the exported parquet contains the same variants as the filtered combined_df (no additional filtering)

### Requirement: Confidence documentation in helper text
The `--help` output SHALL document the confidence tier definitions: HIGH (no soft flags, top-tier caller support, good database evidence), MEDIUM (no soft flags, medium-tier), LOW (any soft flag). The `--min-confidence-tier` and `--export-high-confidence` flags SHALL be documented with their semantics.
