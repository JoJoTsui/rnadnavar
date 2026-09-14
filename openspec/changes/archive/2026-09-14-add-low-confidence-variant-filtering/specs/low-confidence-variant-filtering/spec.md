## ADDED Requirements

### Requirement: Cross-sample recurrence is computed
The system SHALL compute `n_recurrent_samples` for each variant as the number of distinct samples that contain the same (CHROM, POS, REF, ALT) variant. Computation SHALL use `pl.len().over(['CHROM','POS','REF','ALT'])` on `combined_df` after all samples are scanned. The column SHALL be added to the per-variant parquet and included in `_CROSS_SAMPLE_COLS`.

#### Scenario: Variant unique to one sample
- **WHEN** a Somatic variant at chr7:140453136 A>T appears in only 1 of 61 samples
- **THEN** n_recurrent_samples=1

#### Scenario: Recurrent artifact across samples
- **WHEN** a variant at chr17:17383677 A>G appears in 51 of 61 samples
- **THEN** n_recurrent_samples=51
- **AND** the variant is flagged for review as a likely mapping artifact

### Requirement: Strand bias Fisher exact test is computed
The system SHALL compute `max_strand_bias_fisher_p` as the minimum Fisher exact p-value across all callers with F1R2/F2R1 strand data. For each caller, a 2×2 contingency table is constructed: `[[F1R2_ref, F2R1_ref], [F1R2_alt, F2R1_alt]]`. The Fisher exact p-value SHALL be computed using `scipy.stats.fisher_exact` (or a pure-Python fallback). The test SHALL only be computed for variants with alt_reads ≥ 10; otherwise the column is null.

#### Scenario: Extreme strand bias
- **WHEN** a variant has Mutect2 F1R2_alt=2, F2R1_alt=48, F1R2_ref=50, F2R1_ref=50
- **THEN** max_strand_bias_fisher_p < 1e-3 (extreme bias in alt strand)

#### Scenario: Balanced strand
- **WHEN** a variant has F1R2_alt=25, F2R1_alt=25, F1R2_ref=50, F2R1_ref=50
- **THEN** max_strand_bias_fisher_p > 0.05 (no significant bias)

#### Scenario: Insufficient alt reads
- **WHEN** a variant has alt_reads=3 (below threshold of 10)
- **THEN** max_strand_bias_fisher_p=null (not computed)

### Requirement: Stage 0 hard exclusions drop false variants
The system SHALL apply hard exclusion filters that auto-drop variants failing physical, biological, or evidence constraints. Hard exclusions SHALL include:
- `N_SUPPORT_CALLERS == 0` (no caller support, C7 tier)
- `flag_vaf_overflow` (VAF sum > 1.1, physics violation)
- `multiallelic_class == "noise"` (one real allele + noise, after corrected classification)
- `DNA_DP_mean < 5 AND RNA_DP_mean < 5` (no modality has coverage)
- `DNA_ALT_DP_mean < 2 AND RNA_ALT_DP_mean < 2` (no alt evidence in either modality)
- `FILTER == "Germline" AND flag_germline_low_vaf` (Germline with VAF < 0.10)
- `FILTER == "Somatic" AND flag_somatic_high_vaf AND DNA_REF_DP_mean < 0.1 * DNA_DP_mean` (Somatic with LOH-like VAF)
- `FILTER == "Reference" AND DNA_VAF_mean > 0.10 AND DNA_ALT_DP_mean >= 3` (Reference with substantive alt signal)

Dropped variants SHALL be written to `failed.tsv` with a `drop_reason` column. The `--no-hard-filter` CLI flag SHALL disable all hard exclusions.

#### Scenario: Reference with signal is hard-dropped
- **WHEN** a variant has FILTER=Reference, DNA_VAF_mean=0.15, DNA_ALT_DP_mean=12
- **THEN** the variant is dropped with drop_reason="reference_with_signal"

#### Scenario: Germline with low VAF is hard-dropped
- **WHEN** a variant has FILTER=Germline, DNA_VAF_mean=0.03
- **THEN** the variant is dropped with drop_reason="germline_low_vaf"

#### Scenario: RNA-rescued variant with low DNA depth is NOT dropped
- **WHEN** a variant has DNA_DP_mean=3, RNA_DP_mean=50, DNA_ALT_DP_mean=1, RNA_ALT_DP_mean=15
- **THEN** the variant is NOT dropped (RNA has coverage and alt evidence)

#### Scenario: No-hard-filter escape hatch
- **WHEN** --no-hard-filter is set
- **THEN** no hard exclusions are applied and all variants pass to Stage 1

### Requirement: Stage 1 soft flags mark variants for review
The system SHALL apply soft flags that mark variants for review without dropping them. Soft flags SHALL include:
- `flag_category_conflict`
- `flag_multi_allelic_heterogeneity`
- `multiallelic_class == "normalization_artifact"`
- `flag_rna_rescued`
- `FILTER == "Germline" AND DNA_VAF_mean > 0.85` (LOH contamination candidate)
- `modality_evidence_caller == "low_confidence"`
- `n_recurrent_samples > 20` (pan-sample recurrent, likely artifact)
- `N_SUPPORT_CALLERS == 1 AND FILTER in {"Somatic", "Germline"}` (single-caller only)
- `max_strand_bias_fisher_p < 1e-3 AND alt_reads >= 10` (extreme strand bias)
- `BAM_RT_mean_MQ < 2 AND RNA_DP_mean >= 20` (low RNA MAPQ, after MAPQ=255 fix)

Soft-flagged variants SHALL have a `soft_flags` column (comma-separated flag names) in the output parquet.

#### Scenario: Single-caller Somatic is soft-flagged
- **WHEN** a variant has FILTER=Somatic, N_SUPPORT_CALLERS=1
- **THEN** soft_flags includes "single_caller_only"

#### Scenario: High-recurrence variant is soft-flagged
- **WHEN** a variant has n_recurrent_samples=45
- **THEN** soft_flags includes "high_recurrence"

### Requirement: Stage 2 confidence tiers classify variants
The system SHALL assign each variant (that passes Stage 0) to a confidence tier:
- **HIGH**: no soft flags AND `final_tier ∈ {C1D0, C1D1, C2D1, C3D1, C4D1}` AND (`N_DNA_CALLERS_SUPPORT ≥ 1 OR N_RNA_CALLERS_SUPPORT ≥ 2`)
- **MEDIUM**: no soft flags AND `final_tier ∈ {C2D0, C3D0, C4D0}`
- **LOW**: passes Stage 0 but has any soft flag

The `confidence_tier` column SHALL be added to the output parquet. A `confidence_tier_summary.tsv` SHALL be generated showing counts per tier × FILTER.

#### Scenario: High-confidence somatic
- **WHEN** a variant has FILTER=Somatic, final_tier=C2D1, no soft flags, N_DNA_CALLERS_SUPPORT=3
- **THEN** confidence_tier="HIGH"

#### Scenario: Low-confidence reference
- **WHEN** a variant has FILTER=Reference, final_tier=C2D0, soft_flags="reference_with_signal"
- **THEN** confidence_tier="LOW"

### Requirement: Both filtered and unfiltered statistics are emitted
When the unified filter is active, the system SHALL compute and write both filtered and unfiltered versions of major statistics. Unfiltered versions SHALL be computed from `combined_df` before `.filter()`. Filtered versions SHALL be computed from the filtered `combined_df`. Output files:
- `sample_summary.tsv` (unfiltered) + `sample_summary_filtered.tsv` (filtered)
- `set_summary.tsv` (unfiltered) + `set_summary_filtered.tsv` (filtered)
- `dataset_summary.tsv` (filtered) + `dataset_summary_unfiltered.tsv` (unfiltered)

#### Scenario: Filtered run emits both versions
- **WHEN** the pipeline runs with --exclude-filters NoConsensus
- **THEN** both sample_summary.tsv (107.9M variants) and sample_summary_filtered.tsv (9.75M variants) are generated
- **AND** the file names clearly indicate which is filtered

### Requirement: Caller overlap matrix and support distribution are emitted
The system SHALL compute and write:
- `caller_overlap_matrix.tsv`: a symmetric N×N matrix (N = number of callers) where cell (i,j) = count of variants where both caller_i and caller_j called the variant. Diagonal = total calls by that caller.
- `caller_support_distribution.tsv`: histogram of N_SUPPORT_CALLERS (count of variants called by exactly 1, 2, ..., 6 callers).
- `tier_distribution.tsv`: the existing tier count distribution (renamed from the misnamed `caller_overlap.tsv`).

#### Scenario: Caller overlap matrix
- **WHEN** 6 callers are configured (DNA_mutect2, DNA_strelka, DNA_deepsomatic, RNA_mutect2, RNA_strelka, RNA_deepsomatic)
- **THEN** caller_overlap_matrix.tsv is a 6×6 symmetric matrix with caller names as row/column headers

#### Scenario: Support distribution
- **WHEN** 9.75M filtered variants have N_SUPPORT_CALLERS ranging from 0 to 6
- **THEN** caller_support_distribution.tsv has one row per support count with the number of variants
