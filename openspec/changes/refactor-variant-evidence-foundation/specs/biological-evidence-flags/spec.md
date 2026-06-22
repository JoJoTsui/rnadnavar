# biological-evidence-flags

## Purpose

Biologically-grounded variant quality flags that distinguish technical artifacts from real biology using physical constraints (VAF sum ≤ 1.0), expected VAF ranges per FILTER category, RNA caller agreement (not RNA VAF), and multi-allelic co-occurrence patterns. Flags annotate variants for downstream review without automatically filtering them.

## ADDED Requirements

### Requirement: VAF physics violation flag
The system SHALL set `flag_vaf_overflow=True` for all alleles at a multi-variant position where the sum of DNA_VAF_mean across alleles exceeds 1.1. This represents a physical impossibility at a diploid autosomal site — the VAF values are inconsistent with each other due to decomposition error or AD misallocation.

#### Scenario: VAF physics violation detected
- **WHEN** a (CHROM, POS) has two alleles with VAF=0.60 and VAF=0.55
- **THEN** flag_vaf_overflow=True for both alleles (sum=1.15 > 1.1)

#### Scenario: Valid VAF sum
- **WHEN** a (CHROM, POS) has two alleles with VAF=0.30 and VAF=0.20
- **THEN** flag_vaf_overflow=False (sum=0.50 ≤ 1.0)

### Requirement: Multi-allelic heterogeneity flag
The system SHALL set `flag_multi_allelic_heterogeneity=True` for all alleles at positions classified as true multi-allelic biology (two or more alleles with non-zero VAF, non-zero DP, and different ALT bases). These sites are potential subclonal architecture signals.

#### Scenario: Heterogeneity flag set
- **WHEN** chr7:100958135 has T>C (VAF=0.438, DP=64) and T>G (VAF=0.111, DP=27), classified as true multi-allelic
- **THEN** flag_multi_allelic_heterogeneity=True for both alleles
- **AND** allele_balance_ratio=0.438/(0.438+0.111)=0.798

### Requirement: Category conflict flag
The system SHALL set `flag_category_conflict=True` for all alleles at a multi-variant position where the unique FILTER categories among alleles exceed 1 AND both (or all) categories have biological significance (not just NoConsensus/Artifact). Resolution rules SHALL be recorded in `category_conflict_resolution`.

#### Scenario: Germline+NoConsensus conflict (keep Germline)
- **WHEN** chr12:95287821 has G>A (FILTER=Germline) and G>T (FILTER=NoConsensus)
- **THEN** flag_category_conflict=True for both alleles
- **AND** category_conflict_resolution="keep_germline"

#### Scenario: Artifact+NoConsensus conflict (drop both)
- **WHEN** a site has two alleles both with FILTER ∈ {Artifact, NoConsensus}
- **THEN** flag_category_conflict=True but resolution is "drop_both" because neither has biological signal

### Requirement: Germline low VAF flag
The system SHALL set `flag_germline_low_vaf=True` for variants where FILTER="Germline" AND DNA_VAF_mean < 0.10 AND DNA_DP_mean ≥ 10. True heterozygous germline variants are expected at VAF≈0.50. A germline call with VAF < 0.10 is inconsistent with constitutional heterozygosity and may indicate a misclassified somatic variant, low-purity sample, or mosaic variant.

#### Scenario: Suspiciously low germline VAF
- **WHEN** a variant has FILTER=Germline, DNA_VAF_mean=0.03, DNA_DP_mean=50
- **THEN** flag_germline_low_vaf=True

#### Scenario: Valid germline VAF
- **WHEN** a variant has FILTER=Germline, DNA_VAF_mean=0.48, DNA_DP_mean=50
- **THEN** flag_germline_low_vaf=False

### Requirement: Somatic high VAF flag
The system SHALL set `flag_somatic_high_vaf=True` for variants where FILTER="Somatic" AND DNA_VAF_mean > 0.60. In a diploid genome without copy-number alteration, a heterozygous somatic mutation has maximum VAF≈0.50 × purity. VAF > 0.60 may indicate loss of heterozygosity, homozygous somatic mutation, or misclassified germline variant.

#### Scenario: Suspiciously high somatic VAF
- **WHEN** a variant has FILTER=Somatic, DNA_VAF_mean=0.75
- **THEN** flag_somatic_high_vaf=True

### Requirement: Reference with signal flag
The system SHALL set `flag_reference_with_signal=True` for variants where FILTER="Reference" AND DNA_VAF_mean > 0.05. A "Reference" call means no variant was detected, yet measurable alt allele fraction contradicts this. This indicates caller disagreement or a variant at the edge of detection.

#### Scenario: Reference call with alt signal
- **WHEN** a variant has FILTER=Reference, DNA_VAF_mean=0.12, DNA_DP_mean=30
- **THEN** flag_reference_with_signal=True

### Requirement: RNA rescued flag
The system SHALL set `flag_rna_rescued=True` for variants where DNA_VAF_mean < 0.05 AND N_DNA_CALLERS_SUPPORT ≤ 1 AND N_RNA_CALLERS_SUPPORT ≥ 2 AND RNA_DP_mean ≥ 10. RNA caller agreement and adequate RNA depth provide orthogonal biological evidence that is NOT confounded by allele-specific expression (unlike RNA VAF).

#### Scenario: RNA rescues weak DNA evidence
- **WHEN** a variant has DNA_VAF_mean=0.02, N_DNA_CALLERS_SUPPORT=1, N_RNA_CALLERS_SUPPORT=3, RNA_DP_mean=50
- **THEN** flag_rna_rescued=True

#### Scenario: RNA caller support insufficient
- **WHEN** a variant has DNA_VAF_mean=0.02, N_RNA_CALLERS_SUPPORT=1 (only 1 RNA caller)
- **THEN** flag_rna_rescued=False (insufficient RNA caller agreement)

### Requirement: modality_evidence three-way classification
The system SHALL compute a `modality_evidence` column during `process_single_sample()` that classifies each variant into one of four categories based on combined DNA+RNA evidence: `cross_modality` (C1 tier: N_DNA≥2 AND N_RNA≥2, or DNA_DP>30 AND RNA_DP>30), `dna_confident` (C2 tier: N_DNA≥2 AND N_RNA≤1, or DNA_DP>30 AND RNA_DP≤30), `rna_rescued` (C3+C4 tiers: N_RNA≥2 AND N_DNA≤1, or RNA_DP>30 AND DNA_DP≤30), or `low_confidence` (C5-C7 tiers: single or no caller support, or both DP≤30). Two scheme columns SHALL be computed: `modality_evidence_caller` (caller-support-based) and `modality_evidence_dp` (depth-based).

#### Scenario: Cross-modality variant
- **WHEN** a variant has N_DNA_CALLERS_SUPPORT=2, N_RNA_CALLERS_SUPPORT=2
- **THEN** modality_evidence_caller="cross_modality"

#### Scenario: DNA-confident variant
- **WHEN** a variant has N_DNA_CALLERS_SUPPORT=3, N_RNA_CALLERS_SUPPORT=0
- **THEN** modality_evidence_caller="dna_confident"

#### Scenario: RNA-rescued variant (C3)
- **WHEN** a variant has N_DNA_CALLERS_SUPPORT=1, N_RNA_CALLERS_SUPPORT=3
- **THEN** modality_evidence_caller="rna_rescued"

#### Scenario: RNA-rescued variant (C4)
- **WHEN** a variant has N_DNA_CALLERS_SUPPORT=1, N_RNA_CALLERS_SUPPORT=1
- **THEN** modality_evidence_caller="rna_rescued"

#### Scenario: Low-confidence variant
- **WHEN** a variant has N_DNA_CALLERS_SUPPORT=0, N_RNA_CALLERS_SUPPORT=1 (C6)
- **THEN** modality_evidence_caller="low_confidence"

### Requirement: Backward compatible RESCUED/CROSS_MODALITY columns
The system SHALL retain `RESCUED` and `CROSS_MODALITY` as deprecated columns in the per-sample parquet. `RESCUED` SHALL map to `modality_evidence_caller ∈ {"cross_modality", "rna_rescued"}`. A deprecation log message SHALL be emitted once per run. The columns SHALL be excluded from `_CROSS_SAMPLE_COLS` to prevent accidental use in cross-sample aggregation.

#### Scenario: RESCUED column maps from modality_evidence
- **WHEN** modality_evidence_caller="rna_rescued"
- **THEN** RESCUED="YES" (rna_rescued + cross_modality both map to YES)
- **AND** a deprecation warning is logged once per run
