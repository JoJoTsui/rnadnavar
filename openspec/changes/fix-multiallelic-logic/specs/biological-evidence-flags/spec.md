## ADDED Requirements

### Requirement: VAF physics violation flag
The system SHALL set `flag_vaf_overflow=True` for all alleles at a multi-variant position where the sum of DNA_VAF_mean across alleles exceeds 1.1. This represents a physical impossibility at a diploid autosomal site.

#### Scenario: VAF physics violation detected
- **WHEN** a (CHROM, POS) has two alleles with VAF=0.60 and VAF=0.55
- **THEN** flag_vaf_overflow=True for both alleles (sum=1.15 > 1.1)

#### Scenario: Valid VAF sum
- **WHEN** a (CHROM, POS) has two alleles with VAF=0.30 and VAF=0.20
- **THEN** flag_vaf_overflow=False (sum=0.50 ≤ 1.1)

### Requirement: Multi-allelic heterogeneity flag (true_multi_allelic only)
The system SHALL set `flag_multi_allelic_heterogeneity=True` for all alleles ONLY at positions classified as `true_multi_allelic`. The flag SHALL NOT be set for positions classified as `normalization_artifact` or `noise`. These sites represent potential subclonal architecture signals.

#### Scenario: Heterogeneity flag set for true multi-allelic
- **WHEN** chr7:100958135 is classified as true_multi_allelic
- **THEN** flag_multi_allelic_heterogeneity=True for both alleles at that position

#### Scenario: Heterogeneity flag NOT set for normalization artifact
- **WHEN** chr15:34362339 is classified as normalization_artifact
- **THEN** flag_multi_allelic_heterogeneity=False for all alleles at that position

### Requirement: Category conflict flag with biological significance guard
The system SHALL set `flag_category_conflict=True` for all alleles at a multi-variant position where the unique FILTER categories among alleles exceed 1 AND all categories have biological significance. Categories with biological significance are: Somatic, Germline, RNAedit. Categories WITHOUT biological significance (NoConsensus, Artifact, Reference, PASS) SHALL NOT trigger the flag on their own. The system SHALL record `category_conflict_resolution` as one of: `keep_germline`, `keep_somatic`, `keep_rnaedit`, `drop_both`, `keep_higher_tier`.

#### Scenario: Germline + Somatic conflict
- **WHEN** a site has G>A (FILTER=Germline) and G>T (FILTER=Somatic)
- **THEN** flag_category_conflict=True for both alleles
- **AND** category_conflict_resolution="keep_somatic" (somatic is higher biological priority)

#### Scenario: Artifact + NoConsensus no conflict
- **WHEN** a site has two alleles with FILTER=Artifact and FILTER=NoConsensus
- **THEN** flag_category_conflict=False because neither category has biological significance

#### Scenario: Germline + NoConsensus conflict
- **WHEN** chr12:95287821 has G>A (FILTER=Germline) and G>T (FILTER=NoConsensus)
- **THEN** flag_category_conflict=True because Germline has biological significance
- **AND** category_conflict_resolution="keep_germline"

### Requirement: Extreme strand bias flag
The system SHALL set `flag_extreme_strand_bias=True` for alleles where pre-decomposition strand bias data is available AND `strand_balance < 0.10 OR strand_balance > 0.90`. Strand balance is computed as `F1R2_alt / (F1R2_alt + F2R1_alt)`. Extreme strand bias indicates a potential sequencing artifact. The flag SHALL NOT be set when pre-decomposition data is unavailable.

#### Scenario: Extreme strand bias detected
- **WHEN** an allele has F1R2_alt=2, F2R1_alt=48
- **THEN** strand_balance=2/(2+48)=0.04, flag_extreme_strand_bias=True

#### Scenario: Balanced strand
- **WHEN** an allele has F1R2_alt=25, F2R1_alt=25
- **THEN** strand_balance=0.50, flag_extreme_strand_bias=False

### Requirement: Germline low VAF flag
The system SHALL set `flag_germline_low_vaf=True` for variants where FILTER="Germline" AND DNA_VAF_mean < 0.10. The previous DNA_DP_mean ≥ 10 constraint SHALL be removed so that low-coverage false positives are also flagged.

#### Scenario: Suspiciously low germline VAF with adequate depth
- **WHEN** a variant has FILTER=Germline, DNA_VAF_mean=0.03, DNA_DP_mean=50
- **THEN** flag_germline_low_vaf=True

#### Scenario: Suspiciously low germline VAF with low depth
- **WHEN** a variant has FILTER=Germline, DNA_VAF_mean=0.02, DNA_DP_mean=5
- **THEN** flag_germline_low_vaf=True (no DP constraint)

### Requirement: Somatic high VAF flag
The system SHALL set `flag_somatic_high_vaf=True` for variants where FILTER="Somatic" AND DNA_VAF_mean > 0.60.

#### Scenario: Suspiciously high somatic VAF
- **WHEN** a variant has FILTER=Somatic, DNA_VAF_mean=0.75
- **THEN** flag_somatic_high_vaf=True

### Requirement: Reference with signal flag
The system SHALL set `flag_reference_with_signal=True` for variants where FILTER="Reference" AND DNA_VAF_mean > 0.05.

#### Scenario: Reference call with alt signal
- **WHEN** a variant has FILTER=Reference, DNA_VAF_mean=0.12, DNA_DP_mean=30
- **THEN** flag_reference_with_signal=True

### Requirement: RNA rescued flag
The system SHALL set `flag_rna_rescued=True` for variants where DNA_VAF_mean < 0.05 AND N_DNA_CALLERS_SUPPORT ≤ 1 AND N_RNA_CALLERS_SUPPORT ≥ 2 AND RNA_DP_mean ≥ 10.

#### Scenario: RNA rescues weak DNA evidence
- **WHEN** a variant has DNA_VAF_mean=0.02, N_DNA_CALLERS_SUPPORT=1, N_RNA_CALLERS_SUPPORT=3, RNA_DP_mean=50
- **THEN** flag_rna_rescued=True

### Requirement: modality_evidence_caller uses C1-C7 tier mapping
The system SHALL compute `modality_evidence_caller` by mapping `final_tier` to modality categories: C1D0/C1D1 → `cross_modality`, C2D0/C2D1 → `dna_confident`, C3D0/C3D1/C4D0/C4D1 → `rna_rescued`, C5D0/C5D1/C6D0/C6D1/C7D0/C7D1 → `low_confidence`. The system SHALL NOT use raw `N_DNA_CALLERS_SUPPORT >= 1 AND N_RNA_CALLERS_SUPPORT >= 1` for `cross_modality` classification. A second column `modality_evidence_dp` SHALL use depth-based classification: DNA_DP > 30 AND RNA_DP > 30 → `cross_modality`, DNA_DP > 30 AND RNA_DP ≤ 30 → `dna_confident`, RNA_DP > 30 AND DNA_DP ≤ 30 → `rna_rescued`, both ≤ 30 → `low_confidence`.

#### Scenario: Cross-modality variant (C1 tier)
- **WHEN** a variant has final_tier=C1D1 (≥2 DNA + ≥2 RNA concordant callers + database)
- **THEN** modality_evidence_caller="cross_modality"

#### Scenario: DNA-confident variant (C2 tier)
- **WHEN** a variant has final_tier=C2D0 (≥2 DNA concordant callers, no RNA)
- **THEN** modality_evidence_caller="dna_confident"

#### Scenario: Low-confidence variant (single DNA + single RNA)
- **WHEN** a variant has N_DNA_CALLERS_SUPPORT=1, N_RNA_CALLERS_SUPPORT=1, final_tier=C4D0
- **THEN** modality_evidence_caller="rna_rescued" (mapped from C4 tier, NOT "cross_modality")

### Requirement: --exclude-multiallelic-conflict checks flag_category_conflict
The CLI flag `--exclude-multiallelic-conflict` SHALL exclude variants where `flag_category_conflict=True`. It SHALL NOT check `flag_vaf_overflow`.

#### Scenario: Exclude-multiallelic-conflict filters category conflicts
- **WHEN** --exclude-multiallelic-conflict is set and a variant has flag_category_conflict=True
- **THEN** the variant is excluded from the filtered dataset

### Requirement: Cyvcf2 fallback removed from seq2neo stats
The seq2neo statistics subsystem SHALL use the Rust VCF parser as the only backend. The cyvcf2 import and all Python fallback code paths in `caller_parser.py` and `rescue_parser.py` SHALL be removed. The `--vcf-parser` CLI flag SHALL be removed. If the Rust parser (`stats_core.so`) is unavailable, the pipeline SHALL fail with a clear error message indicating the Rust extension is required.

#### Scenario: Rust parser unavailable
- **WHEN** stats_core.so cannot be loaded
- **THEN** the pipeline exits with error: "Rust VCF parser (stats_core) is required but could not be loaded. Ensure the compiled extension is available."

#### Scenario: Cyvcf2 not imported
- **WHEN** caller_parser.py or rescue_parser.py is imported
- **THEN** no `from cyvcf2 import VCF` statement exists in the module
