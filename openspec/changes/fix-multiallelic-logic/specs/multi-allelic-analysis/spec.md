## ADDED Requirements

### Requirement: Multi-variant positions are detected with per-allele aggregation
The system SHALL detect positions where more than one variant record exists at the same (CHROM, POS) using polars `group_by(["CHROM", "POS"])`. For each multi-variant position, the system SHALL aggregate per-allele data: FILTER category, ALT allele sequence, REF allele sequence, DNA_VAF_mean, RNA_VAF_mean, DNA_DP_mean, DNA_ALT_DP_mean, per-caller support counts, and pre-decomposition allele registry data (`pre_norm_*_af`, `pre_norm_*_f1r2`, `pre_norm_*_f2r1`, `pre_norm_*_gt_alleles`) when available.

#### Scenario: Position with two alleles detected
- **WHEN** chr7:100958135 has T>C (VAF=0.438, DP=64) and T>G (VAF=0.111, DP=27)
- **THEN** n_alleles_at_site=2, vaf_sum=0.549, both alleles' per-allele data is aggregated

#### Scenario: Position with single allele
- **WHEN** a (CHROM, POS) has exactly one variant record
- **THEN** n_alleles_at_site=1, multiallelic_class="single", no multi-allelic metrics computed

### Requirement: Multi-variant positions classified using per-allele criteria
The system SHALL classify each multi-variant position into one of four categories using per-allele (not aggregate) criteria:

- **normalization_artifact**: Alleles at the position have different REF/ALT string lengths (e.g., REF="AT" ALT="A" vs REF="A" ALT="T" at same position) OR at least one allele has VAF ≤ 0.001 AND DP ≤ 1 (zero signal). These represent the same biological event encoded differently after vt decompose.
- **noise**: Exactly one allele has VAF ≥ 0.01 AND DP ≥ 5; all other alleles have VAF < 0.01 OR DP < 5. The signal allele is the biological variant; others are technical noise.
- **true_multi_allelic**: Two or more alleles have VAF ≥ 0.01 AND DP ≥ 5 AND different ALT base sequences (not just indel/SNV representation variants of the same event).
- **single**: Only one allele at the position (not a multi-variant position).

The classification SHALL check normalization_artifact first, then noise, then true_multi_allelic.

#### Scenario: Normalization artifact by length difference
- **WHEN** chr15:34362339 has AT>A (REF_len=2, ALT_len=1, VAF=0.066) and A>T (REF_len=1, ALT_len=1, VAF=0.0)
- **THEN** multiallelic_class="normalization_artifact" because REF/ALT lengths differ between alleles

#### Scenario: Noise classification
- **WHEN** a position has T>C (VAF=0.30, DP=50) and T>G (VAF=0.002, DP=2)
- **THEN** multiallelic_class="noise" because only one allele has sufficient signal

#### Scenario: True multi-allelic biology
- **WHEN** chr7:100958135 has T>C (VAF=0.438, DP=64) and T>G (VAF=0.111, DP=27)
- **THEN** multiallelic_class="true_multi_allelic" because both alleles have VAF≥0.01, DP≥5, and different ALT bases

#### Scenario: Missing depth data does not default to true_multi_allelic
- **WHEN** a position has two alleles where one allele has null DNA_ALT_DP_mean
- **THEN** the allele with null depth is treated as VAF=0, DP=0 for classification purposes
- **AND** the position is classified as noise (not true_multi_allelic) unless both alleles have non-null signal

### Requirement: allele_balance_ratio computed as max_alt_dp / total_alt_dp
The system SHALL compute `allele_balance_ratio` as `max(DNA_ALT_DP_mean across alleles) / sum(DNA_ALT_DP_mean across alleles)`. When pre-decomposition per-allele AF data is available, the system SHALL prefer `max(pre_norm_af) / sum(pre_norm_af)` for the ratio. The ratio SHALL be in range [0, 1], where 1.0 means one allele dominates completely and 0.5 means equal balance between two alleles. For positions with 3+ alleles, the denominator SHALL include all alleles (not just top-2).

#### Scenario: Two-allele balance ratio
- **WHEN** a site has two alleles with ALT_DP=15 and ALT_DP=5
- **THEN** allele_balance_ratio=15/(15+5)=0.75

#### Scenario: Three-allele balance ratio
- **WHEN** a site has three alleles with ALT_DP=40, ALT_DP=30, ALT_DP=20
- **THEN** allele_balance_ratio=40/(40+30+20)=0.444

#### Scenario: Equal-balance site
- **WHEN** a site has two alleles with ALT_DP=15 and ALT_DP=15
- **THEN** allele_balance_ratio=15/(15+15)=0.50

### Requirement: Pre-decomposition data is consumed by downstream metrics
The system SHALL add derived scalar columns from pre-decomposition parsing to `_CROSS_SAMPLE_COLS` so they survive the cross-sample scan: `gt_cooccurrence` (string encoding which alleles coexist in the same sample, computed from `pre_norm_*_gt_alleles`), `strand_balance` (per-allele F1R2/(F1R2+F2R1) ratio from `pre_norm_*_f1r2` and `pre_norm_*_f2r1`). When pre-decomposition data is unavailable, these columns SHALL be null-filled. The `vaf_sum` metric SHALL prefer pre-norm per-allele AF over `DNA_VAF_mean` when available.

#### Scenario: GT co-occurrence computed
- **WHEN** pre-norm Mutect2 VCF shows GT="0/1/2" at a multi-allelic position
- **THEN** gt_cooccurrence="0/1/2" is stored for all decomposed alleles at that position

#### Scenario: Strand balance computed
- **WHEN** pre-norm data shows F1R2_alt=10, F2R1_alt=30 for an allele
- **THEN** strand_balance=10/(10+30)=0.25 is stored for that allele

#### Scenario: Pre-norm data unavailable
- **WHEN** pre-decomposition VCFs do not exist for a sample
- **THEN** gt_cooccurrence=null, strand_balance=null, and vaf_sum uses DNA_VAF_mean fallback

### Requirement: OLD_MULTIALLELIC validates pre-norm parsing
The system SHALL use the `OLD_MULTIALLELIC` INFO field from the rescue VCF to validate pre-decomposition parsing. When `OLD_MULTIALLELIC` indicates a position was multi-allelic before decomposition, the system SHALL verify that `parse_pre_norm_multiallelic` found a corresponding multi-allelic record. Mismatches SHALL be logged as warnings but SHALL NOT halt the pipeline.

#### Scenario: OLD_MULTIALLELIC matches pre-norm
- **WHEN** OLD_MULTIALLELIC="T,C,G" at chr7:100958135 and pre-norm Mutect2 VCF has ALT="C,G" at that position
- **THEN** validation passes silently

#### Scenario: OLD_MULTIALLELIC mismatch
- **WHEN** OLD_MULTIALLELIC indicates multi-allelic but pre-norm VCF is missing
- **THEN** a warning is logged: "OLD_MULTIALLELIC indicates multi-allelic at chr7:100958135 but pre-norm VCF not found"
- **AND** processing continues with null-filled pre-norm columns
