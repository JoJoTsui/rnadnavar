# multi-allelic-analysis

## Purpose

Detect multi-variant positions by (CHROM, POS) grouping, classify them into normalization artifacts, noise, or true multi-allelic biology, compute per-site metrics (n_alleles, VAF sum, allele balance ratio, category conflict), and flag tumor heterogeneity signals for preservation rather than filtering.

## ADDED Requirements

### Requirement: Multi-variant positions are detected
The system SHALL detect positions where more than one variant record exists at the same (CHROM, POS). Detection SHALL use polars group_by on the combined or per-sample DataFrame. For each such position, the system SHALL aggregate allele-level data: FILTER categories, ALT alleles, DNA_VAF_mean, RNA_VAF_mean, DNA_DP_mean, per-caller support counts, and pre-decomposition allele registry data when available.

#### Scenario: Position with two alleles detected
- **WHEN** chr7:100958135 has T>C (NoConsensus, VAF=0.438) and T>G (NoConsensus, VAF=0.111)
- **THEN** n_alleles_at_site=2, vaf_sum=0.549, category_conflict=False (both NoConsensus)

#### Scenario: Position with single allele
- **WHEN** a (CHROM, POS) has exactly one variant record
- **THEN** n_alleles_at_site=1, no multi-allelic metrics are computed (null-filled)

### Requirement: Multi-variant positions are classified into three categories
The system SHALL classify each multi-variant position into one of three categories based on allele characteristics:

- **Normalization artifact**: Alleles have different REF/ALT string lengths (e.g., AT>A deletion vs A>T SNV at same position) OR one allele has VAF≈0 and DP≈0. These represent the same biological event encoded differently after vt decompose.
- **One real allele + noise**: One allele has non-zero VAF and DP, all other alleles have VAF≈0 or DP≈0. The signal allele is the biological variant; others are technical noise.
- **True multi-allelic biology**: Two or more alleles have non-zero VAF (≥0.01) AND non-zero DP (≥5) AND different ALT base sequences (not just indel/SNV representation variants of the same event).

#### Scenario: Normalization artifact classification
- **WHEN** chr15:34362339 has AT>A (VAF=0.066) and A>T (VAF=0.0)
- **THEN** the position is classified as normalization artifact because one allele has zero VAF and the REF/ALT lengths differ (2>1 vs 1>1)

#### Scenario: True multi-allelic biology classification
- **WHEN** chr7:100958135 has T>C (VAF=0.438, DP=64) and T>G (VAF=0.111, DP=27)
- **THEN** the position is classified as true multi-allelic biology because both alleles have non-zero VAF, non-zero DP, and different ALT bases

### Requirement: Per-site metrics are computed
For each multi-variant position, the system SHALL compute: `n_alleles_at_site` (count of variant records at this position), `vaf_sum` (sum of DNA_VAF_mean across all alleles, replaces per-allele AF from pre-norm when available), `allele_balance_ratio` (max_alt_dp / total_alt_dp, reflecting dominance of the major allele), `category_conflict` (True if unique FILTER categories among alleles > 1), and `multiallelic_class` (one of: "normalization_artifact", "noise", "true_multi_allelic").

#### Scenario: Allele balance ratio computation
- **WHEN** a site has two alleles with ALT_DP=15 and ALT_DP=5
- **THEN** total_alt_dp=20, allele_balance_ratio=15/20=0.75

#### Scenario: Category conflict detection
- **WHEN** chr12:95287821 has G>A (FILTER=Germline) and G>T (FILTER=NoConsensus)
- **THEN** category_conflict=True (Germline vs NoConsensus)

### Requirement: True multi-allelic sites are preserved during filtering
The system SHALL NOT filter out variants at true multi-allelic sites solely because of category conflict or low minor-allele VAF. These sites represent potential tumor heterogeneity signals. Instead, the system SHALL set `flag_multi_allelic_heterogeneity=True` for all alleles at true multi-allelic sites.

#### Scenario: Multi-allelic site survives filtering
- **WHEN** a unified filter is configured to drop NoConsensus variants
- **AND** chr7:100958135 has two NoConsensus alleles (T>C at VAF=0.438, T>G at VAF=0.111) classified as true multi-allelic
- **THEN** both alleles SHALL survive filtering despite being NoConsensus
- **AND** flag_multi_allelic_heterogeneity=True is set for both

### Requirement: Strand bias analysis at multi-allelic sites
When pre-decomposition strand bias data (F1R2/F2R1) is available from the allele registry, the system SHALL compute per-allele strand balance and flag alleles where `strand_balance < 0.10 OR strand_balance > 0.90` as potentially artifactual. For multi-allelic positions where one allele has extreme strand bias and the other does not, the biased allele SHALL be flagged for review but NOT automatically filtered.

#### Scenario: Strand-biased allele at true multi-allelic site
- **WHEN** a multi-allelic site has allele A (strand_balance=0.45, balanced) and allele B (strand_balance=0.02, extreme bias)
- **THEN** allele B gets flag_extreme_strand_bias=True
- **AND** both alleles are preserved (not auto-filtered) because the site is true multi-allelic
