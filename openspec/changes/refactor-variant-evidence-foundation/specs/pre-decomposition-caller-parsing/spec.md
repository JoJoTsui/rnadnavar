# pre-decomposition-caller-parsing

## Purpose

Parse pre-normalization (pre-`vt decompose`) Mutect2 and DeepSomatic caller VCFs from `variant_calling/` directories to recover full multi-allelic allele depth (AD) arrays, per-allele allele fraction (AF), per-allele strand bias (F1R2/F2R1), and genotype (GT) co-occurrence information that is discarded when vt decompose reduces FORMAT to GT:DP only.

## ADDED Requirements

### Requirement: Pre-decomposition caller VCF paths are discoverable
The system SHALL construct pre-decomposition caller VCF paths from the sample output directory by replacing `normalized/{caller}/.../*.dec.norm.vcf.gz` with `variant_calling/{caller}/.../*.vcf.gz` for Mutect2 and DeepSomatic callers. For Strelka callers, the system SHALL skip pre-decomposition parsing because Strelka output is always biallelic. Path discovery SHALL use the existing `base_output_dir` and `dir_name` fields from the sample manifest.

#### Scenario: Mutect2 pre-norm VCF exists
- **WHEN** a sample has `variant_calling/mutect2/{prefix}DT_vs_{prefix}DN/{prefix}DT_vs_{prefix}DN.mutect2.vcf.gz`
- **THEN** the pre-decomposition parser extracts multi-allelic data from that file
- **AND** the file is NOT required — if absent, parsing falls back to normalized-only mode with a warning

#### Scenario: Strelka caller has no pre-norm multi-allelic data
- **WHEN** a Strelka caller is processed
- **THEN** pre-decomposition parsing is skipped because Strelka VCFs contain zero multi-allelic records

#### Scenario: RNA caller pre-norm VCF path
- **WHEN** an RNA caller's pre-decomposition VCF is at `variant_calling/mutect2/{prefix}RT_vs_{prefix}DN/{prefix}RT_vs_{prefix}DN.mutect2.vcf.gz`
- **THEN** the system discovers it using the same path transformation as DNA callers

### Requirement: Full AD arrays are extracted from multi-allelic records
The system SHALL parse FORMAT/AD from pre-decomposition VCF records where ALT contains a comma, extracting the full allele depth array `AD=[ref, alt1, alt2, ...]` for Mutect2 and `AD=[ref, alt1, alt2, ...]` for DeepSomatic. The system SHALL store per-allele AD values keyed by `(CHROM, POS, REF, ALT_allele)` in an allele registry.

#### Scenario: Mutect2 multi-allelic record with 2 ALT alleles
- **WHEN** a pre-norm Mutect2 record has REF=CA, ALT=C,CAA and FORMAT/AD=30,5,5
- **THEN** AD_REF=30, AD_ALT_C=5, AD_ALT_CAA=5 are stored in the registry

#### Scenario: Biallelic record (single ALT, no comma)
- **WHEN** a pre-norm record has REF=A, ALT=G (single allele, no comma in ALT)
- **THEN** the record is parsed as-is with AD=[ref, alt] without multi-allelic registry entry

### Requirement: Per-allele AF and strand bias are extracted
For Mutect2 pre-decomposition records, the system SHALL extract FORMAT/AF (per-allele allele fraction), FORMAT/F1R2 (forward-strand read counts per allele), and FORMAT/F2R1 (reverse-strand read counts per allele). For DeepSomatic, the system SHALL extract FORMAT/VAF. The system SHALL compute per-allele strand balance as `F1R2_alt / (F1R2_alt + F2R1_alt)` for each alternate allele.

#### Scenario: Per-allele strand bias extraction
- **WHEN** a Mutect2 multi-allelic record has F1R2=9,1,2 and F2R1=20,4,3
- **THEN** for alt allele 1: strand_balance = 1/(1+4) = 0.20
- **AND** for alt allele 2: strand_balance = 2/(2+3) = 0.40

#### Scenario: DeepSomatic per-allele VAF
- **WHEN** a DeepSomatic multi-allelic record has VAF=0.095,0.053
- **THEN** VAF is stored per allele: alt1_VAF=0.095, alt2_VAF=0.053

### Requirement: GT co-occurrence is captured
The system SHALL parse FORMAT/GT from pre-decomposition multi-allelic records. When GT contains multiple non-zero allele indices (e.g., GT=0/1/2), the system SHALL record which alleles co-occur in the same sample, enabling tumor heterogeneity analysis.

#### Scenario: Co-occurring alleles detected
- **WHEN** a Mutect2 multi-allelic record has GT=0/1/2 in the tumor sample
- **THEN** alleles 1 and 2 are flagged as co-occurring in the same sample
- **AND** the `gt_cooccurrence` field for both decomposed alleles records the other allele's ALT string

### Requirement: Allele registry joins with normalized parse
The system SHALL build an allele registry as a polars DataFrame from pre-decomposition data, keyed on `(CHROM, POS, REF, ALT)`, and left-join it with the existing normalized parse during `process_single_sample()`. For positions that exist in the normalized data but not in the pre-decomposition registry, multi-allelic columns SHALL be null-filled.

#### Scenario: Join enriches decomposed records
- **WHEN** two decomposed normalized records at chr1:1318647 (CA>C and C>CA) are joined with pre-norm registry entry (REF=CA, ALT=C,CAA, AD=30,5,5)
- **THEN** the CA>C record gets allele_rank=1, n_alleles_at_site=2, original_alt_alleles=["C","CAA"]
- **AND** the C>CA record similarly gets allele_rank=2 and the same n_alleles_at_site and original_alt_alleles
