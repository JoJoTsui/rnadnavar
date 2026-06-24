## ADDED Requirements

### Requirement: Germline high VAF flag
The system SHALL set `flag_germline_high_vaf=True` for variants where FILTER="Germline" AND DNA_VAF_mean > 0.85 AND GT is 0/1 or 0/2 in at least 2 GT-bearing callers. A heterozygous germline variant at VAF > 0.85 is inconsistent with constitutional heterozygosity and may indicate LOH or sample contamination.

#### Scenario: Germline with near-homozygous VAF
- **WHEN** a variant has FILTER=Germline, DNA_VAF_mean=0.92, GT=0/1 in 3 callers
- **THEN** flag_germline_high_vaf=True

### Requirement: Somatic LOH flag
The system SHALL set `flag_somatic_loh=True` for variants where FILTER="Somatic" AND DNA_VAF_mean > 0.60 AND DNA_REF_DP_mean < 0.10 × DNA_DP_mean. This pattern indicates loss of heterozygosity where the reference allele is nearly absent, which is a technical artifact of LOH regions rather than a true somatic mutation.

#### Scenario: Somatic with LOH pattern
- **WHEN** a variant has FILTER=Somatic, DNA_VAF_mean=0.85, DNA_REF_DP_mean=5, DNA_DP_mean=100
- **THEN** flag_somatic_loh=True (REF_DP is 5% of DP, VAF is 85%)

### Requirement: No caller support flag
The system SHALL set `flag_no_caller_support=True` for variants where N_SUPPORT_CALLERS == 0. This replaces manual C7D0 tier filters with a single boolean flag.

#### Scenario: Variant with no caller support
- **WHEN** a variant has N_SUPPORT_CALLERS=0
- **THEN** flag_no_caller_support=True

### Requirement: Low RNA MAPQ flag
The system SHALL set `flag_low_rna_mapq=True` for variants where BAM_RT_mean_MQ < 2 AND RNA_DP_mean ≥ 20. This indicates RNA support from poorly-mapped reads. The flag requires the MAPQ=255 fix (MAPQ=255 excluded from mean_mapq) to produce meaningful values.

#### Scenario: RNA support from low-MAPQ reads
- **WHEN** a variant has BAM_RT_mean_MQ=1.5, RNA_DP_mean=50
- **THEN** flag_low_rna_mapq=True

#### Scenario: RNA support from well-mapped reads
- **WHEN** a variant has BAM_RT_mean_MQ=30, RNA_DP_mean=50
- **THEN** flag_low_rna_mapq=False

### Requirement: Strand bias and BAM per-position metrics wired into flags
The system SHALL use existing but currently unused strand bias data (`*_SB`, `pre_norm_*_f1r2/f2r1`, `BAM_{DT,RT}_F1R2*/F2R1*`) and BAM per-position metrics (`BAM_*_mean_BQ`, `BAM_*_mean_MQ`) in biological flag computation. The `max_strand_bias_fisher_p` metric SHALL be computed from these fields and used in Stage 1 soft flagging.

#### Scenario: Strand bias data used for flagging
- **WHEN** a variant has Mutect2 F1R2_alt=2, F2R1_alt=48, alt_reads=50
- **THEN** max_strand_bias_fisher_p is computed and < 1e-3
- **AND** the variant receives a "strand_bias" soft flag
