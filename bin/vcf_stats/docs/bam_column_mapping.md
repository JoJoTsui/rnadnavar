# BAM Pileup Column to Caller Field Mapping

## Column Naming Convention

BAM pileup columns use the naming pattern: `BAM_{bam_type}_{metric}`

Where:
- `bam_type` = `DN` (DNA Normal), `DT` (DNA Tumor), or `RT` (RNA Tumor)
- `metric` = one of the pileup output columns

### Pileup Output Columns

| Column | Description |
|--------|-------------|
| `DP` | Total depth at position |
| `REF_DP` | Reference allele depth |
| `ALT_DP` | Alternate allele depth |
| `F1R2_ref` | Forward read 1, reverse read 2 — ref allele |
| `F2R1_ref` | Forward read 2, reverse read 1 — ref allele |
| `F1R2_alt` | Forward read 1, reverse read 2 — alt allele |
| `F2R1_alt` | Forward read 2, reverse read 1 — alt allele |
| `mean_BQ` | Mean base quality at position |
| `mean_MQ` | Mean mapping quality at position |

### Renamed Columns in Parquet

After pileup, columns are renamed for joining into per-sample parquet:

```
DP       → BAM_DT_DP, BAM_RT_DP
REF_DP   → BAM_DT_REF_DP, BAM_RT_REF_DP
ALT_DP   → BAM_DT_ALT_DP, BAM_RT_ALT_DP
F1R2_alt → BAM_DT_F1R2_alt, BAM_RT_F1R2_alt
...
```

## Validation Mapping

In `bam_validation.py`, BAM pileup columns are compared against caller VCF fields:

| BAM Column | Caller Field | Comparison |
|-----------|-------------|------------|
| `BAM_DT_DP` | `DNA_mutect2_DP`, `DNA_deepsomatic_DP`, `DNA_strelka_DP` | DP correlation, MAD |
| `BAM_RT_DP` | `RNA_mutect2_DP`, `RNA_deepsomatic_DP`, `RNA_strelka_DP` | DP correlation, MAD |
| `BAM_DT_ALT_DP / BAM_DT_DP` | `DNA_mutect2_VAF`, `DNA_deepsomatic_VAF` | VAF correlation |
| `BAM_RT_ALT_DP / BAM_RT_DP` | `RNA_mutect2_VAF`, `RNA_deepsomatic_VAF` | VAF correlation |
| `BAM_DT_F1R2_alt` | `DNA_mutect2_SB` | Strand bias check |
| `BAM_RT_F1R2_alt` | `RNA_mutect2_SB` | Strand bias check |

## Strelka VAF Note

Strelka VAF can exceed 1.0 — this is **correct by design**. Strelka uses
tier-1 filtered depth (`DP_tier1`) as its VAF denominator, which excludes
low-quality reads. The alternate allele count can exceed this filtered
denominator. Visualization clamps Strelka VAF to [0, 1] for display but
the raw value is preserved in the data.
