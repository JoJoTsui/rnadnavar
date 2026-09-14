## Why

The filtered dataset (9.75M variants after NoConsensus exclusion) still contains 4.57M low-confidence variants (46.9%): 3.31M Reference-with-signal (VAF>0.05), 1.53M low-evidence tier variants (C5/C6), 189K category conflicts, 35K Germline with paradoxically low VAF, 2K Somatic with suspiciously high VAF, and 54K multi-allelic heterogeneity flags. The strongest single-caller artifact discriminator — strand bias (F1R2/F2R1) — is already in the parquet but drives ZERO flags or filters. Cross-sample recurrence is not computed at all. A unified 3-stage filtering strategy (hard exclusions, soft flags, confidence tiers) is needed to distinguish true variants from false positives.

## What Changes

### New Metrics to Compute

- **`n_recurrent_samples`**: Cross-sample recurrence count per (CHROM, POS, REF, ALT) — how many samples share this exact variant. Computed via `pl.len().over(['CHROM','POS','REF','ALT'])` on `combined_df`. Variants found in >20 samples are likely mapping artifacts or recurrent germline, not somatic.
- **`max_strand_bias_fisher_p`**: Fisher exact test p-value on F1R2/F2R1 ref/alt counts per caller. Computed from existing `*_SB` (Mutect2) and `pre_norm_*_f1r2/f2r1` columns. Values < 1e-3 with alt_reads ≥ 10 indicate extreme strand bias artifact.
- **`flag_germline_high_vaf`**: Germline with VAF > 0.85 AND GT=0/1 in 2+ callers — LOH or contamination, not germline.
- **`flag_somatic_loh`**: Somatic with VAF > 0.60 AND DNA_REF_DP_mean < 0.10 × DNA_DP_mean — LOH-driven artifact.
- **`flag_no_caller_support`**: N_SUPPORT_CALLERS == 0 — one-step replacement for manual C7D0 filters.
- **`flag_low_rna_mapq`**: RNA support from poorly-mapped reads (BAM_RT_mean_MQ < 2 AND RNA_DP_mean ≥ 20). Requires the MAPQ=255 fix from `fix-bam-stats-rust`.

### 3-Stage Filtering Strategy

- **Stage 0 — Hard Exclusions (auto-drop ~1.45M, ~15% of filtered)**: Variants that fail physical, biological, or evidence constraints. Includes: N_SUPPORT_CALLERS==0, flag_vaf_overflow, multiallelic_class=="noise", insufficient depth/alt across both modalities, paradoxical FILTER+VAF combinations (Germline VAF<0.10, Somatic VAF>0.60 with LOH, Reference VAF>0.10 with ALT_DP≥3).
- **Stage 1 — Soft Flags (mark for review)**: Variants that pass hard filters but have suspicious characteristics. Includes: flag_category_conflict, flag_multi_allelic_heterogeneity, flag_rna_rescued, low_confidence modality, high cross-sample recurrence, single-caller-only Somatic/Germline, strand bias, low RNA MAPQ.
- **Stage 2 — Confidence Tiers**: HIGH (no soft flags + C1/C2/C3/C4 with database + minimum caller support), MEDIUM (no soft flags + C2/C3/C4 without database), LOW (passes hard filters but has any soft flag).

### Both Filtered and Unfiltered Statistics

- When running in filtered mode, the system SHALL emit both filtered and unfiltered versions of major statistics: `sample_summary.tsv` + `sample_summary_filtered.tsv`, `set_summary.tsv` + `set_summary_filtered.tsv`, `dataset_summary.tsv` + `dataset_summary_unfiltered.tsv`, `rescue_validation_summary.tsv` (both).

### Caller Overlap Output

- Add a real caller overlap statistic: `caller_overlap_matrix.tsv` (pairwise co-occurrence matrix: N variants where caller_i AND caller_j both call) and `caller_support_distribution.tsv` (histogram: how many variants called by exactly 1, 2, ..., 6 callers). Keep the existing tier distribution as `tier_distribution.tsv` (renamed from the misnamed `caller_overlap.tsv`).

## Capabilities

### New Capabilities

- `low-confidence-variant-filtering`: 3-stage variant filtering system with hard exclusions, soft flags, and confidence tiering. Computes cross-sample recurrence, strand bias Fisher test, and paradoxical VAF flags. Emits both filtered and unfiltered statistics.

### Modified Capabilities

- `variant-tiering-stats`: Add `caller_overlap_matrix.tsv` and `caller_support_distribution.tsv` outputs. Rename `caller_overlap.tsv` to `tier_distribution.tsv`. Emit both filtered and unfiltered sample/set/dataset summaries.
- `biological-evidence-flags`: Add `flag_germline_high_vaf`, `flag_somatic_loh`, `flag_no_caller_support`, `flag_low_rna_mapq` flags. Wire strand bias and BAM per-position metrics into flag computation.

## Impact

- `bin/vcf_stats/seq2neo/statistics.py` — ~120 lines: add `n_recurrent_samples` computation (~15 lines), add `max_strand_bias_fisher_p` computation (~25 lines), add 4 new biological flags (~40 lines), add `compute_confidence_tier` function (~20 lines), rename `caller_overlap_distribution` + add `caller_overlap_matrix` + `caller_support_distribution` (~20 lines)
- `bin/vcf_stats/seq2neo/cli.py` — ~80 lines: add `build_low_confidence_filter` function (~30 lines), emit both filtered+unfiltered summaries (~25 lines), wire new outputs (~15 lines), add `--confidence-tier` CLI flag (~10 lines)
- `bin/vcf_stats/seq2neo/visualizer.py` — ~20 lines: add confidence tier color registry, add recurrence plot
- Parquet schema: +6 new columns (n_recurrent_samples, max_strand_bias_fisher_p, flag_germline_high_vaf, flag_somatic_loh, flag_no_caller_support, flag_low_rna_mapq, confidence_tier)
- New output files: `caller_overlap_matrix.tsv`, `caller_support_distribution.tsv`, `tier_distribution.tsv`, `sample_summary_filtered.tsv`, `set_summary_filtered.tsv`, `dataset_summary_unfiltered.tsv`, `confidence_tier_summary.tsv`
- **Depends on**: `fix-bam-stats-rust` (MAPQ=255 fix needed for `flag_low_rna_mapq`), `fix-multiallelic-logic` (corrected classification needed for Stage 0 noise exclusion, corrected flags needed for Stage 1)
