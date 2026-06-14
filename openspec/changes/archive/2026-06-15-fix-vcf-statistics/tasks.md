## 1. Bug Fixes — Dashboard and Color Legend

- [x] 1.1 Fix dashboard.html rendering: extract `<body>` inner content from each chart's `to_html()` output, assemble into single valid HTML page with one `<html>`/`<head>`/`<body>` and shared vega-embed CDN loaded once
- [x] 1.2 Fix VC distribution chart color legend: verify altair color scale is properly encoding VC field and legend is visible; test against real data to confirm Somatic/Germline/Reference/Artifact labels appear
- [x] 1.3 Verify all 12+ charts render in the fixed dashboard.html by running e2e test against 4 real samples

## 2. Tiering Integration

- [x] 2.1 Create `bin/vcf_stats/seq2neo/tiering_stats.py` bridge module that imports `TieringEngine` from `bin/vcf_stats/tiering_engine.py` and computes CxDy tiers from parsed FILTERS_NORMALIZED / FILTER_NORMALIZED_* fields already present in rescue VCF INFO data
- [x] 2.2 Add `"FILTERS_NORMALIZED"` and `"FILTER_NORMALIZED_DNA_TUMOR"` pattern fields to `RESCUE_STRING_FIELDS` in `rescue_parser.py` to capture tiering-relevant INFO fields
- [x] 2.3 Compute tier columns (`caller_tier`, `database_tier`, `final_tier`, `tier_quality`) in `process_single_sample()` of `cli.py` after joining caller columns, store in `variant_details.parquet`
- [x] 2.4 Add `tier_summary()` function to `tiering_stats.py`: per-tier aggregate statistics (count, mean VAF, mean DP, mean REF_DP, mean ALT_DP, variant type distribution, Ti/Tv)
- [x] 2.5 Write tests in `test_seq2neo_stats.py` for tier computation and per-tier statistics

## 3. REF_DP and ALT_DP Statistics

- [x] 3.1 Add `{caller}_REF_DP` and `{caller}_ALT_DP` columns in `caller_parser.py` `join_caller_columns()` — rename existing AD_REF→REF_DP and AD_ALT→ALT_DP for clarity, keep backward-compatible aliases
- [x] 3.2 Add `DNA_REF_DP_mean`, `DNA_ALT_DP_mean`, `RNA_REF_DP_mean`, `RNA_ALT_DP_mean` to `statistics.py` `compute_mean_columns()`
- [x] 3.3 Add `mean_dna_ref_dp`, `mean_dna_alt_dp`, `mean_rna_ref_dp`, `mean_rna_alt_dp` to `sample_summary()` in `statistics.py`
- [x] 3.4 Add REF_DP and ALT_DP scatter plots to `visualizer.py` (DNA vs RNA mean REF_DP, DNA vs RNA mean ALT_DP, colored by VC or tier)
- [x] 3.5 Write tests for REF_DP/ALT_DP computation and visualization

## 4. Multi-Level Aggregation

- [x] 4.1 Add `dataset_summary()` function to `statistics.py`: whole-dataset aggregates across all samples (total variants, variant type distribution, Ti/Tv, VC distribution, caller support histogram, tier distribution)
- [x] 4.2 Generate `dataset_summary.csv` in `cli.py` after all samples are processed
- [x] 4.3 Generate `tier_summary.csv` from per-tier aggregates across all samples
- [x] 4.4 Write tests for dataset-level and tier-level aggregation outputs

## 5. Visualization Updates — Violin Plots, Per-Tier Charts

- [x] 5.1 Replace VAF boxplot with violin-style plot (layered density area charts) in `visualizer.py` — `plot_vaf_distribution_violin()` function per caller, faceted by tier
- [x] 5.2 Add `plot_tiered_dp_violin()` — DP distribution violin plots per caller, faceted by tier
- [x] 5.3 Add `plot_tiered_gt_concordance()` — GT concordance faceted by tier
- [x] 5.4 Add `plot_tiered_caller_overlap()` — N_SUPPORT_CALLERS histogram faceted by tier
- [x] 5.5 Add `plot_tiered_variant_types()` — SNV/INS/DEL/MNV distribution faceted by tier
- [x] 5.6 Update `generate_dashboard()` to include all new charts
- [x] 5.7 Write tests for all new chart functions

## 6. Per-Sample Visualization Fixes

- [x] 6.1 Rewrite `plot_per_sample_violin()` → `plot_per_sample_distribution()`: horizontal bar chart with sample_id on y-axis, total_variants on x-axis, sorted descending
- [x] 6.2 Add percentage annotations to `plot_cross_modality()` showing proportion relative to total variants per set
- [x] 6.3 Fix `plot_per_sample_violin()` → `plot_per_sample_distribution()` to use sample ID labels on y-axis (not disease)
- [x] 6.4 Write tests for updated per-sample and cross-modality charts

## 7. BAM Statistics

- [x] 7.1 Create `bin/vcf_stats/seq2neo/bam_stats.py` module using pysam for BAM file parsing
- [x] 7.2 Implement per-sample per-modality BAM statistics: locate BAM files in sample output directories, compute total reads, mapped reads, mapping rate, mean coverage, mean insert size
- [x] 7.3 Generate `bam_stats.csv` in CLI pipeline with one row per sample per modality
- [x] 7.4 Extract per-variant caller FORMAT fields (SB, FAD from Mutect2; AU/CU/GU/TU from Strelka) and append to `variant_details.parquet`
- [x] 7.5 Write tests for BAM statistics module (use a small test BAM or skip if BAM unavailable)

## 8. Integration and Final Testing

- [x] 8.1 Wire all new modules into `cli.py`: tiering_stats, bam_stats, new charts, updated output files
- [x] 8.2 Run full e2e test on 4 real samples (1 per set) — verify all new outputs generated correctly
- [x] 8.3 Run full test suite and confirm all tests pass
- [x] 8.4 Verify dashboard.html renders all charts correctly by checking HTML structure
- [x] 8.5 Run on all 65 complete samples to validate production readiness
