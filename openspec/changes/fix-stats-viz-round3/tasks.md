## 1. Chromosome Ordering Fix (P0)

- [x] 1.1 Create `_chrom_sort_list(df, col="CHROM")` helper that returns sorted chromosome list from `_sort_chromosomes` output
- [x] 1.2 Fix `plot_vc_distribution`: pass `sort=chrom_order` to `alt.X()` when `group_col=="CHROM"`
- [x] 1.3 Fix `plot_variant_type_distribution`: same pattern
- [x] 1.4 Fix `plot_ti_tv_ratio`: same pattern
- [x] 1.5 Fix `plot_cross_modality`: same pattern
- [x] 1.6 Fix `plot_filter_distribution`: same pattern
- [x] 1.7 Fix `plot_redi_evidence`: same pattern
- [x] 1.8 Fix `plot_cosmic_gnomad_annotation`: same pattern (uses "group" column rename)
- [x] 1.9 Fix chromosome summary TSV: apply `_sort_chromosomes()` before `write_tsv()` in cli.py chromosome wise section

## 2. DP Threshold Sweep Classification (P1)

- [x] 2.1 Add per-FILTER classification inner loop to `compute_dp_threshold_sweep()` matching VAF sweep pattern
- [x] 2.2 Skip per-classification breakdown for BAM pileup DP columns (they're not caller-specific)
- [ ] 2.3 Verify DP sweep TSV now has `classification` column (verification — needs pipeline run)

## 3. Per-Sample Charts Grid by Set (P1)

- [x] 3.1 Update `plot_per_sample_distribution` (chart 08): add `set_number` faceting via `alt.Row("set_number:N")`
- [x] 3.2 Add `set_number` to `sample_tier_summary()` group-by in statistics.py
- [x] 3.3 Update `plot_per_sample_tier_distribution` (chart 22): add `set_number` faceting

## 4. BAM Coverage Violin Fix (P1)

- [x] 4.1 Clip depth values to [0, 2000] before `transform_density` in `plot_bam_coverage_violin`
- [x] 4.2 Add subtitle annotation showing count of values > 2000

## 5. DP Scale and Coloring (P1)

- [x] 5.1 Add `scale=alt.Scale(type="symlog")` to `plot_bam_dp_distribution` y-axis
- [x] 5.2 Fix violin/boxplot color mismatch in `plot_vaf_distribution`: add `color` to boxplot `enc`
- [x] 5.3 Fix violin/boxplot color mismatch in `plot_dp_distribution`: add `color` to boxplot `enc`
- [x] 5.4 Fix violin/boxplot color mismatch in `plot_bam_dp_distribution`: add `color` to boxplot `enc`

## 6. Depth Capping at 2000 (P1)

- [x] 6.1 Cap `plot_dp_boxplot_per_tier`: add `scale=alt.Scale(domain=[0, 2000])`
- [x] 6.2 Cap `plot_dp_distribution`: add domain cap
- [x] 6.3 Cap `plot_dna_vs_rna_dp`: add domain cap on both axes
- [x] 6.4 Cap `plot_ref_alt_dp_scatter`: add domain cap on both axes
- [x] 6.5 Cap `plot_per_tier_dp_boxplot`: add domain cap
- [x] 6.6 Add outlier count annotation to each capped chart (subtitle: "N values > 2000 excluded")

## 7. ML Threshold Guidance Statistics (P2)

- [x] 7.1 Add `partition` computed column in cli.py: chr1→test, chr21+22→val, rest→train
- [x] 7.2 Add `compute_filter_vaf_dp_cross_tab(df, partition_col)` to statistics.py
- [x] 7.3 Output `cross_tab_filter_vaf_dp.tsv` with FILTER × VAF_bin × DP_bin × partition × count
- [x] 7.4 Add `per_chrom_partition_summary.tsv`: per-partition n_variants, n_somatic, disease distribution, mean VAF/DP
- [x] 7.5 Add `disease_x_chrom_partition.tsv`: disease × partition cross-tab (validates zero-shot isolation)
- [x] 7.6 Add `low_vaf_rna_support.tsv`: variants with VAF < 0.05 where N_RNA_CALLERS_SUPPORT >= 2
- [x] 7.7 Add heatmap visualization: FILTER × VAF_bin per partition (chart 43)
- [x] 7.8 Add stacked bar: disease × partition variant counts (chart 44)

## 8. FP Cross-Tabulation (P2)

- [x] 8.1 Add `compute_fp_cross_tab(df)` to statistics.py: non-Somatic variants × N_SUPPORT_CALLERS × VAF_bin × DP_bin
- [x] 8.2 Output `fp_cross_tab.tsv`
- [x] 8.3 Add heatmap: FILTER × N_SUPPORT_CALLERS counts (chart 45)

## 9. Somatic Modality Sub-Classification (P2)

- [x] 9.1 Add `compute_somatic_modality(df)` to statistics.py: derive from `caller_tier` — C1→Multi, C2/C5→DNA_only, C3/C6→RNA_only, C4/C7→Weak
- [x] 9.2 Add `somatic_modality_summary.tsv`: modality × count, mean_vaf, mean_dp, n_cosmic, n_gnomad
- [x] 9.3 Add `somatic_modality_x_disease.tsv`: modality × disease cross-tab
- [x] 9.4 Add boxplot: VAF distribution per somatic_modality (chart 46)
- [x] 9.5 Add stacked bar: somatic_modality × disease (chart 47)
- [x] 9.6 Add `N_DNA_CALLERS_SUPPORT` and `N_RNA_CALLERS_SUPPORT` to `_CROSS_SAMPLE_COLS`

## 10. Scientific Publishing Theme (P2)

- [x] 10.1 Add `_register_publishing_theme()` to visualizer.py: white bg, 11pt labels, 13pt titles, Arial, minimal grid
- [x] 10.2 Add `--theme` CLI flag with choices `["default", "publishing"]`
- [x] 10.3 Call `alt.themes.enable("publishing")` before chart generation when `--theme publishing`

## 11. Chart ID Collision Fix

- [x] 11.1 Renumber chart IDs to avoid 35_, 36_, 37_ collisions — renumbered to 40_, 41_, 42_; new charts use 43-47

## 12. Documentation

- [x] 12.1 Document somatic modality sub-classification derivation from C-tier (docs/somatic_modality.md)
- [x] 12.2 Document ML partition strategy (chr1=test, chr21-22=val, rest=train) (docs/ml_partition_strategy.md)
- [x] 12.3 Document FP cross-tabulation methodology (docs/fp_cross_tabulation.md)

## 13. Verification

- [ ] 13.1 Run pipeline on 8 samples, verify chromosome ordering in all chromosome charts
- [ ] 13.2 Verify DP sweep TSV has classification column
- [ ] 13.3 Verify per-sample charts are gridded by set
- [ ] 13.4 Verify coverage violin is readable (not blank/white)
- [ ] 13.5 Verify DP boxplots are capped at 2000
- [ ] 13.6 Verify somatic modality summary TSV is generated
- [ ] 13.7 Verify ML cross-tab TSVs are generated
- [ ] 13.8 Run test suite, confirm no regressions
