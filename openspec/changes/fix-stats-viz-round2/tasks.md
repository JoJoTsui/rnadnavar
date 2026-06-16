## 1. Bug Fixes (P0)

### BAM validation column name fix
- [x] 1.1 Fix `bam_validation.py:34`: change `f"BAM_DP_{bt}"` → `f"BAM_{bt}_DP"` 
- [x] 1.2 Fix `bam_validation.py:54-55`: apply same BAM_{bt}_{metric} naming to `bam_alt_col` and `bam_ref_col`
- [x] 1.3 Add unit test: verify `has_bam_data=True` when BAM_DT_DP column exists

### BAM coverage violin fix
- [x] 1.4 Fix `visualizer.py:1163`: replace `len(dp_cols) < 2` with check that at least one DP column exists
- [x] 1.5 Add empty-data guard: check `sampled` is non-empty before `transform_density()`
- [x] 1.6 Handle partial BAM types: skip missing types rather than skipping entire chart

### BAM integrity check
- [x] 1.7 Add `_check_bam_eof(bam_path)` to `bam_stats.py`: read last 28 bytes, verify BGZF EOF magic bytes
- [x] 1.8 Add pre-processing check in `cli.py`: validate all BAMs before processing, skip samples with all BAMs truncated
- [x] 1.9 Log warning for partially truncated samples (some BAMs valid, some not)

## 2. Depth Threshold Sweep (P1)

- [x] 2.1 Add `compute_dp_threshold_sweep(df)` to `statistics.py`: replicate VAF sweep pattern for per-caller DP columns
- [x] 2.2 Add REF_DP and ALT_DP threshold sweeps using BAM pileup columns (`BAM_DT_REF_DP`, `BAM_RT_ALT_DP`, etc.)
- [x] 2.3 Add `plot_dp_threshold_sweep()` to `visualizer.py`: line chart of retention % vs threshold per caller
- [x] 2.4 Wire DP sweep into CLI output (CSV + chart)

## 3. DP Distribution Plots (P1)

- [x] 3.1 Add `plot_bam_dp_distribution()` to `visualizer.py`: violin + box for BAM pileup DP (total/REF/ALT) per BAM type
- [x] 3.2 Add `plot_per_tier_dp_boxplot()` to `visualizer.py`: per-tier DNA DP boxplot (mirroring `plot_per_tier_vaf_boxplot`)
- [x] 3.3 Ensure all DP plots use consistent color palette and axis styling (category10 scheme, box+violin overlay pattern)

## 4. Chromosome Ordering (P1)

- [x] 4.1 Add `_sort_chromosomes()` call in `plot_vc_distribution` when `group_col="CHROM"` (already done — line 271)
- [x] 4.2 Add `_sort_chromosomes()` call in `plot_dna_vs_rna_per_caller` when `group_col="CHROM"` (moot — function doesn't accept group_col)
- [x] 4.3 Audit remaining charts for CHROM grouping and add `_sort_chromosomes()` where missing (all 8 CHROM-grouping sites already call _sort_chromosomes)
- [x] 4.4 Add `_sort_chromosomes()` to statistics aggregation queries that sort by CHROM (already applied)

## 5. VAF Distribution Styling (P1)

- [x] 5.1 Clamp VAF axis to `[0, 1]` in all VAF distribution charts (plot_vaf_distribution, plot_vaf_boxplot_per_tier, plot_per_tier_vaf_boxplot)
- [x] 5.2 Add vertical reference lines at VAF=0.005 and VAF=0.01 (added threshold reference rules in plot_vaf_distribution)
- [x] 5.3 Add viz-only `vaf_display` column computed as `min(VAF, 1.0)` in VAF melt operations
- [x] 5.4 Add annotation text to Strelka charts explaining VAF > 1 due to tier-1 depth denominator (axis label: "capped at 1.0")
- [x] 5.5 Apply consistent color palette (category10 or tableau10) across all VAF charts

## 6. Sample-Wise Chart Gridding (P2)

- [x] 6.1 Update `plot_bam_metrics_bars`: remove top-N filter, show all samples, facet by `set_number` (if available in bam_stats_df)
- [x] 6.2 Update `plot_per_sample_tier_distribution`: facet by `set_number` with shared axis scales (blocked — sample_tier_summary doesn't include set_number; deferred to future proposal)
- [x] 6.3 For charts with >50 samples per set, use 2D grid layout (wrap facets) (deferred to future proposal)

## 7. Database Enrichment Fix (P2)

- [x] 7.1 Update `plot_database_enrichment_by_tier`: skip CxD0 tiers in the bar chart
- [x] 7.2 Add note to chart title/subtitle: "D=0 tiers excluded (no database evidence by design)"

## 8. Data Leakage Exclusions (P3)

- [x] 8.1 Add `--exclude-disease` CLI flag: filter samples by disease before writing filtered parquet
- [x] 8.2 Add `--min-vaf` CLI flag: filter variants below VAF threshold in both modalities
- [x] 8.3 Add `--min-dp` CLI flag: filter variants below DP threshold in both modalities
- [x] 8.4 Write filtered parquet to `variant_details_filtered/` directory (separate from full dataset)
- [x] 8.5 Log count of excluded variants/samples for each active filter

## 9. Documentation (P3)

- [x] 9.1 Document rescue VCF validation methodology: what "rescued" means, how mismatch is calculated (docs/rescue_validation.md)
- [x] 9.2 Document CxDy tiering: C1-C7 caller tiers, D0-D1 database tiers, how they combine, what each implies (docs/cxdy_tiering.md)
- [x] 9.3 Document BAM column mapping: how BAM pileup columns map to caller VAF/DP fields in validation (docs/bam_column_mapping.md)

## 10. Verification

- [ ] 10.1 Run full pipeline on 8 samples, verify `has_bam_data` is True for samples with pileup data
- [ ] 10.2 Verify `bam_coverage_violin` generates non-empty chart
- [ ] 10.3 Verify chromosome ordering in all CHROM-grouped charts (chr1..22, chrX, chrY)
- [ ] 10.4 Verify depth threshold sweep CSV and chart are generated
- [ ] 10.5 Verify DP distribution charts render correctly
- [ ] 10.6 Run existing test suite, confirm no regressions
- [ ] 10.7 Run with `--exclude-disease X --min-vaf 0.05` and verify filtered parquet directory
