## Phase 1: Data Bug Fixes

### 1.1 Insert size
- [x] 1.1.1 Add `!is_supplementary() && !is_secondary()` filter in bam.rs `whole_genome_stats()`
- [x] 1.1.2 Add `not read.is_supplementary and not read.is_secondary` in bam_stats.py `_compute_bam_stats_pysam()`

### 1.2 WES coverage
- [x] 1.2.1 Add `--bed` CLI flag accepting path to BED file
- [x] 1.2.2 Implement BED region length summation for coverage denominator (Python side, post-Rust)
- [x] 1.2.3 Default to whole-genome if no BED provided

### 1.3 Sample summary classification
- [x] 1.3.1 Replace `FILTER == "PASS"` with FILTER-based classification in `sample_summary()`
- [x] 1.3.2 Add all 6 categories: Somatic, Germline, Reference, Artifact, RNAedit, NoConsensus
- [x] 1.3.3 Remove `pass_variants`/`pass_pct` (PASSES_CONSENSUS meaningless)
- [x] 1.3.4 Update `set_summary()` to aggregate all 6 categories

### 1.4 VC → FILTER migration
- [x] 1.4.1 Replace `VC_DOMAIN`/`VC_COLORS` with `CLASSIFICATION_DOMAIN`/`CLASSIFICATION_COLORS` (6 categories)
- [x] 1.4.2 Replace all `"VC" in df.columns` checks with `"FILTER"` in visualizer.py (5 functions)
- [x] 1.4.3 Update `sample_summary()` in statistics.py to use FILTER
- [x] 1.4.4 Update `dataset_summary()` in statistics.py to use FILTER

### 1.5 BAM pileup integration
- [x] 1.5.1 Add pileup call to `process_single_sample()` after tiering
- [x] 1.5.2 Add `--no-pileup` flag (pileup enabled by default)
- [x] 1.5.3 Restore `plot_bam_coverage_violin` chart with pileup data

### 1.6 Strelka VAF documentation
- [x] 1.6.1 Add `vaf_denominator` column to per-caller VAF output: "tier1_depth" for Strelka, "total_depth" for others
- [x] 1.6.2 Add comment in statistics.py documenting the Strelka VAF difference

## Phase 2: CSV → TSV Migration

- [x] 2.1 Add `write_tsv(df, path)` helper in statistics.py
- [x] 2.2 Replace all `.write_csv()` calls with `write_tsv()` in cli.py
- [x] 2.3 Replace `.write_csv()` in bam_stats.py (no .write_csv calls in bam_stats.py — writes are in cli.py)
- [x] 2.4 Add `test_write_tsv_roundtrip` test verifying TSV output reads back correctly

## Phase 3: Statistics Redesign

### 3.1 Wise kernel
- [x] 3.1.1 Define `_WISE_METRICS` dict with all 30+ shared metrics
- [x] 3.1.2 Implement `compute_wise_summary(df, group_cols)` generic kernel
- [x] 3.1.3 Generate set-wise, disease-wise, sample-wise, tier-wise, caller-wise, chromosome-wise summaries
- [x] 3.1.4 Output to `stats/{wise}/` directory structure

### 3.2 Threshold analysis
- [x] 3.2.1 Implement VAF threshold sweep per caller per classification
- [x] 3.2.2 Implement filter effectiveness matrix (FILTER × classification)
- [x] 3.2.3 Output to `stats/threshold/`

### 3.3 Remove PASS-based stats
- [x] 3.3.1 Remove `pass_variants`/`pass_pct` from sample_summary, set_summary, dataset_summary
- [x] 3.3.2 Remove PASS-based charts or replace with classification-based

### 3.4 Per-caller column fix (NEW — root cause of callers-wise zeros + VAF sweep empty)
- [x] 3.4.1 Add 24 per-caller VAF/DP/AD columns to `_CROSS_SAMPLE_COLS` (from 40 → 64 columns)
- [x] 3.4.2 Fix `_ensure_eager` PerformanceWarning: use `df.collect_schema().names()` for LazyFrame column check
- [x] 3.4.3 Verify `compute_caller_wise_summary` returns non-zero counts
- [x] 3.4.4 Verify `compute_vaf_threshold_sweep` returns non-empty DataFrame

## Phase 4: Visualization Redesign

### 4.1 Chart factories
- [x] 4.1.1 Implement `_plot_bar_wise()` — bar chart factory (vc, filter, variant_type, cross_modality)
- [x] 4.1.2 Implement `_plot_box_violin_wise()` — box+violin overlay (vaf, dp)
- [x] 4.1.3 Implement `_plot_scatter_wise()` — scatter factory (dna_vs_rna, ref_alt_dp)
- [x] 4.1.4 Implement `_plot_pie_wise()` — pie chart factory (cosmic_gnomad)

### 4.2 Per-chart fixes — Category A: add `group_col` parameter (10 functions)
- [x] 4.2.1 `plot_vc_distribution`: group_col param already exists, fix title from "Set" to dynamic
- [x] 4.2.2 `plot_caller_overlap`: add group_col param, dynamic title (currently hardcodes "set_number" in 4 places)
- [x] 4.2.3 `plot_variant_type_distribution`: add group_col param, % marks, dynamic title
- [x] 4.2.4 `plot_ti_tv_ratio`: add group_col param, dynamic title
- [x] 4.2.5 `plot_cross_modality`: add group_col param, dynamic title
- [x] 4.2.6 `plot_dna_vs_rna_dp`: add group_col param, dynamic title
- [x] 4.2.7 `plot_ref_alt_dp_scatter`: add group_col param, dynamic title
- [x] 4.2.8 `plot_filter_distribution`: add group_col param, % marks, dynamic title
- [x] 4.2.9 `plot_redi_evidence`: add group_col param, dynamic title
- [x] 4.2.10 `plot_gt_concordance`: add group_col param, per-group concordance computation

### 4.2 Per-chart fixes — Category B: add `color_col` parameter (5 functions)
- [x] 4.2.11 `plot_vaf_distribution`: add color_col param for set/disease/tier coloring
- [x] 4.2.12 `plot_dna_vs_rna_vaf`: add color_col param (currently hardcodes FILTER color)
- [x] 4.2.13 `plot_dna_vs_rna_per_caller`: add color_col param
- [x] 4.2.14 `plot_bam_coverage_violin`: add color_col param for set/disease/tier coloring
- [x] 4.2.15 `plot_caller_concordance_vs_vaf`: add color_col param for set/disease coloring

### 4.2 Per-chart fixes — Category C: add `facet_col` parameter (6 functions)
- [x] 4.2.16 `plot_vaf_boxplot_per_tier`: add facet_col param (facets on top of caller_tier)
- [x] 4.2.17 `plot_dp_boxplot_per_tier`: add facet_col param
- [x] 4.2.18 `plot_gt_concordance_per_tier`: add facet_col param
- [x] 4.2.19 `plot_tiered_caller_overlap`: add facet_col param + % marks + log scale y-axis
- [x] 4.2.20 `plot_tiered_variant_types`: add facet_col param + % marks
- [x] 4.2.21 `plot_tier_quality_distribution`: add facet_col param for set/disease

### 4.2 Per-chart fixes — Category D: quality fixes (no new params)
- [x] 4.2.22 `plot_chromosome_density`: natural sort via when/then chain
- [x] 4.2.23 `plot_tier_quality_distribution`: verified _sample_if_large works
- [x] 4.2.24 `plot_per_sample_distribution`: remove top_n=30 limit, add horizontal scroll for all 65 samples
- [x] 4.2.25 `plot_cosmic_gnomad_annotation`: add group_col param for per-group annotation %

### 4.2 Per-chart fixes — Category E: box+violin overlay (5 functions)
- [x] 4.2.26 `plot_vaf_distribution`: replace mark_boxplot with _plot_box_violin_wise
- [x] 4.2.27 `plot_dp_distribution` (new chart 12): use _plot_box_violin_wise
- [x] 4.2.28 `plot_vaf_boxplot_per_tier`: replace mark_boxplot with _plot_box_violin_wise
- [x] 4.2.29 `plot_dp_boxplot_per_tier`: replace mark_boxplot with _plot_box_violin_wise
- [x] 4.2.30 `plot_bam_coverage_violin`: replace mark_area with _plot_box_violin_wise

### 4.3 New threshold charts (already done)
- [x] 4.3.1 VAF threshold sweep: multi-line retention% vs threshold per caller
- [x] 4.3.2 Caller concordance vs VAF: box plot by # supporting callers
- [x] 4.3.3 Filter effectiveness heatmap: FILTER × Classification
- [x] 4.3.4 Database enrichment by tier: COSMIC/gnomAD % per CxDy tier

### 4.4 New chart functions (5 new)
- [x] 4.4.1 Chart 12: `plot_dp_distribution(df, output_dir, color_col=None)` — per-caller DP box+violin, 4 wises
- [x] 4.4.2 Chart 35: `plot_mean_vaf_per_group(summary_df, output_dir, group_col)` — bar from wise summary, 6 wises
- [x] 4.4.3 Chart 36: `plot_mean_dp_per_group(summary_df, output_dir, group_col)` — bar from wise summary, 6 wises
- [x] 4.4.4 Chart 37: `plot_n_support_callers_dist(df, output_dir, group_col)` — N_SUPPORT_CALLERS histogram, 4 wises
- [x] 4.4.5 Chart 38: `plot_caller_tier_heatmap(df, output_dir, facet_col=None)` — caller×tier matrix, 2 wises
- [x] 4.4.6 Chart 39: `plot_sample_overview_scatter(sample_stats_df, output_dir)` — VAF vs DP per sample, 1 wise

### 4.5 Chart numbering cleanup
- [x] 4.5.1 Consolidate charts 23 and 25 (both "per_tier_vaf") into chart 25, remove chart 23
- [x] 4.5.2 Re-assign chart 12 (was gap) → `plot_dp_distribution`

## Phase 5: CLI Wise Loop & Code Organization

- [x] 5.1 Restructure output directories: `stats/{wise}/`
- [x] 5.2 Remove dead code: duplicate functions, unused aliases, old VC references
- [x] 5.3 Add `--wise` flag to generate specific wises on demand
- [x] 5.4 Define `WISE_CHART_REGISTRY`: wise_name → [(chart_fn, kwargs), ...] with ~100 entries
- [x] 5.5 Implement per-wise chart generation loop in cli.py
- [x] 5.6 Chart functions write to `plots/{wise}/` directory (pass wise plot dir as output_dir)
- [x] 5.7 Generate master dashboard + per-wise dashboards
- [x] 5.8 Create `examples/seq2neo/run_stats.sh` with documented --bed, --wise, --no-pileup usage
- [x] 5.9 Update `_save_chart` to handle wise-specific output directories
- [x] 5.10 Clean up stale `.csv` duplicates (remove, keep only `.tsv` files) — update resume path to only read `.tsv`

## Phase 6: Tests

### 6.1 Existing test updates
- [x] 6.1.1 Update tests for FILTER-based classification (was VC-based)
- [x] 6.1.2 Remove pass_variants assertions
- [x] 6.1.3 Fix plot_per_sample_violin → plot_per_sample_distribution import

### 6.2 Data foundation tests
- [x] 6.2.1 `TestCrossSampleCols::test_per_caller_vaf_columns_included` — 24 columns present
- [x] 6.2.2 `TestCrossSampleCols::test_caller_wise_summary_nonzero` — n_with_vaf > 0 for all callers
- [x] 6.2.3 `TestCrossSampleCols::test_vaf_threshold_sweep_nonempty` — returns non-empty DataFrame
- [x] 6.2.4 `TestCrossSampleCols::test_ensure_eager_no_warning` — no PerformanceWarning

### 6.3 Chart group_col parameterization tests
- [x] 6.3.1 `TestChartGroupCol::test_caller_overlap_with_disease` — group_col="disease_normalized"
- [x] 6.3.2 `TestChartGroupCol::test_variant_type_with_tier` — group_col="final_tier"
- [x] 6.3.3 `TestChartGroupCol::test_ti_tv_with_disease` — group_col="disease_normalized"
- [x] 6.3.4 `TestChartGroupCol::test_cross_modality_with_disease` — group_col="disease_normalized"
- [x] 6.3.5 `TestChartGroupCol::test_filter_dist_with_tier` — group_col="final_tier"
- [x] 6.3.6 `TestChartGroupCol::test_gt_concordance_with_disease` — group_col="disease_normalized"
- [x] 6.3.7 `TestChartGroupCol::test_cosmic_gnomad_with_tier` — group_col="final_tier"
- [x] 6.3.8 `TestChartGroupCol::test_redi_with_tier` — group_col="final_tier"

### 6.4 Chart color_col parameterization tests
- [x] 6.4.1 `TestChartColorCol::test_vaf_distribution_with_disease_color` — color_col="disease_normalized"
- [x] 6.4.2 `TestChartColorCol::test_dna_vs_rna_vaf_with_tier_color` — color_col="final_tier"
- [x] 6.4.3 `TestChartColorCol::test_per_caller_with_set_color` — color_col="set_number"
- [x] 6.4.4 `TestChartColorCol::test_bam_coverage_with_disease_color` — color_col="disease_normalized"

### 6.5 Chart facet_col parameterization tests
- [x] 6.5.1 `TestChartFacetCol::test_vaf_per_tier_with_set_facet` — facet_col="set_number"
- [x] 6.5.2 `TestChartFacetCol::test_dp_per_tier_with_disease_facet` — facet_col="disease_normalized"
- [x] 6.5.3 `TestChartFacetCol::test_gt_per_tier_with_disease_facet` — facet_col="disease_normalized"
- [x] 6.5.4 `TestChartFacetCol::test_caller_overlap_tier_with_set_facet` — facet_col="set_number"
- [x] 6.5.5 `TestChartFacetCol::test_variant_types_tier_with_disease_facet` — facet_col="disease_normalized"

### 6.6 Chart quality tests
- [x] 6.6.1 `TestChartQuality::test_percentage_marks_present` — chart has text marks
- [x] 6.6.2 `TestChartQuality::test_log_scale_on_caller_overlap` — y-axis scale is log
- [x] 6.6.3 `TestChartQuality::test_per_sample_no_top_n_limit` — all 65 samples shown
- [x] 6.6.4 `TestChartQuality::test_box_violin_overlay` — chart has both box and violin layers

### 6.7 New chart function tests
- [x] 6.7.1 `TestNewCharts::test_dp_distribution_returns_chart` — chart 12
- [x] 6.7.2 `TestNewCharts::test_mean_vaf_per_group_returns_chart` — chart 35
- [x] 6.7.3 `TestNewCharts::test_mean_dp_per_group_returns_chart` — chart 36
- [x] 6.7.4 `TestNewCharts::test_n_support_callers_dist_returns_chart` — chart 37
- [x] 6.7.5 `TestNewCharts::test_caller_tier_heatmap_returns_chart` — chart 38
- [x] 6.7.6 `TestNewCharts::test_sample_overview_scatter_returns_chart` — chart 39

### 6.8 CLI wise loop tests
- [x] 6.8.1 `TestCLIWiseLoop::test_wise_flag_parses_correctly` — --wise set disease
- [x] 6.8.2 `TestCLIWiseLoop::test_wise_directories_created` — plots/set/, plots/disease/, etc.
- [x] 6.8.3 `TestCLIWiseLoop::test_wise_flag_default_all` — no --wise → all 6 wises
- [x] 6.8.4 `TestCLIWiseLoop::test_wise_chart_count` — verify chart file count per wise

### 6.9 TSV tests
- [x] 6.9.1 `TestTSV::test_write_tsv_roundtrip` — write→read produces identical DataFrame
- [x] 6.9.2 `TestTSV::test_resume_reads_tsv_not_csv` — resume path reads .tsv files

### 6.10 Regression
- [x] 6.10.1 Run full test suite, verify 0 failures, 0 errors (114 passed, 0 failures)

## Phase 7: Verification

- [ ] 7.1 Run 12-sample pipeline, verify sample_summary counts use 6-category FILTER
- [ ] 7.2 Verify all 6 wise directories created with charts
- [ ] 7.3 Verify ~100 chart files generated (count across all wise dirs)
- [ ] 7.4 Verify caller-wise stats non-zero (per-caller columns present in _CROSS_SAMPLE_COLS)
- [ ] 7.5 Verify VAF threshold sweep non-empty
- [ ] 7.6 Verify BAM pileup columns present in parquet
- [ ] 7.7 Verify all charts render (HTML inspection — no white pages)
- [ ] 7.8 Verify Strelka VAF documentation in output
- [ ] 7.9 Verify no stale .csv duplicates (only .tsv files in output)
- [ ] 7.10 Verify shell script `run_stats.sh` works with --bed, --wise, --no-pileup flags
