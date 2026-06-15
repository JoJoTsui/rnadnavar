## Why

The statistics and visualization pipeline has 10 confirmed data bugs that produce incorrect results, and 30+ chart functions that only generate set-wise views instead of providing insights across all six relevant aggregation dimensions (set, disease, sample, tier, caller, chromosome). After running the pipeline on 65 samples, deep review identified three categories of gaps:

1. **Silent data failures**: `_CROSS_SAMPLE_COLS` prunes 24 per-caller VAF/DP/AD columns, causing caller-wise stats (all zeros), VAF threshold sweep (empty), and multiple charts to silently produce no data.

2. **Single-dimension charts**: All 30 chart types only generate set-wise views. The design calls for ~105 chart×wise combinations (6 wises × ~18 chart types each), but currently only ~30 flat set-wise charts exist.

3. **Chart quality gaps**: No percentage text marks, no box+violin overlays, hardcoded "set_number" in 12 functions, missing chart types (DP distribution, caller-tier heatmap, etc.), and a numbering gap at chart 12.

This update captures every required change — no deferrals, no "Phase 5 later" items. Every task must be implemented.

## What Changes

### Part 1: Data Foundation Fix
- **Per-caller columns**: Add 24 per-caller VAF/DP/AD columns to `_CROSS_SAMPLE_COLS` (from 40 → 64 columns). This single fix unblocks caller-wise stats, VAF threshold sweep, and 6+ chart functions that currently produce empty output.

### Part 2: Chart Function Parameterization (14 functions)
- **10 functions**: Add `group_col` parameter to replace hardcoded "set_number" (caller_overlap, variant_type_distribution, ti_tv_ratio, cross_modality, dna_vs_rna_dp, ref_alt_dp_scatter, filter_distribution, redi_evidence, gt_concordance, cosmic_gnomad_annotation)
- **5 functions**: Add `color_col` parameter for wise-based coloring (vaf_distribution, dna_vs_rna_vaf, dna_vs_rna_per_caller, bam_coverage_violin, caller_concordance_vs_vaf)
- **6 functions**: Add `facet_col` parameter for wise-based faceting on top of existing dimensions (vaf_boxplot_per_tier, dp_boxplot_per_tier, gt_concordance_per_tier, tiered_caller_overlap, tiered_variant_types, tier_quality_distribution)

### Part 3: Chart Quality Fixes (5 functions)
- **Percentage text marks**: tiered_caller_overlap, tiered_variant_types, filter_distribution, variant_type_distribution
- **Log scale**: tiered_caller_overlap y-axis
- **Per-sample**: Remove top_n=30 limit, add horizontal scroll for all 65 samples
- **Box+violin overlays**: Apply `_plot_box_violin_wise` factory to VAF distribution, DP distribution, per-tier VAF, per-tier DP, BAM coverage violin

### Part 4: New Chart Functions (5 functions)
- **Chart 12**: `plot_dp_distribution` — per-caller DP box+violin (4 wises)
- **Chart 35**: `plot_mean_vaf_per_group` — mean VAF bar from wise summary (6 wises)
- **Chart 36**: `plot_mean_dp_per_group` — mean DP bar from wise summary (6 wises)
- **Chart 37**: `plot_n_support_callers_dist` — N_SUPPORT_CALLERS histogram per group (4 wises)
- **Chart 38**: `plot_caller_tier_heatmap` — caller × tier detection matrix (2 wises)
- **Chart 39**: `plot_sample_overview_scatter` — mean VAF vs mean DP per sample, colored by set/disease (1 wise)

### Part 5: CLI Wise Loop
- `WISE_CHART_REGISTRY`: wise_name → [(chart_fn, kwargs), ...] with ~105 entries mapping every chart to every applicable wise
- Per-wise chart generation loop creating `plots/{wise}/` directory structure
- `--wise` flag controls which wises are generated (default: all 6)
- Charts 23 and 25 consolidated (both "per_tier_vaf", keep chart 25)
- Chart 12 filled (DP distribution)

### Part 6: Shell Script & Documentation
- Create `examples/seq2neo/run_stats.sh` with documented `--bed`, `--wise`, `--no-pileup` usage
- Update `README.usage.md` with all new CLI flags

### Part 7: Tests (~55 new methods)
- `TestCrossSampleCols`: verify 64 columns include per-caller VAF/DP
- `TestChartGroupCol`: each chart × each applicable group_col
- `TestChartColorCol`: coloring with each color_col
- `TestChartFacetCol`: faceting with each facet_col
- `TestChartQuality`: % marks, log scale, horizontal scroll
- `TestNewCharts`: 5 new chart functions
- `TestCLIWiseLoop`: --wise flag, per-wise directories
- `TestTSVRoundtrip`: write_tsv + read_csv(separator="\t")
- All existing tests pass, 0 regressions

### Part 8: Verification
- Run 12-sample pipeline, verify all outputs correct
- Verify ~105 chart files generated across 6 wise directories
- Verify caller-wise stats non-zero
- Verify VAF threshold sweep non-empty
- Verify BAM pileup columns present
- Verify all charts render (HTML inspection)
- Verify Strelka VAF documentation in output

## Impact

- **Python**: `statistics.py` (per-caller columns, wise kernel), `visualizer.py` (14 chart function signatures + 5 new functions + quality fixes), `cli.py` (wise registry + per-wise loop), `bam_stats.py` (already done)
- **Tests**: ~55 new test methods, updated fixtures
- **Docs**: New shell script, updated README
- **Output**: ~105 chart×wise combinations across `plots/{set,disease,sample,tier,caller,chromosome}/`
