## Why

The statistics and visualization pipeline has 10 confirmed data bugs that produce incorrect results, and 24+ visualization issues that produce broken or incomplete charts. The fundamental architecture is wrong: statistics blindly count without purpose, and charts are single-wise (set-only) instead of providing views across all relevant dimensions. This fix addresses all discovered issues and reorganizes the pipeline around the true goal — identifying thresholds for filtering low-quality variants for downstream ML model training.

## What Changes

### Phase 1: Data bug fixes
- **Insert size**: Add `!is_supplementary() && !is_secondary()` filter in both Rust and pysam BAM stats
- **Coverage**: Accept `--bed` flag for WES exome region denominator
- **Sample summary**: Fix classification counts — use FILTER column (not null VC), include all 6 categories (Somatic/Germline/Reference/Artifact/RNAedit/NoConsensus), remove meaningless PASS stats
- **VC → FILTER**: Replace all "VC" column references with "FILTER" across 5 files; expand domain from 4 to 6 categories with distinct colors
- **BAM pileup**: Integrate into `process_single_sample`, enabled by default, `--no-pileup` to disable
- **Strelka VAF**: Document that Strelka uses tier1-filtered depth denominator (not total depth), add `vaf_denominator` column

### Phase 2: Format migration (CSV → TSV)
- All `.csv` outputs → `.tsv` with tab separator; write_tsv() helper

### Phase 3: Statistics redesign — threshold-focused
- Unified `compute_wise_summary(df, group_cols)` kernel for 6 wises (set/disease/sample/tier/caller/chromosome) + per-caller VAF threshold sweep + filter effectiveness matrix
- Remove meaningless PASS-based stats
- All 6 variant categories (not just Somatic/Germline)

### Phase 4: Visualization redesign — wise-aware factories
- `_plot_bar_wise()`, `_plot_box_violin_wise()`, `_plot_scatter_wise()` factories applied to all 6 wises
- Per-chart fixes: box+violin overlay, percentage annotations, chromosome natural sort, unique colors per tier, handle empty/null data gracefully
- New charts: VAF threshold sweep, caller concordance vs VAF, filter effectiveness heatmap, database enrichment by tier

### Phase 5: Code reorganization
- Output directory restructured by wise (`stats/{wise}/`, `plots/{wise}/`)
- Shared kernel functions in statistics.py, factory functions in visualizer.py
- Remove dead code (duplicate functions, unused aliases)

## Capabilities

### New Capabilities
- `wise-statistics`: Unified multi-dimension aggregation (set/disease/sample/tier/caller/chromosome) with shared metric kernel
- `threshold-analysis`: VAF threshold sweep, filter effectiveness matrix, caller concordance vs quality metrics
- `wise-visualizations`: All chart types parameterized by aggregation dimension, produced for all applicable wises

### Modified Capabilities
- `variant-visualization`: VC→FILTER migration, 6-category domain, box+violin overlays, percentage annotations, chromosome natural sort
- `bam-statistics`: Insert size supplementary filter, WES coverage denominator via --bed

## Impact

- **Python**: `statistics.py` (kernel + threshold), `visualizer.py` (factories + fixes), `cli.py` (--bed, --no-pileup, wise loops), `bam_stats.py` (supplementary filter), `rescue_parser.py` (minor)
- **Rust**: `bam.rs` (supplementary filter, BED denominator), `caller.rs` (document Strelka VAF)
- **Tests**: Extensive updates to cover new wises, threshold analysis, fixed classifications
