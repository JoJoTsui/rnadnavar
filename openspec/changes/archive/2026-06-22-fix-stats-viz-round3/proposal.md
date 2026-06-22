## Why

The first full 61-sample pipeline run with fix-stats-viz-round2 revealed 11 remaining issues across statistics computation, visualization rendering, and downstream ML guidance. Three independent review agents confirmed all findings with code-level, data-level, and design-level evidence. The most critical issue (P0) is that chromosome ordering is broken in 7 of 8 chart sites despite the sorting function existing — altair ignores DataFrame row order unless explicitly told via `sort=`. Five P1 bugs make charts unreadable (DP outliers compress the IQR, violin/boxplot color mismatch, coverage violin appears blank). Four P2 features are needed to guide downstream variant filtering and ML dataset splitting.

## What Changes

### Bug Fixes (P0)
- **Chromosome ordering — pass `sort=` to altair X encoding**: 7 chart functions call `_sort_chromosomes()` but don't pass `sort=<ordered_list>` to altair's `alt.X()`. Also fix chromosome summary TSV output ordering.

### Bug Fixes (P1)
- **DP threshold sweep — add per-classification breakdown**: Replicate VAF sweep's per-FILTER inner loop for the DP sweep
- **Per-sample charts — grid by set_number**: Add set faceting to charts 08 and 22; add `set_number` to `sample_tier_summary` group-by
- **BAM coverage violin — clip depth outliers**: Cap depth at 2000 before `transform_density` (only 0.45% of values exceed 2000)
- **BAM DP distribution — log/symlog scale + boxplot color**: Add `type="symlog"` to y-axis; fix violin/boxplot color mismatch in 3 functions
- **VAF/DP coloring — systematic violin/boxplot mismatch**: Boxplot encoding lacks `color` in `plot_vaf_distribution`, `plot_dp_distribution`, `plot_bam_dp_distribution`
- **Depth capping at 2000**: Apply `scale=alt.Scale(domain=[0, 2000])` or clip data in 7 DP chart functions; annotate excluded outlier count

### New Features (P2)
- **ML threshold guidance statistics**: Combined FILTER × VAF_bin × DP_bin cross-tabulation per chromosome partition (chr1=test, chr21-22=val, rest=train); disease × partition cross-tab; low-VAF RNA support stats
- **FP cross-tabulation**: Non-Somatic variants: FILTER × N_SUPPORT_CALLERS for removal guidance
- **Somatic modality sub-classification**: Derive `somatic_modality` from existing `caller_tier` (C1→MultiModality, C2/C5→DNA_only, C3/C6→RNA_only, C4/C7→Weak) — NOT a parallel system
- **Scientific publishing altair theme**: Global theme via `alt.themes.register()` with `--theme` CLI flag

## Capabilities

### New Capabilities
- `ml-threshold-guidance`: Combined cross-tabulation statistics and visualizations to guide variant filtering thresholds and train/val/test chromosome partition decisions
- `somatic-modality-subclass`: Derive somatic sub-classification (DNA_only, RNA_only, MultiModality, Weak) from existing CxDy C-tier caller consensus
- `scientific-publishing-theme`: Global altair theme registration for publication-quality chart output

### Modified Capabilities
- `variant-visualization`: Fix chromosome ordering (pass `sort=` to altair), DP capping, log scale, violin/boxplot color mismatch, per-sample set gridding, depth outlier annotation
- `depth-threshold-sweep`: Add per-classification (FILTER) breakdown to DP threshold sweep, matching VAF sweep pattern

## Impact

- `bin/vcf_stats/seq2neo/visualizer.py` — chromosome sort fix (7 sites), DP capping (7 functions), log scale, color mismatch fix (3 functions), per-sample gridding, coverage violin clip, theme registration (~200 lines)
- `bin/vcf_stats/seq2neo/statistics.py` — DP sweep classification, somatic_modality column, ML cross-tabulation, FP cross-tab (~150 lines)
- `bin/vcf_stats/seq2neo/cli.py` — wire new stats/charts, `--theme` flag, chromosome TSV sort (~50 lines)
- `bin/vcf_stats/docs/` — somatic modality docs, ML guidance docs
