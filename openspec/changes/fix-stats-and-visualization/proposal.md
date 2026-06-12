## Why

The current statistics and visualization codebase has 28 confirmed issues across BAM stats (incorrect sampling, missing filters), statistics computation (duplicate columns, missing metrics, wrong caller overlap), and visualization (dead charts, wrong sectioning, CDN dependency, unreadable plots). Additionally, key biological analyses are missing: cross-modality DNA↔RNA comparison, per-caller agreement matrix, chromosome-level variant density, REDIportal visualization. The goal is to produce statistics that enable determining quality thresholds for a 3-category (Germline/Reference/Somatic) variant prediction model.

## What Changes

### 1. BAM Statistics — Fix All Metrics
- Remove 1M-read sampling cap; report actual total reads
- Add `read.is_proper_pair` filter for insert_size calculation
- Compute actual mean coverage from BAM header reference lengths × total bases
- Add DNA normal (DN) BAM to whole-genome stats
- Document RNA 100% mapping rate (STAR behavior)

### 2. Statistics — Fix + Extend
- Replace `caller_overlap` (N_SUPPORT_CALLERS) with C1-C7 tiering distribution throughout
- Remove duplicate `DNA_REF_mean`/`DNA_ALT_mean` columns (keep `DNA_REF_DP_mean`/`DNA_ALT_DP_mean`)
- Add per-FILTER, per-tier, per-variant-type breakdowns to `set_summary()` and `disease_summary()`
- Include MNV/INDEL counts alongside Ti/Tv ratio
- Replace `.to_dicts()` iterations with polars-native operations for GT concordance
- Add per-caller pairwise agreement matrix (6×6)
- Extract caller FORMAT strand bias (SB) and use for cross-validation against BAM pileup

### 3. Visualization — Fix + Add New Charts
- **Fix blank dashboard**: bundle vega libraries inline instead of loading from CDN
- **Fix sectioning**: dynamic section assignment based on chart category tags
- **Fix dead charts**: `plot_bam_coverage_violin` needs BAM_DP_* columns wired; `plot_per_tier_cross_sample_vaf` needs VAF columns in sample_tier_summary
- **Fix duplicates**: remove second `plot_ref_alt_dp_scatter` call, remove dead `_to_pandas()`
- **Fix names**: rename "violin" functions to "boxplot" or implement actual violin
- **Fix readability**: top-N filtering for per-sample charts, facet by set instead of 65 bars
- **New charts**: per-caller agreement matrix heatmap, chromosome variant density (Manhattan), REDIportal evidence distribution, DNA vs RNA per-caller VAF at shared positions, FILTER distribution pie, BAM insert size distribution

### 4. Cross-Modality DNA↔RNA Comparison
- Replace rescue validation (DP_MEAN vs caller mean) with DNA↔RNA modality comparison
- Per-variant: DNA caller VAF vs RNA caller VAF at positions detected by both modalities
- Per-sample: DNA vs RNA mean DP/VAF correlation
- New visualizations: DNA vs RNA VAF scatter per caller, modality concordance heatmap

### 5. Code Cleanup
- Merge duplicate `validate_all_samples` functions
- Wire `bam_validation.py` into CLI
- Remove unused `_to_pandas()` helper
- Remove DNA_REF_mean/DNA_ALT_mean duplicate columns

## Capabilities

### New Capabilities
- `cross-modality-comparison`: DNA↔RNA VAF/DP comparison at shared variant positions
- `caller-agreement-matrix`: 6×6 pairwise caller agreement heatmap
- `chromosome-density`: Manhattan-style variant density per chromosome
- `bam-quality-metrics`: Insert size distribution, proper pair rate, coverage depth

### Modified Capabilities
- `bam-statistics`: Fix sampling, insert_size, coverage; add DN BAM type
- `four-level-statistics`: Replace N_SUPPORT_CALLERS with C1-C7 tiering; extend set/disease summaries
- `variant-visualization`: Fix all 28 confirmed issues; add 8 new chart types
- `rescue-validation`: Replace INFO mean validation with DNA↔RNA cross-modality comparison

## Impact
- `bam_stats.py`: sampling cap removal, proper_pair filter, coverage fix
- `statistics.py`: tiering replacement, remove duplicates, extend summaries
- `visualizer.py`: CDN inline, section fix, 8 new charts, dead code removal
- `rescue_validator.py`: replace validation metrics with cross-modality comparison
- `bam_validation.py`: wire into CLI
- `cli.py`: fix duplicate chart, wire new outputs, replace .to_dicts()
- `rust_vcf.py`: no changes needed (type casting already fixed)
