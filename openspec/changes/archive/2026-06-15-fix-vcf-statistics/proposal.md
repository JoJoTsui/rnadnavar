## Why

The current seq2neo VCF statistics pipeline has multiple bugs (dashboard showing only one plot, missing color legends) and analytical gaps (no variant tiering integration, insufficient per-modality BAM statistics, missing per-tier VAF/DP/GT analysis). This needs a comprehensive fix to produce correct, actionable variant statistics with proper tiering integration from the existing `bin/vcf_stats/tiering.py` module.

## What Changes

### Bug Fixes
- Fix dashboard.html concatenating full HTML documents (each `to_html()` includes doctype/head/body) — only first chart renders. Switch to extracting body content or using `to_html(inline=True)` style embedding.
- Fix VC distribution chart missing color bar legend

### Tiering Integration (from bin/vcf_stats/tiering.py)
- Integrate the CxDy hybrid tiering system (C1-C7 caller support × D0-D1 database evidence) into variant statistics
- Re-express caller support distribution as tier-aware analysis using FILTERS_NORMALIZED / FILTER_NORMALIZED_* fields and the TieringEngine
- Compute per-tier VAF, DP, GT statistics with violin plots
- Add REF_DP and ALT_DP statistics (currently only DP and VAF are computed)

### Analysis Levels
- Statistics at three levels: **per-sample**, **per-tier**, and **whole-dataset**
- Overall statistics aggregating across all tiers alongside per-tier breakdowns

### Visualization Fixes & New Charts
- Per-sample variant plot: x-axis = sample ID (not disease)
- VAF distribution: replace box plots with violin plots, split by tier
- GT concordance: per-tier analysis
- Cross-modality plot: add percentage annotations
- Add per-sample variant distribution plot

### BAM Statistics
- Per-sample per-modality BAM statistics (read counts, mapping rates, coverage)
- Per-variant BAM stats appended to `variant_details.parquet` (per-variant read depth, strand bias, etc. from caller FORMAT fields)

### Performance (Future)
- Evaluate polars-based approach; if speed insufficient, consider maturin + Rust for heavy computation loops

## Capabilities

### New Capabilities
- `variant-tiering-stats`: Variant statistics organized by CxDy tiering levels with per-tier VAF/DP/GT violin plots
- `bam-statistics`: Per-sample per-modality BAM statistics and per-variant BAM stats in parquet output
- `ref-alt-dp-stats`: REF_DP and ALT_DP statistics alongside existing DP/VAF metrics
- `multi-level-aggregation`: Statistics at per-sample, per-tier, and whole-dataset levels

### Modified Capabilities
- `variant-visualization`: Dashboard fix, VC color bar fix, violin plots, per-tier breakdown, sample-ID x-axis, cross-modality percentages, per-sample distribution

## Impact

- **bin/vcf_stats/seq2neo/visualizer.py**: Restructure dashboard to work with altair's full-page output; fix color bar; add violin plots; tier-aware charts
- **bin/vcf_stats/seq2neo/statistics.py**: Add tier-aware statistics, REF_DP/ALT_DP stats, multi-level aggregation
- **bin/vcf_stats/seq2neo/caller_parser.py**: Extract REF_DP and ALT_DP (already partially done via AD_REF/AD_ALT, needs REF_DP)
- **bin/vcf_stats/seq2neo/rescue_parser.py**: Parse FILTERS_NORMALIZED and FILTER_NORMALIZED_* fields for tiering
- **New file: bin/vcf_stats/seq2neo/tiering_stats.py**: Integration layer between tiering.py and statistics module
- **New file: bin/vcf_stats/seq2neo/bam_stats.py**: BAM statistics computation module
- **bin/vcf_stats/seq2neo/cli.py**: Wire new modules into pipeline
- **bin/vcf_stats/tests/test_seq2neo_stats.py**: Update tests for all changes
