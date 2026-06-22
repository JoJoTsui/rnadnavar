## Context

Round 3 of stats/viz fixes based on full 61-sample pipeline run review by 3 independent agents. The architecture is stable (Rust BAM, no pysam fallback, manifest-first paths). These fixes target chart readability, data completeness, and downstream ML guidance.

## Goals / Non-Goals

**Goals:**
- Fix all 7 chart rendering bugs (chromosome ordering, DP capping, color mismatch, coverage violin, log scale, set gridding, DP sweep classification)
- Add ML threshold guidance cross-tabulations for train/val/test partition decisions
- Add somatic modality sub-classification derived from existing C-tier (not a parallel system)
- Add global scientific publishing altair theme

**Non-Goals:**
- Refactoring the wise-based chart architecture
- Adding new variant callers or tiering logic
- Automatic variant filtering (statistics guide thresholds; users decide)
- Switching from altair to seaborn/matplotlib

## Decisions

### D1: Chromosome sort — pass explicit `sort=` list to altair
**Choice:** After calling `_sort_chromosomes(counts, group_col)`, extract `counts[group_col].to_list()` and pass as `sort=` parameter to `alt.X()`. Also apply to TSV output.
**Rationale:** Altair ignores DataFrame row order for `:N` (nominal) type. The `sort=` parameter is the only way to control axis order. This matches the working pattern in `plot_chromosome_density` (line 449).

### D2: DP capping — domain=[0, 2000] with outlier annotation
**Choice:** Use `scale=alt.Scale(domain=[0, 2000])` on DP y-axes. Add subtitle annotation showing count of excluded outliers.
**Alternative:** Log scale everywhere. Rejected for boxplots because IQR interpretation is harder on log scale.
**Exception:** `plot_bam_dp_distribution` uses `symlog` scale since it's a violin+box overlay where log-density distortion matters.

### D3: Somatic modality — derive from C-tier, not parallel system
**Choice:** Map existing `caller_tier` to `somatic_modality` for Somatic variants: C1→MultiModality, C2/C5→DNA_only, C3/C6→RNA_only, C4/C7→Weak.
**Rationale:** Agent 3 identified the C-tier already encodes per-modality caller consensus. Building a parallel system from DP non-null counts would create confusion and potential inconsistency.

### D4: ML cross-tabulation — partition-aware statistics
**Choice:** Add a `partition` computed column (chr1→test, chr21-22→val, rest→train). Compute FILTER × VAF_bin × DP_bin × partition cross-tabs. Output as TSV.
**Rationale:** The user needs statistics to guide filtering thresholds, not automatic filtering. These tables answer: "If I keep only Somatic AND require VAF≥X AND DP≥Y, how many variants remain per partition?"

### D5: Scientific theme — global altair registration
**Choice:** Register a "publishing" theme via `alt.themes.register()`, enable via `--theme publishing` CLI flag.
**Rationale:** Applies globally to all charts with zero per-chart changes. Theme includes: white background, 11pt labels, 13pt titles, minimal grid, Arial font.

## Risks / Trade-offs

- **[Risk]** DP capping at 2000 hides 0.45% of values → mitigated by outlier count annotation
- **[Risk]** `symlog` scale for BAM DP violin may not render correctly in all browsers → fallback to linear with capping
- **[Risk]** Somatic modality derived from C-tier may not exactly match DP-based presence → documented as design choice
