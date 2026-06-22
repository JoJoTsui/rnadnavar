## Phase 1: Infrastructure Helpers

- [x] 1.1 Define 10 global color constant pairs at module top of visualizer.py:
  - `CALLER_DOMAIN` / `CALLER_COLORS` (6 callers — DNA/RNA paired colors)
  - `RESCUE_DOMAIN` / `RESCUE_COLORS` (YES/NO — green/red)
  - `VARIANT_TYPE_DOMAIN` / `VARIANT_TYPE_COLORS` (SNV/INS/DEL/MNV)
  - `SOMATIC_MODALITY_DOMAIN` / `SOMATIC_MODALITY_COLORS` (MultiModality/DNA_only/RNA_only/Weak/Unknown)
  - `BAM_TYPE_DOMAIN` / `BAM_TYPE_COLORS` (DN/DT/RT)
  - `MODALITY_DOMAIN` / `MODALITY_COLORS` (DNA/RNA)
  - `AGREEMENT_DOMAIN` / `AGREEMENT_COLORS` (2/3/4)
  - Keep existing `CLASSIFICATION_DOMAIN` / `CLASSIFICATION_COLORS`
  - Remove `_RESCUE_COLORS` dict (replaced by RESCUE_DOMAIN/RESCUE_COLORS lists)
  - Remove inline `modality_domain`/`modality_colors` from `plot_somatic_modality_pie`

- [ ] 1.2 Create `_color_scale(entity)` helper:
  - Accepts entity string (e.g., "caller", "FILTER", "RESCUED", "variant_type", "bam_type", "somatic_modality", "modality")
  - Returns `alt.Scale(domain=..., range=...)` from global constants
  - Fallback for unknown entities: `alt.Scale(scheme="category10")`

- [ ] 1.3 Create `_apply_faceting(chart, group_col, columns=None)` helper:
  - `set_number` → `.facet(columns=2).resolve_scale(x="independent", y="shared")`
  - `disease_normalized` → `.facet(columns=4).resolve_scale(y="shared")`
  - `sample_id` → `.facet(columns=2).resolve_scale(x="independent", y="shared")`
  - `FILTER` → `.facet(columns=3).resolve_scale(y="shared")`
  - Override `columns` if explicitly provided
  - Default → `.facet(columns=3)`

- [ ] 1.4 Create `_add_heatmap_text(base, x_enc, y_enc, text_col, pdf, fontSize=8)` helper:
  - Computes median of text_col for conditional coloring (white on dark, black on light)
  - Returns `mark_text()` chart layer ready to be combined with `rect + text`
  - Default fontSize=8, format=",d" for integers

- [ ] 1.5 Create `_add_bar_labels(chart, pdf, x_enc, y_enc, label_col, fmt=",d", dy=-8, fontSize=9)` helper:
  - Returns `mark_text()` layer for stacked or grouped bar charts
  - Handles stack offset encoding

- [ ] 1.6 Create `_clip_dp(pdf, col, cap=2000)` helper:
  - Clips column values to [0, cap]
  - Returns `(clipped_pdf, n_over)` tuple for subtitle annotation

- [ ] 1.7 Create `_count_scale()` helper:
  - Returns `alt.Scale(type="symlog", constant=1)` for count axes

## Phase 2: Color Unification

- [ ] 2.1 Replace all caller color encodings with `_color_scale("caller")`:
  - `plot_vaf_distribution` (lines 650, 656, 666 — 3 sites)
  - `plot_vaf_boxplot_per_tier` (line 762)
  - `plot_dp_boxplot_per_tier` (line 798)
  - `plot_dp_distribution` (lines 1572, 1577)
  - `plot_vaf_threshold_sweep` (line 1381)
  - `plot_dp_threshold_sweep` (line 1411)

- [ ] 2.2 Replace all classification/FILTER color encodings with `_color_scale("FILTER")`:
  - `plot_filter_distribution` (line 473 — currently uses default, NO scale)
  - Verify existing uses in plot_vc_distribution, plot_dna_vs_rna_vaf, plot_dna_vs_rna_per_caller, plot_low_vaf_rna_support already use constants

- [ ] 2.3 Replace all rescue color encodings with `_color_scale("RESCUED")`:
  - `plot_rescue_breakdown` (line 2033-2034)
  - `plot_rescue_by_filter` (line 2056-2057)
  - `plot_rescue_vaf_boxplot` (line 2116)
  - `plot_rescue_dp_boxplot` (line 2145)
  - `plot_rescue_sample_distribution` (line 2166-2167)
  - `plot_rescue_by_tier` (line 2188-2189)
  - `plot_rescue_caller_support` (line 2224-2225)
  - `plot_cross_modality` (line 437 — currently uses default, NO scale)

- [ ] 2.4 Replace remaining ad-hoc color encodings:
  - `plot_variant_type_distribution` (line 372, 376) → `_color_scale("variant_type")`
  - `plot_tiered_variant_types` (line 575) → `_color_scale("variant_type")`
  - `plot_somatic_modality_pie` (lines 1863-1864) → `_color_scale("somatic_modality")`
  - `plot_somatic_modality_bars` (line 1895) → `_color_scale("somatic_modality")`
  - `plot_bam_metrics_bars` (line 1235) → `_color_scale("bam_type")`
  - `plot_bam_coverage_violin` (line 1323) → `_color_scale("bam_metric")` or category10
  - `plot_bam_dp_distribution` (line 1988) → same as above
  - `plot_per_sample_tier_distribution` (line 1272) → `_color_scale("final_tier")` or category10
  - `plot_per_tier_vaf_boxplot` (line 851) → `_color_scale("final_tier")`
  - `plot_per_tier_dp_boxplot` (line 2009) → `_color_scale("final_tier")`
  - `plot_mean_vaf_per_group` (line 1603) → `_color_scale("modality")`
  - `plot_mean_dp_per_group` (line 1636) → `_color_scale("modality")`

## Phase 3: Faceting Migration

- [ ] 3.1 Migrate 11 functions from `alt.Column()` in `.encode()` to `_apply_faceting()`:
  - `plot_caller_overlap` (line 347)
  - `plot_tiered_caller_overlap` (line 538)
  - `plot_tiered_variant_types` (line 568)
  - `plot_tier_quality_distribution` (line 597)
  - `plot_vaf_boxplot_per_tier` (line 754)
  - `plot_dp_boxplot_per_tier` (line 792)
  - `plot_gt_concordance` (line 1077)
  - `plot_gt_concordance_per_tier` (line 1141)
  - `plot_n_support_callers_dist` (line 1666)
  - `plot_caller_tier_heatmap` (line 1720)
  - `plot_bam_dp_distribution` (line 1984)

- [ ] 3.2 Update 6 existing `.facet()` calls to use `_apply_faceting()`:
  - `plot_per_sample_distribution` (line 1195) — add `x="independent"`
  - `plot_bam_metrics_bars` (line 1245) — add `x="independent"`
  - `plot_per_sample_tier_distribution` (line 1279) — change to `x="independent"`
  - `plot_rescue_sample_distribution` (line 2172) — add `x="independent"`
  - `plot_vaf_distribution` faceted path (line 674) — add `columns=` for disease wrapping
  - `plot_caller_concordance_vs_vaf` (line 1440) — standardize

- [ ] 3.3 Add disease-wise 4×3 grid support:
  - Ensure `_apply_faceting` applies `columns=4` when `group_col="disease_normalized"`
  - Test with disease-wise charts that currently produce single-row horizontal layouts

## Phase 4: Violin → Boxplot

- [ ] 4.1 Remove `transform_density` + `mark_area` from `_plot_box_violin_wise` — keep boxplot only
- [ ] 4.2 Remove violin from `plot_vaf_distribution` (both non-faceted and faceted paths) — keep boxplot only
- [ ] 4.3 Remove violin from `plot_bam_coverage_violin` — rename to `plot_bam_coverage_boxplot`, keep boxplot
- [ ] 4.4 Remove violin from `plot_dp_distribution` (non-faceted path) — keep boxplot only

## Phase 5: Per-FILTER Sub-Plots + Variant-Category-Wise

- [ ] 5.1 In `plot_vaf_threshold_sweep`: remove `classification.is_null()` filter, add `color=alt.Color("classification:N", scale=_color_scale("FILTER"))`, use row faceting for per-classification panels
- [ ] 5.2 In `plot_dp_threshold_sweep`: same pattern as 5.1
- [ ] 5.3 In `plot_dna_vs_rna_per_caller`: add optional FILTER faceting via `alt.Row("FILTER:N")` when FILTER column present
- [ ] 5.4 In `plot_rescue_vaf_boxplot`: add optional FILTER faceting
- [ ] 5.5 In `plot_rescue_dp_boxplot`: add optional FILTER faceting
- [ ] 5.6 Add `("variant-category", ["FILTER"])` to `wise_configs` in cli.py
- [ ] 5.7 Add `"variant-category"` to `all_wise_names` in cli.py
- [ ] 5.8 Add `"variant-category"` registry entry to `_WISE_CHART_REGISTRY` with 11 charts:
  - plot_vc_distribution, plot_variant_type_distribution, plot_ti_tv_ratio, plot_cross_modality, plot_filter_distribution, plot_redi_evidence, plot_cosmic_gnomad_annotation, plot_vaf_distribution, plot_caller_concordance_vs_vaf, plot_caller_agreement_matrix, plot_tier_quality_distribution

## Phase 6: Symlog Count Axes

- [ ] 6.1 Apply `_count_scale()` to all count-bearing chart y-axes (~20 functions):
  - `plot_vc_distribution`, `plot_variant_type_distribution`, `plot_filter_distribution`
  - `plot_chromosome_density`, `plot_redi_evidence`, `plot_tiered_variant_types`
  - `plot_tier_quality_distribution`, `plot_per_sample_distribution`
  - `plot_per_sample_tier_distribution`, `plot_gt_concordance`, `plot_gt_concordance_per_tier`
  - `plot_low_vaf_rna_support`, `plot_n_support_callers_dist`, `plot_cross_modality`
  - `plot_rescue_breakdown`, `plot_rescue_by_filter`, `plot_rescue_sample_distribution`
  - `plot_rescue_by_tier`, `plot_rescue_caller_support`, `plot_somatic_modality_bars`

## Phase 7: Bug Fixes

- [ ] 7.1 Fix `plot_rescue_vaf_boxplot`: normalize RESCUED to YES/NO BEFORE `drop_nulls()`, then `drop_nulls(subset=vaf_cols)` only
- [ ] 7.2 Fix `plot_rescue_dp_boxplot`: same pattern as 7.1
- [ ] 7.3 Fix `plot_bam_dp_distribution`: add `_clip_dp()` before sampling (missed in round 4)

## Phase 8: Text Labels

- [ ] 8.1 Add `_add_heatmap_text()` to `plot_caller_agreement_matrix` (chart 28) — text_col="pct", format=".1f"
- [ ] 8.2 Add `_add_heatmap_text()` to `plot_validation_heatmap` (chart 13) — text_col="mismatch_pct", format=".1f"
- [ ] 8.3 Add `_add_heatmap_text()` to `plot_filter_effectiveness_heatmap` (chart 33) — text_col="pct_flagged", format=".1f"
- [ ] 8.4 Add `_add_heatmap_text()` to `plot_caller_tier_heatmap` (chart 38) — text_col="pct_detected", format=".1f"
- [ ] 8.5 Update `plot_fp_cross_tab_heatmap` (chart 45): reduce fontSize from 9→7
- [ ] 8.6 Update `plot_rescue_cross_tab_heatmap` (chart 50): reduce fontSize from 9→8
- [ ] 8.7 Include Somatic in `compute_fp_cross_tab()`: remove `FILTER != "Somatic"` filter, rename chart to "Classification × Caller Support"
- [ ] 8.8 Add `_add_bar_labels()` to `plot_rescue_by_tier` (chart 54)
- [ ] 8.9 Add `_add_bar_labels()` to `plot_rescue_caller_support` (chart 56)

## Phase 9: Resume Fallback + Registry Cleanup

- [ ] 9.1 In cli.py `--resume` path: if `sample_summary.tsv` doesn't exist, recompute from parquet files
- [ ] 9.2 In cli.py `--resume` path: if `bam_stats.tsv` doesn't exist, skip BAM charts gracefully (don't fail)
- [ ] 9.3 Remove duplicate chart entries from global chart list that are already in per-wise registry (plot_caller_agreement_matrix, plot_tier_quality_distribution)

## Phase 10: Verification

- [ ] 10.1 Run pipeline with `--resume`, verify no `Facet has no parameter named 'columns'` errors
- [ ] 10.2 Verify all callers use the same colors across set-wise, disease-wise, and global charts
- [ ] 10.3 Verify all classification/FILTER uses same colors across all charts
- [ ] 10.4 Verify rescue YES/NO uses same green/red across all rescue charts AND chart 11
- [ ] 10.5 Verify no violin plots have empty density areas — all should show boxplots
- [ ] 10.6 Verify threshold sweep charts show per-classification lines (Somatic, Germline, Reference)
- [ ] 10.7 Verify variant-category-wise directory exists with ~11 chart outputs
- [ ] 10.8 Verify all count y-axes use symlog scale (check a few representative SVGs)
- [ ] 10.9 Verify rescue boxplots show both YES and NO data
- [ ] 10.10 Verify all 7 heatmaps have text labels on cells
- [ ] 10.11 Verify FP heatmap includes Somatic classification
- [ ] 10.12 Verify per-sample charts show only within-set samples in each facet
- [ ] 10.13 Verify disease-wise charts use 4×3 grid layout
- [ ] 10.14 Verify set-wise charts use 2×2 grid layout
- [ ] 10.15 Verify BAM DP distribution is clipped at 2000
- [ ] 10.16 Run test suite, confirm no regressions
