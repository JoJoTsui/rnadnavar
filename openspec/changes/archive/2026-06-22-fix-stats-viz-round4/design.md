## Context

Round 4 of stats/viz fixes based on 3 rounds of adversarial review (9 agents) against 61-sample pipeline output. The most critical finding is a Strelka field mapping inversion (`AD_ALT←TAR` instead of `AD_ALT←TIR`) that makes every Strelka VAF the complement of the true value. Confirmed in production data: `AD_ALT == TAR` for 100.0% of 980 non-null rows. Raw TAR/TIR/TOR columns are preserved in parquet files, enabling post-load repair without re-parsing VCFs.

## Goals / Non-Goals

**Goals:**
- Fix the P0 Strelka TAR/TIR inversion and enable `--resume` repair from existing parquet
- Fix all 7 visualization rendering bugs (DP overflow, violin drop, faceting, text labels, sweep filter)
- Add comprehensive rescue analytics (7 stats functions + 9 charts) covering set-wise, sample-wise, per-FILTER, per-tier breakdowns
- Add text labels to all heatmaps and pie charts for static-render readability

**Non-Goals:**
- Refactoring the wise chart registry architecture (save for round 5)
- Adding new variant callers or tiering logic
- Implementing the streaming/OOM architecture (separate change)
- Switching from altair to another charting library
- Per-sample filtered mode for distribution charts (architectural change)

## Decisions

### D1: Strelka field mapping fix — two-line fix + post-load repair

**Choice:** Fix the alias mapping at `caller_parser.py:447` and `cli.py:181` from `("AD_REF", "TOR"), ("AD_ALT", "TAR")` to `("AD_REF", "TAR"), ("AD_ALT", "TIR")`. Add a `_repair_strelka_columns()` function called during `--resume` that detects and fixes the inversion from raw columns in parquet.

**Detection logic:** Check if `AD_ALT == TAR` for non-null rows. If 100% match → old (buggy) parquet → repair needed. If `AD_ALT == TIR` → already fixed → skip.

**Repair scope:** After remapping AD_ALT/AD_REF, recompute:
1. `{mod}_{caller}_VAF` = TIR / DP (clipped to [0, 1])
2. `DNA_VAF_mean` / `RNA_VAF_mean` = mean_horizontal of per-caller VAFs
3. `DNA_ALT_DP_mean` / `RNA_ALT_DP_mean` = mean_horizontal of per-caller AD_ALTs
4. `DNA_REF_DP_mean` / `RNA_REF_DP_mean` = mean_horizontal of per-caller AD_REFs

**DP unchanged:** FORMAT/DP is read directly from the Strelka VCF and is correct per spec (tier-1 filtered depth). No recomputation needed.

**Defensive cap:** Add `.clip(0.0, 1.0)` in `compute_vaf_columns()` after VAF = AD_ALT / DP for all callers, not just Strelka.

### D2: DP data clipping — clip data, not just axis domain

**Choice:** Add `.clip(0, 2000)` on pandas DP columns after converting from polars, creating `*_display` columns (same pattern as `VAF_display`). Keep `scale=alt.Scale(domain=[0, 2000])` as belt-and-suspenders. Add outlier count annotation as subtitle.

**Affected functions (5):**
- `plot_dna_vs_rna_dp` — clip both axes
- `plot_dp_boxplot_per_tier` — clip y
- `plot_ref_alt_dp_scatter` — clip both axes
- `plot_dp_distribution` — clip y (dead code, wire into registry)
- `plot_per_tier_dp_boxplot` — clip y

### D3: Sample charts — `.facet(columns=2)` with shared value axis

**Choice:** Replace `enc["row"] = alt.Row("set_number:N")` with `.facet(facet='set_number:N', columns=2)` applied after `.properties()`. Use `resolve_scale()` to share the value axis across facets. Set fixed per-panel dimensions.

| Chart | Value Axis | Shared | Per-Panel Size |
|-------|-----------|--------|----------------|
| 08 per_sample_distribution | x (counts) | `resolve_scale(x='shared')` | height=300 |
| 20 bam_metrics | y (reads) | `resolve_scale(y='shared')` | width=350 |
| 22 sample_tier_dist | x (variants) | `resolve_scale(x='shared')` | height=300 |

### D4: Rescue analytics — statistics + charts following existing patterns

**Choice:** Add 7 statistics functions and 9 chart functions following the existing `compute_*()` → TSV → `plot_*()` pattern. Output to `stats/rescue/` and `plots/rescue/` subdirectories. Wire into cli.py as a new "rescue" section after threshold analysis.

**Statistics functions:**
1. `compute_rescue_breakdown(df, group_col="set_number")` — per-group rescued/non-rescued counts + proportions
2. `compute_rescue_by_filter(df)` — rescued × FILTER cross-tab
3. `compute_rescue_cross_tab(df)` — rescued × FILTER × set_number 3-way
4. `compute_rescue_vaf_dp(df)` — VAF/DP quartiles by rescue status
5. `sample_rescue_summary(df)` — per-sample rescued counts with FILTER breakdown
6. `compute_rescue_by_tier(df)` — rescue rate per CxDy tier
7. `compute_rescue_by_caller_support(df)` — N_SUPPORT_CALLERS × RESCUED distribution

**Chart functions (48-56):**
48. `plot_rescue_breakdown` — stacked bar per set with % labels
49. `plot_rescue_by_filter` — grouped bar: rescued/non-rescued per FILTER
50. `plot_rescue_cross_tab_heatmap` — 3-way heatmap with text labels
51. `plot_rescue_vaf_boxplot` — VAF distributions rescued vs non-rescued
52. `plot_rescue_dp_boxplot` — DP distributions rescued vs non-rescued
53. `plot_rescue_sample_distribution` — per-sample rescue counts, faceted by set (2×2)
54. `plot_rescue_by_tier` — stacked bar: rescue rate per tier
55. `plot_rescue_rate_trend` — line chart: rescue proportion across sets
56. `plot_rescue_caller_support` — grouped bar: N_SUPPORT_CALLERS by rescue status

### D5: Heatmap text labels — `mark_text()` overlay with conditional color

**Choice:** Add a `mark_text(baseline="middle", fontSize=9)` layer to heatmap charts. Use conditional color: white text on dark cells (high count), black text on light cells (low count). Threshold at median count value.

```python
text = alt.Chart(pdf).mark_text(baseline="middle", fontSize=9).encode(
    x=..., y=...,
    text=alt.Text("count:Q", format=",d"),
    color=alt.condition(
        alt.datum.count > median_count,
        alt.value("white"), alt.value("black")
    )
)
chart = rect + text
```

For the faceted `plot_filter_vaf_dp_heatmap`, layer rect+text BEFORE applying `.facet()`.

### D6: Pie chart labels — `mark_text(radius=80)` with pre-computed labels

**Choice:** Pre-compute percentage and label string in pandas, then add `mark_text(radius=80)` layer:

```python
pdf["pct"] = (pdf["n_variants"] / pdf["n_variants"].sum() * 100).round(1)
pdf["label"] = pdf.apply(lambda r: f"{r['n_variants']:,}\n({r['pct']}%)", axis=1)
text = alt.Chart(pdf).mark_text(size=10, radiusOffset=20).encode(
    theta=alt.Theta("n_variants:Q", stack=True),
    text="label:N"
)
chart = arc + text
```

### D7: DP sweep filter — mirror VAF sweep pattern

**Choice:** Add `classification.is_null()` filter in `plot_dp_threshold_sweep` before `to_pandas()`, identical to `plot_vaf_threshold_sweep:1333`.

### D8: VAF violin — restore in faceted branch

**Choice:** In the faceted branch (line 649+), build both `transform_density` violin and boxplot layers, then combine with `alt.layer()` before applying `facet`. Currently only boxplot is rendered.

### D9: Caller concordance — `.facet(columns=3)` wrap

**Choice:** Replace `enc["column"] = alt.Column(f"{color_col}:N")` with `.facet(facet=f"{color_col}:N", columns=3)`. Add `width=200` per panel.

## Risks / Trade-offs

- **[Risk]** Post-load Strelka repair changes parquet data in memory but doesn't overwrite files → next `--resume` re-repairs. Mitigated: repair is idempotent (detect→skip if already fixed).
- **[Risk]** `.clip(0, 2000)` on DP data hides 0.45% of values → mitigated by outlier count annotation.
- **[Risk]** Layered `mark_text` on faceted heatmaps may have Altair rendering order issues → test with both SVG and HTML exports.
- **[Risk]** 9 new rescue charts increase total chart count to ~56 → minor impact on pipeline runtime (~seconds per chart).
- **[Trade-off]** `--resume` repair recomputes VAF means but does NOT recompute tiering (N_*_CALLERS_SUPPORT depends on VAF thresholds). Full accuracy requires a fresh run. Document this limitation.
