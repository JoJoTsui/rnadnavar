## Context

The seq2neo visualizer (`bin/vcf_stats/seq2neo/visualizer.py`, ~2400 lines) generates 100+ SVG charts across 9 wise dimensions. Round 5 introduced a `_count_scale()` helper (line 274) used by 20 chart functions to standardize count axes. It returns `alt.Scale(type="symlog", constant=1)`, which was chosen to handle zero values. However, symlog compresses the visual range so severely that 88% of chart height covers only the near-zero region, making bar comparisons across decades impossible.

Six agents reviewed every SVG and TSV in the 61-sample output and traced each issue to code-level patterns in the visualizer. The 9 user-reported issues collapse into 7 systematic patterns, all fixable at shared code points.

## Goals / Non-Goals

**Goals:**
- Replace symlog with log scale on all variant count axes (20 functions via 1 helper)
- Fix chart 43 data correctness bug (tripled heatmap text from un-aggregated partition data)
- Enable sample-wise charts to facet by set_number
- Swap chart 22 axes (sample→x, count→y) and apply log scale
- Convert threshold sweep charts (31, 35) from vertical stack to 3-column grid
- Remove redundant BAM chart legends
- Increase heatmap cell size / abbreviate large numbers

**Non-Goals:**
- Changing the underlying statistics computation (TSV data is correct)
- Migrating boxplot functions from `alt.Column()` to `.facet()` (boxplot is composite — this is an Altair limitation, not a bug)
- Re-running BAM pileup or VCF parsing

## Decisions

### D1: Log scale with `domainMin=1` instead of pure `log`

**Decision:** `_count_scale()` returns `alt.Scale(type="log", domainMin=1)`.

**Rationale:** Pure `log` maps zero to negative infinity, which would cause Vega-Lite rendering failures for any bar with count=0. Adding `domainMin=1` ensures the axis starts at 1. Count data from `group_by` aggregations omits zero-count groups, so this is a safety net rather than a common case.

**Alternative considered:** `alt.Scale(type="log")` without domainMin — rejected because edge cases with zero counts would produce silent rendering failures.

**Alternative considered:** `log1p` transform (log(1+x)) — rejected because Vega-Lite doesn't support it natively, and it would distort the tick labels.

### D2: Aggregate chart 43 data by summing across `partition`

**Decision:** Add `pdf = pdf.group_by(["classification", "vaf_bin", "dp_bin"]).agg(pl.col("count").sum())` before charting in `plot_filter_vaf_dp_heatmap`.

**Rationale:** The partition column (train/val/test) is a data-split artifact, not a visualization dimension. The chart facets by classification but not partition, causing 3 text labels to overlap at each cell position. Summing across partitions gives the correct total count per cell.

### D3: Join set_number to sample-wise data at the CLI level

**Decision:** In `cli.py`, when generating sample-wise charts, join `set_number` from the sample manifest (already available as `all_stats` DataFrame) to the sample-level aggregation data.

**Rationale:** The root-level charts (08, 22) already have set_number because they use `sample_summary` which includes it. The `sample/` wise dimension charts use a different data path that strips set_number during aggregation. The fix is to re-join it before passing to chart functions.

### D4: Threshold sweep `.facet(columns=3)` replacing `row=alt.Row()`

**Decision:** For `plot_vaf_threshold_sweep`, remove `row=alt.Row("classification:N")` from encode and append `.facet(facet="classification:N", columns=3)`. Keep all 7 classifications (3×3 grid with 2 empty cells). Same for `plot_dp_threshold_sweep`.

**Rationale:** `mark_line` is NOT a composite mark, so `.facet()` works correctly (unlike `mark_boxplot`). Using `columns=3` gives a compact 3×3 grid instead of a 2333px vertical stack. Keeping "Overall" preserves the current information; excluding it would require a data filter and might surprise users.

### D5: Abbreviated number format for dense heatmaps

**Decision:** Use Vega-Lite `format=".2s"` for text labels in `plot_fp_cross_tab_heatmap` (yields "94M" instead of "94,003,313"). Increase cell step to 35px.

**Rationale:** 20px cells cannot display 10-character formatted numbers at any readable font size. Abbreviated SI format (`.2s`) is standard in scientific visualization. Combined with slightly larger cells (35px), text fits without overflow.

### D6: BAM legend suppression via `legend=None`

**Decision:** When color encoding matches x-axis encoding (same field `metric:N`), add `legend=None` to the color specification.

**Rationale:** The legend is 100% redundant with x-axis labels. Removing it saves horizontal space and eliminates the user-reported "duplication" confusion.

## Risks / Trade-offs

**[Log scale hides zero-count bars]** → Mitigated by `domainMin=1`. Zero-count groups are almost never present in aggregated data; if they appear, the bar is simply absent (correct behavior — a bar of height 0 has no visual representation anyway).

**[`.facet(columns=3)` may change axis sharing behavior in threshold sweeps]** → Mitigated by testing. The `.facet()` method with `mark_line` should share scales by default, matching the current `row=` behavior. If not, add `.resolve_scale(y="shared")`.

**[Sample-wise set_number join may increase data size]** → Negligible. set_number is a single integer column added to 61-row DataFrames.

**[Abbreviated number format loses precision]** → Acceptable for heatmap overview. Full values remain available in the interactive tooltip (HTML version) and in the TSV files.
