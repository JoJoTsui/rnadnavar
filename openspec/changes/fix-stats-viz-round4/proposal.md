## Why

Three rounds of adversarial review (9 agents) against 61-sample pipeline output uncovered 1 critical data-correctness bug (Strelka TAR/TIR field mapping inverted — every Strelka VAF is the complement of the true value), 7 visualization rendering bugs (DP overflow, dropped violin, wrong faceting, missing text labels, noisy sweep lines), and 1 major feature gap (no systematic rescue analytics). The Strelka bug affects 2 of 6 callers and corrupts all cross-caller mean VAFs, threshold sweeps, and tiering decisions. The raw TAR/TIR/TOR columns are preserved in parquet files, enabling a post-load repair path compatible with `--resume`.

## What Changes

### Bug Fixes (P0)
- **Strelka TAR/TIR field inversion**: Fix `AD_ALT←TIR` (was TAR), `AD_REF←TAR` (was TOR) in `caller_parser.py:447` and `cli.py:181`. Add post-load repair step for `--resume` to recompute Strelka VAF/means from raw TAR/TIR columns already in parquet. Add defensive `.clip(0, 1)` in `compute_vaf_columns()`.

### Bug Fixes (P1)
- **DP data overflow**: Replace `scale=alt.Scale(domain=[0, 2000])` with actual data `.clip(0, 2000)` in 5 DP chart functions — domain only sets axis viewport, data >2000 still renders and overflows.
- **DP threshold sweep anomaly**: Add `classification.is_null()` filter in `plot_dp_threshold_sweep` to match `plot_vaf_threshold_sweep` — per-classification rows are currently mixed with overall rows.
- **VAF violin dropped in faceted mode**: Restore violin overlay in the faceted branch of `plot_vaf_distribution` (lines 649-661 only render boxplot).
- **Sample charts vertical stack**: Replace `alt.Row("set_number:N")` with `.facet(columns=2).resolve_scale()` in charts 08, 20, 22 for 2×2 set grid with shared value axis.
- **Caller concordance too wide**: Add `columns=3` wrap to `plot_caller_concordance_vs_vaf` faceted mode.
- **Heatmap/pie missing labels**: Add `mark_text()` overlay to charts 43, 45 (count values on heatmap cells) and chart 46 (count + percentage on pie slices).

### New Features (P1)
- **Rescue analytics system**: 7 new statistics functions + 9 new charts for systematic rescued vs non-rescued analysis — per-set breakdowns, rescued×FILTER cross-tabs, VAF/DP distributions, per-sample rescue counts, rescue by tier, rescue by caller support, rescue rate trends.

### Cleanup
- **Dead code**: Wire `plot_dp_distribution` into chart registry or remove it.
- **Test gaps**: Add Strelka mapping correctness test.

## Capabilities

### New Capabilities
- `rescue-analytics`: Systematic statistics and visualizations for rescued vs non-rescued variants — per-set, per-sample, per-FILTER, per-tier breakdowns with cross-tabulations, VAF/DP distributions, and rescue rate trends.
- `strelka-repair`: Post-load repair mechanism for `--resume` to fix Strelka TAR/TIR inversion from raw columns preserved in parquet files without re-parsing VCFs.

### Modified Capabilities
- `variant-visualization`: Fix DP data clipping (5 functions), restore VAF violin in faceted mode, 2×2 set grid for sample charts, caller concordance column wrap, heatmap/pie text labels, DP sweep classification filter.
- `depth-threshold-sweep`: Fix classification filter bug (per-classification rows mixed with overall rows in plot).

## Impact

- `bin/vcf_stats/seq2neo/caller_parser.py` — Strelka field mapping fix (1 line)
- `bin/vcf_stats/seq2neo/cli.py` — Strelka field mapping fix (1 line) + post-load repair step + rescue chart wiring (~100 lines)
- `bin/vcf_stats/seq2neo/statistics.py` — Defensive VAF clip, docstring fix, 7 new rescue statistics functions (~200 lines)
- `bin/vcf_stats/seq2neo/visualizer.py` — DP clipping (5 functions), VAF violin restore, sample chart faceting (3 functions), concordance wrap, heatmap text (2 functions), pie labels, DP sweep filter, 9 new rescue charts (~400 lines)
- `bin/vcf_stats/tests/test_seq2neo_stats.py` — Strelka mapping test, rescue stats tests (~50 lines)
