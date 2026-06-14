## Context

After the first full pipeline run completed (65 samples), visualization generation crashed at `plot_bam_coverage_violin` due to missing columns from one broken sample. Root cause analysis confirmed the broken sample had truncated BAM files — not a code bug. The remaining issues are code-quality and robustness problems discovered during systematic code review.

## Goals / Non-Goals

**Goals:**
- Fix `_save_chart` PNG/SVG crash when vl-convert missing
- Fix `plot_gt_concordance_per_tier` row slicing bug
- Optimize count-based charts (eliminate redundant lazy scans)
- Add `--exclude-sample-ids` CLI flag
- Harden helpers against `ColumnNotFoundError`
- Remove dead code

**Non-Goals:**
- Changing chart layouts, colors, or visual design
- Changing dashboard HTML generation
- Adding new chart types
- Fixing BAM pileup panic (addressed in `fix-bam-pileup-performance`)

## Decisions

### Decision 1: Fix `_save_chart` — HTML always, PNG/SVG best-effort

```python
def _save_chart(chart, name, output_dir):
    # HTML always works (no external deps)
    chart.save(str(html_path))

    # PNG/SVG require vl-convert — optional
    try:
        chart.save(str(png_path), format="png", scale_factor=3)
    except (ImportError, ModuleNotFoundError):
        pass  # vl-convert not installed, skip PNG
    try:
        chart.save(str(svg_path), format="svg")
    except (ImportError, ModuleNotFoundError):
        pass  # vl-convert not installed, skip SVG
```

### Decision 2: Fix `plot_gt_concordance_per_tier` row slicing

The bug: when `facet_col` is present, `cols_to_collect` is `[GT1, GT2, GT3, GT4, caller_tier, facet_val]` and `row[:-1]` returns `[GT1, GT2, GT3, GT4, caller_tier]` — including `caller_tier` as a GT value.

The fix: use explicit indexing on the GT columns only.

```python
# BEFORE (broken):
gts = [g for g in row[:-1] if g is not None and g not in ("./.", "./.", ".")]

# AFTER (fixed):
n_gt = len(existing_gt)
gts = [g for g in row[:n_gt] if g is not None and g not in ("./.", "./.", ".")]
tier = row[n_gt]  # caller_tier is at index n_gt
```

### Decision 3: Batch count queries in slow charts

Instead of per-group lazy scans, pre-aggregate with a single `group_by`:

```python
# BEFORE (slow — O(n_groups) lazy scans):
for group in groups:
    n = _count_rows(df.filter(pl.col(group_col) == group))
    has_cosmic = _count_rows(df.filter(
        (pl.col(group_col) == group) & pl.col("COSMIC_ID").is_not_null()
    ))

# AFTER (fast — single lazy scan):
counts = df.group_by(group_col).agg([
    pl.len().alias("n_total"),
    pl.col("COSMIC_ID").is_not_null().sum().alias("n_cosmic"),
    pl.col("GNOMAD_AF").is_not_null().sum().alias("n_gnomad"),
]).collect()
```

This eliminates redundant parquet scans and makes chart generation O(1) instead of O(n_groups).

### Decision 4: Harden `_sample_if_large` and `_maybe_collect`

```python
def _sample_if_large(df, max_rows=5000):
    try:
        if isinstance(df, pl.LazyFrame):
            df = df.collect()
    except pl.exceptions.ColumnNotFoundError:
        return pl.DataFrame()  # graceful degradation
    if df.height > max_rows:
        return df.sample(max_rows)
    return df
```

### Decision 5: Add `--exclude-sample-ids` CLI flag

```python
parser.add_argument("--exclude-sample-ids", nargs="*", default=None,
                    help="Exclude specific sample IDs from processing")

# After existing filters:
if args.exclude_sample_ids:
    manifest = manifest.filter(~pl.col("sample_id").is_in(args.exclude_sample_ids))
```

### Decision 6: Remove dead code

- `_chromosome_sort_key` (line 77-83): unused function — `plot_chromosome_density` uses inline `when/then` instead
- `plot_per_sample_violin` alias: backward compatibility alias, never called

## Risks / Trade-offs

- **[Risk] Batch count queries change output format** → `plot_cosmic_gnomad_annotation` and `plot_database_enrichment_by_tier` currently build results row-by-row. Switching to `group_by().agg()` changes intermediate DataFrames. → **Mitigation:** Keep the existing per-group iteration for correctness, just replace `_count_rows` calls with indexed lookups from a pre-computed DataFrame.
- **[Trade-off] `_sample_if_large` returns empty DataFrame on error** → Charts silently produce nothing instead of crashing. → **Acceptable:** An empty chart is better than a crashed pipeline. The underlying issue (missing data) is logged elsewhere.
