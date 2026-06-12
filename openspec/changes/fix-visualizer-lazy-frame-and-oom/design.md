## Context

`visualizer.py` was originally written for eager `pl.DataFrame`. When `combined_df` was changed to `pl.scan_parquet()` (a `LazyFrame`) in the streaming architecture, none of the chart functions were updated. The OOM kills during cross-sample aggregation masked these bugs — the pipeline never reached the visualization phase. Now that aggregation is fixed, the first chart call crashes.

## Goals / Non-Goals

**Goals:**
- Fix all 17 LazyFrame crash bugs
- Ensure no chart loads >50K rows into memory without explicit sampling
- Standardize data loading pattern across all 23 chart functions
- Remove dead code (duplicate `_maybe_collect`, broken `BAM_DP_*` chart)

**Non-Goals:**
- Changing chart layout, colors, or visual design
- Changing the dashboard HTML generation
- Adding new chart types
- Parallelizing chart generation (separate optimization)

## Decisions

### Decision 1: Two shared helpers instead of `_maybe_collect`

Replace the single `_maybe_collect()` with two focused helpers:

```python
def _count_rows(df):
    """Efficient row count on LazyFrame — reads no data columns."""
    if isinstance(df, pl.LazyFrame):
        return df.select(pl.len()).collect().item()
    return df.height

def _sample_if_large(df, max_rows=5000):
    """Sample LazyFrame if it exceeds max_rows, then collect."""
    n = _count_rows(df)
    if n > max_rows:
        df = df.sample(max_rows)
    return df.collect()
```

Every chart that needs sampling uses `_sample_if_large`. Every chart that only needs counts uses `_count_rows`. No more `.height` on lazy frames — it doesn't exist.

### Decision 2: Fix pattern per bug category

**Category A — height-before-sampling (12 functions):**
```python
# BEFORE (broken):
pdf = df.select(cols).drop_nulls()
if pdf.height > 5000:           # ← CRASH
    pdf = pdf.sample(5000)
pdf_pd = pdf.pipe(_maybe_collect).to_pandas()

# AFTER (fixed):
pdf = df.select(cols).drop_nulls()
pdf = _sample_if_large(pdf, max_rows=5000)
pdf_pd = pdf.to_pandas()
```

**Category B — count queries (cosmic_gnomad, caller_agreement):**
```python
# BEFORE (broken):
n = df.height                              # ← CRASH
has_cosmic = df.filter(...).height         # ← CRASH

# AFTER (fixed):
n = _count_rows(df)
has_cosmic = _count_rows(df.filter(pl.col("COSMIC_ID").is_not_null()))
```

**Category C — iter_rows (gt_concordance × 2):**
```python
# BEFORE (broken):
for row in df.select(cols).iter_rows():     # ← CRASH

# AFTER (fixed):
pdf = df.select(cols).collect()
for row in pdf.iter_rows():
```

### Decision 3: Remove `plot_bam_coverage_violin` or fix data source

`plot_bam_coverage_violin` accesses `BAM_DP_*` columns from `combined_df`. These columns exist only in `bam_validation.csv` (written by `validate_bam_one`), NOT in the per-sample parquet files. The function silently returns today. Fix: remove the chart from the pipeline (it's the only chart that references non-existent columns).

### Decision 4: Reorganize function layout

Current layout is ad-hoc. Reorganize into clear sections:
1. Helpers (`_count_rows`, `_sample_if_large`, `_save_chart`)
2. Aggregate charts (group_by → small result, no sampling needed)
3. Sampled charts (use `_sample_if_large`)
4. Count-based charts (use `_count_rows` only)
5. Per-sample/small-data charts (use eager DataFrames, no lazy concerns)
6. Dashboard assembly

This makes it clear which charts have lazy frame concerns and which don't.

## Risks / Trade-offs

- **[Risk] `_count_rows` on 58M rows with filters is a full scan** → For `plot_caller_agreement_matrix` with 36 pairwise counts, that's 36 sequential scans of 24 parquet files. Slow (potentially minutes). → **Mitigation:** Acceptable — this chart was already broken, so any working version is an improvement. Can be optimized later with single-pass computation.
- **[Trade-off] `iter_rows()` replacement loads all GT data** → `plot_gt_concordance` loads 58M × 4 string columns (~2 GB) into memory for Python iteration. → **Acceptable:** GT columns are short strings, 2 GB fits in 200 GB limit. Alternative (vectorized computation) would be faster but more complex.
