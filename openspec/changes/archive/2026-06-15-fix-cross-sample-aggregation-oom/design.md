## Context

The cross-sample aggregation phase runs after per-sample variant processing and BAM statistics. It calls `pl.scan_parquet("*.parquet")` to create a lazy frame over all per-sample parquet files, then passes it to 8+ aggregation functions. Each function calls `_ensure_eager(df)` which does `df.collect()` — loading ALL 165 columns from ALL 58M+ rows without column projection. Repeated `.collect()` calls cause glibc heap fragmentation and monotonic RSS growth until the cgroup OOM-killer triggers.

Additionally, the `fork` context for `multiprocessing.Pool` deadlocked because polars' rayon thread pool is killed by `fork()` but leaves internal locks held. `spawn` avoids this (fresh process via fork+exec) but produces a cosmetic `resource_tracker` warning about 6 leaked POSIX semaphores — these are cleaned up by the kernel on process exit.

## Goals / Non-Goals

**Goals:**
- Reduce peak RSS during cross-sample aggregation from ~200 GB to ~20 GB (at current 58M-variant scale)
- Scale to 175M+ variants (3×) without exceeding 50 GB peak RSS
- Maintain identical CSV output (same aggregation results)
- Eliminate the `fork` deadlock by reverting to `spawn`
- Reduce test runtime by ~80s

**Non-Goals:**
- Changing visualization behavior (already uses column-projected lazy queries)
- Changing per-sample processing or BAM stats
- Implementing true streaming in visualization (already efficient)

## Decisions

### Decision 1: Column-projected lazy aggregation (streaming)

**Choice:** Replace `_ensure_eager(df)` with column-projected lazy queries in every cross-sample aggregation function.

Each function declares the columns it needs, pushes selection into the lazy scan, groups/aggregates, then collects only the aggregated result:

```python
# BEFORE (current — OOM at scale):
def disease_summary(df):
    df = _ensure_eager(df)  # .collect() → 58M rows × 165 cols → ~80 GB
    return df.group_by("disease").agg([...])  # returns 20 rows

# AFTER (streaming — scales to any size):
def disease_summary(df):
    needed = ["disease", "disease_normalized", "final_tier", ...]
    return (
        df.select(needed)           # push column selection into parquet scan
          .group_by(["disease", "disease_normalized"])
          .agg([...])
          .collect()                # collect AFTER aggregation → 20 rows
    )
```

polars' query optimizer pushes the `.select()` into the parquet reader, so only the needed columns are read from disk. The `.group_by().agg()` runs in streaming mode — rows are aggregated as they're read, not accumulated. Memory usage is proportional to the number of unique group keys (tens to hundreds), not the number of input rows (58M+).

**Column requirements per function (determined by reading current implementations):**

| Function | Columns needed |
|----------|---------------|
| `disease_summary` | disease, disease_normalized, DNA_VAF_mean, RNA_VAF_mean, DNA_DP_mean, RNA_DP_mean, DNA_REF_DP_mean, DNA_ALT_DP_mean, RNA_REF_DP_mean, RNA_ALT_DP_mean, final_tier, caller_tier, database_tier, tier_quality, variant_type, ti_tv |
| `tier_summary` (tiering_stats) | final_tier, DNA_VAF_mean, RNA_VAF_mean, DNA_DP_mean, RNA_DP_mean, DNA_REF_DP_mean, DNA_ALT_DP_mean, RNA_REF_DP_mean, RNA_ALT_DP_mean, variant_type, ti_tv |
| `dataset_summary` | sample_id, disease, disease_normalized, set_number, final_tier, caller_tier, database_tier, DNA_VAF_mean, RNA_VAF_mean, variant_type, ti_tv |
| `sample_tier_summary` | sample_id, set_number, disease, final_tier, DNA_VAF_mean, RNA_VAF_mean, DNA_DP_mean, RNA_DP_mean |
| `caller_overlap_distribution` | final_tier |
| `filter_distribution` | FILTER, final_tier, set_number |
| `variant_type_distribution` | variant_type, final_tier, set_number |
| `gt_concordance` | GT columns (6 caller GT cols) |

**Total unique columns needed across all functions: ~25 columns** (out of 165).

**Alternatives considered:**
- **Materialize once with column pruning, reuse DataFrame**: Works at current scale (~15 GB for 25 cols × 58M rows) but breaks at 3× scale (~160 GB). Rejected — doesn't scale.
- **Use polars' new streaming engine (`collect(engine='streaming')`)**: Would work but requires polars 1.42+. Rejected — not available in 1.41.2.
- **Process per-sample aggregates, merge**: More code, more I/O. Rejected — lazy aggregation in polars already streams.

### Decision 2: Revert fork → spawn with warning suppression

**Choice:** Change `mp.get_context("fork")` back to `mp.get_context("spawn")`. Add `warnings.filterwarnings('ignore', message='resource_tracker')` to suppress the cosmetic warning.

**Rationale:**
- `fork` deadlocks with polars' rayon thread pool (threads killed at fork, locks held)
- `spawn` creates a fresh process via fork+exec, no inherited lock state
- The 6 leaked semaphore warning is cosmetic — POSIX named semaphores (`/dev/shm/sem.*`) are cleaned up by the kernel on process exit
- `spawn` worker startup is ~3s (module re-import) vs fork's ~1ms, but with `maxtasksperchild=1`, this overhead is paid once per worker per round (max 6 rounds for 24 samples → ~18s total, negligible vs hours of processing)

### Decision 3: Class-scoped pytest fixture for rescue VCF

**Choice:** Add a `scope="class"` fixture to `TestColumnOrientedRescueParser` that parses the rescue VCF once and caches both column-oriented and row-oriented results.

```python
@pytest.fixture(scope="class")
def rescue_data(self, request):
    # ... parse once ...
    return cols, rows
```

The 7 tests currently parse the VCF independently (5 parses × ~19s = ~95s overhead). With a class-scoped fixture, only 1 parse (~19s), saving ~80s.

### Decision 4: Remove duplicate `_maybe_collect`

The function is defined twice identically (lines 16-20 and 23-27) in `visualizer.py`. Delete the second copy.

## Risks / Trade-offs

- **[Risk] Some aggregation functions need columns from caller data that aren't in the ~25 shared columns** → **Mitigation:** Audit each function's column usage before implementing. Add any missing columns to the selection list.
- **[Risk] `spawn` overhead grows with sample count** → **Mitigation:** The 3s per-worker startup is amortized across minutes of sample processing. At 24 samples with 4 workers, that's ~18s total overhead. Acceptable.
- **[Trade-off] Lazy group_by with many unique keys** → polars streaming group_by buffers one row per unique key. For low-cardinality keys (disease: ~10, tier: ~100, sample_id: ~100), this is negligible. For high-cardinality group_by, this could be problematic — but none of our cross-sample aggregations use high-cardinality keys.

## Memory Budget (58M variants, 24 samples)

```
Step                                    Memory        Result rows
──────────────────────────────────────────────────────────────────
disease_summary                         ~1 GB  →      20 rows
tier_summary                            ~1 GB  →      100 rows
dataset_summary                         ~1 GB  →      1 row (dict)
sample_tier_summary                     ~1 GB  →      2,000 rows
caller_overlap_distribution             ~1 GB  →      50 rows
filter_distribution                     ~1 GB  →      20 rows
variant_type_distribution               ~1 GB  →      4 rows
gt_concordance                          ~1 GB  →      10 rows
validation (read 1 parquet at a time)   ~15 GB →      0 rows (just reports)
visualizations (lazy, column-projected) ~2 GB  →      pre-aggregated charts
──────────────────────────────────────────────────────────────────
Peak RSS:                              ~20 GB
At 3× scale (175M variants):           ~25 GB  (same streaming buffers)
```

Memory is now bounded by streaming buffer size (~1-2 GB per aggregation, freed between calls), not by input data size.
