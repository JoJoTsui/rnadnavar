## Why

The 24-sample pipeline is OOM-killed after BAM stats during cross-sample aggregation. The root cause: `combined_df = pl.scan_parquet("*.parquet")` is `.collect()`'d **8 times without column projection** by the aggregation functions (`disease_summary`, `tier_summary`, `dataset_summary`, etc.). Each `.collect()` loads all 58M variants × 165 columns (~80 GB each). Even though each DataFrame is freed between calls, glibc retains every freed page — RSS climbs 60→80→100→120→...→200 GB until the cgroup OOM-killer triggers. This breaks at any scale: with 3× the data (175M variants), even a single full `.collect()` would allocate 200+ GB alone.

Additionally, the `fork` context for `multiprocessing.Pool` deadlocks due to polars' rayon thread pool being killed by `fork()` but leaving internal locks held. Revert to `spawn` and suppress the cosmetic `resource_tracker` warning.

## What Changes

- **Streaming aggregation**: Replace `_ensure_eager(df)` (which does `df.collect()` loading all columns and rows) with column-projected lazy aggregation: `df.select(needed_cols).group_by(...).agg(...).collect()`. Each aggregation function only reads the 5-10 columns it needs and returns a small aggregated result (tens of rows), not a materialized 58M-row DataFrame.
- **Back to spawn + warning suppression**: Revert from `fork` to `spawn` context for `multiprocessing.Pool`. Add `warnings.filterwarnings('ignore', message='resource_tracker')` — the 6 leaked semaphores are POSIX named semaphores cleaned up by the kernel on process exit; the warning is cosmetic.
- **Reduce test time**: Use `pytest.fixture(scope="class")` in `TestColumnOrientedRescueParser` to parse the rescue VCF once and share the result across all 7 tests, reducing 5 parses to 1 (~80s savings).
- **Remove duplicate code**: Delete the duplicate `_maybe_collect` definition in `visualizer.py` (lines 23-27 are identical to lines 16-20).

## Capabilities

### New Capabilities
- `streaming-cross-sample-aggregation`: Cross-sample aggregation functions use polars lazy queries with column projection and streaming group_by, loading only needed columns and returning small aggregated DataFrames instead of materializing all 58M+ rows.

### Modified Capabilities
- `process-isolated-sample-processing`: Revert from `fork` to `spawn` context with warning suppression (no behavior change, same per-process isolation).

## Impact

- **Python**: `statistics.py` (replace `_ensure_eager` with column-projected lazy aggregation in all cross-sample functions), `tiering_stats.py` (same), `cli.py` (revert to spawn, suppress warning, remove `_ensure_eager` calls, add memory logging), `visualizer.py` (remove duplicate `_maybe_collect`)
- **Tests**: `test_seq2neo_stats.py` — add class-scoped fixture to `TestColumnOrientedRescueParser`
- **Memory at 3× scale**: Peak RSS ~20 GB (streaming buffers) vs current ~200+ GB (repeated full collects). Scales linearly with file count (I/O) but constant with row count (memory).
