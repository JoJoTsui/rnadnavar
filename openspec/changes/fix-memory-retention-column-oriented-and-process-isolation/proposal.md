## Why

The 65-sample stats pipeline is consistently OOM-killed at 8-12 samples despite 200 GB cgroup
limits and prior fixes (streaming join, column pruning, large-sample throttling). Analysis with
`/proc/self/status` memory logging and controlled simulations reveals the root cause:
**allocator retention across samples**. Both pymalloc (CPython) and glibc malloc retain freed
memory pages indefinitely. The rescue VCF parser creates 7M PyDicts per large sample (~43 GB of
Python objects), and after they're freed the pages are never returned to the OS. As increasingly
large samples are processed, RSS grows monotonically from 11 GB → 33 GB → 45 GB → 124 GB → 204 GB
until the cgroup OOM-killer triggers.

## What Changes

- **Column-oriented rescue VCF parsing (Rust)**: Rewrite `stats_core::parse_rescue` to return
  `{column_name: [values]}` instead of row-oriented `[{field: value}]`. Eliminates 7M PyDicts
  per large sample — the single largest source of pymalloc retention (~43 GB).
- **Process isolation via multiprocessing**: Replace `ThreadPoolExecutor` with
  `multiprocessing.Pool(maxtasksperchild=1, context='spawn')`. Each sample runs in a fresh Python
  process; when the process exits the kernel reclaims ALL memory (pymalloc, glibc, polars pool,
  Rust allocator). Zero cross-sample memory creep.
- **polars-native derived columns**: Replace `.to_list()` + Python list comprehension for
  `variant_type`/`ti_tv` with polars vectorized `when/then/otherwise` expressions. Eliminates
  materializing 7M Python strings and is 6× faster.
- **`malloc_trim(0)` after each sample**: Force glibc to release free pages back to the kernel
  via `madvise(MADV_DONTNEED)`. Modest additional benefit on top of process isolation.
- **Cross-process large-sample throttle**: Replace `threading.Semaphore` with
  `multiprocessing.Manager().Semaphore` for the >2M-variant exclusive-access gate.

## Capabilities

### New Capabilities
- `column-oriented-rescue-parsing`: Rust rescue VCF parser returns column-oriented data directly,
  bypassing row-oriented PyDict construction and eliminating per-field type casting.
- `process-isolated-sample-processing`: Each sample is processed in an independent Python process
  that exits after completion, guaranteeing all memory is reclaimed by the kernel.

### Modified Capabilities
- None — output behavior is identical. Same parquet files, same stats, same charts.

## Impact

- **Rust**: `vcf.rs` (new `RescueColumns` struct + column-oriented parser), `lib.rs`
  (new `parse_rescue_columns` pyfunction, keep `parse_rescue` for backward compat during
  transition)
- **Python**: `rust_vcf.py` (Series-based DF construction, remove `_cast_columns`, remove
  `.to_list()`), `cli.py` (multiprocessing.Pool, Manager().Semaphore, malloc_trim, polars-native
  derived columns)
- **Speed**: Column-oriented parsing is 6-10× faster for Python post-processing. polars-native
  derived columns are 6× faster than current `.to_list()` loop. Process spawn overhead is ~2-4s
  per worker, amortized by parallelism (~1 min total wall-time overhead for 65 samples).
- **Tests**: `tests/test_seq2neo_stats.py` — add memory isolation tests, column-oriented format
  tests, update existing tests for new parse function signature.
