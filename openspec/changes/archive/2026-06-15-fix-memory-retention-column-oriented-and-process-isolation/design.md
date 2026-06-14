## Context

The seq2neo stats pipeline processes 65 cancer samples (RNA+DNA) through a rescue VCF parser
(Rust noodles-vcf via pyo3), streaming caller joins, tiering, statistics, and visualization.
Despite prior memory optimizations (streaming join, column pruning, large-sample throttle),
the pipeline is OOM-killed at 8-12 samples on a 200 GB cgroup.

Memory logging (`/proc/self/status VmRSS`) at every step boundary reveals monotonic RSS growth
even after `del df; gc.collect()` at the end of each sample. Controlled simulations confirmed
four compounding causes:

1. **Row-oriented rescue parsing** creates 7M PyDicts × ~80 keys = ~560M Python objects per
   large sample. After `pl.DataFrame(records)`, pymalloc retains the arena pages (~43 GB).
2. **polars Arrow arrays** are the largest single consumer (~55 GB for the joined+stacked DF
   at 7M variants). When freed, glibc malloc retains the pages.
3. **glibc malloc heap fragmentation** prevents `malloc_trim` from lowering the program break
   below fragmented pages. Each allocation-free cycle leaves permanent RSS scars.
4. **Python `.to_list()` calls** materialize 7M strings for derived column computation,
   adding pressure to pymalloc.

Simulations showed RSS baseline climbing 4.2× across 4 sample cycles (0.40 → 1.68 GB after
`malloc_trim`). At production scale (7M variants), this drives RSS from 11 GB to 200+ GB.

## Goals / Non-Goals

**Goals:**
- Eliminate the 7M-PyDict bottleneck in rescue VCF parsing via column-oriented Rust return
- Guarantee memory reclamation after each sample via process isolation
- Replace Python-loop derived columns with polars vectorized expressions
- Maintain identical output: same parquet files, same stats CSVs, same charts
- No performance regression; fixes should be faster or neutral

**Non-Goals:**
- Changing the statistics computation or visualization logic
- Modifying the BAM pileup or BAM stats (separate concern)
- Optimizing the Rust VCF record parsing itself (noodles-vcf is already fast)
- Removing the large-sample throttle (still needed even with process isolation)

## Decisions

### Decision 1: Column-oriented rescue parsing in Rust

**Choice:** New `parse_rescue_columns()` pyfunction returns `{column_name: PyList}` — same
pattern already proven by `parse_caller_vcf()`.

**Rust side** (`vcf.rs` + `lib.rs`):
- New struct `RescueColumns` with `Vec` per fixed column + `HashMap<String, Vec<Option<String>>>`
  for INFO columns
- Parse loop populates column vectors directly instead of building `Vec<RescueRecord>`
- `variant_type` and `ti_tv` computed in Rust during the parse pass (no extra iteration)
- INFO fields: iterate `info_keys` for every record, push `Some(value)` or `None`
- GIL released during parse (`py.detach()`); Python object construction re-acquires GIL

**Python side** (`rust_vcf.py`):
- Replace `pl.DataFrame(records)` + `_cast_columns()` + `.to_list()` with:
  ```python
  cols = stats_core.parse_rescue_columns(path)
  series = [pl.Series(name, vals, dtype=...) for name, vals in cols.items()]
  df = pl.DataFrame(series)
  ```
- No `_cast_columns` needed — Rust returns correct types (Int64/Float64/Bool/String)
- No `.to_list()` needed — `variant_type` and `ti_tv` are already columns from Rust
- Remove `_cast_columns` function entirely

**Alternatives considered:**
- Arrow IPC from Rust: Eliminates Python lists entirely but adds `pyo3-arrow` dependency,
  complex build setup. Rejected — Python lists → polars Series is fast enough.
- Keep row-oriented but add `malloc_trim` only: Tested — only 30% memory recovered per cycle.
  Insufficient for 65 samples.
- polars-native derived columns in Python: Works for variant_type/ti_tv but doesn't solve
  the PyDict problem. Computing in Rust during parse is cleaner.

### Decision 2: Process isolation via multiprocessing.Pool

**Choice:** `multiprocessing.Pool(processes=N, maxtasksperchild=1, context='spawn')` with
`pool.imap_unordered()`.

`spawn` context ensures each worker is a **fresh Python interpreter** — no memory inheritance
from parent (unlike `fork` which COW-shares parent pages, leading to the same retention issue).

`maxtasksperchild=1` ensures each worker processes **exactly one sample** then exits. The
kernel reclaims ALL memory: pymalloc arenas, glibc heap, polars pool, Rust allocator.

**Cross-process large-sample throttle:**
- Current: `threading.Semaphore(1)` — doesn't work across processes
- New: `multiprocessing.Manager().Semaphore(1)` — shared via manager process
- Worker function acquires semaphore before processing samples >2M variants
- Manager process overhead: ~0.1ms per acquire/release, negligible

**Worker function signature:**
```python
def _process_one_worker(args):
    """Top-level function for multiprocessing (must be picklable)."""
    row, max_workers, use_rust, output_dir, large_sem = args
    # ... process_single_sample + write_parquet + malloc_trim ...
    return {"sample_id": sid, "stats": stats}  # small dict, no DF
```

**Speed impact:** `spawn` + module imports = 2-4s per worker creation. With 4 workers and
65 samples, ~17 process creations, each overlapping with other workers' processing.
Net wall-time overhead: ~35-70 seconds. Acceptable given total runtime (minutes to hours
depending on BAM pileup enabled).

**Alternatives considered:**
- `fork` context with `maxtasksperchild=1`: Fork is faster but COW means child inherits
  parent's dirty pages. CPython reference counting dirties pages on every object access.
  Rejected — doesn't guarantee clean memory slate.
- `subprocess.run()` per sample: Cleanest isolation but slowest (Python startup per sample
  with no worker reuse). Rejected — more overhead, same result.
- Keep ThreadPoolExecutor + aggressive `malloc_trim`: Tested — residual creep 0.40→1.68 GB
  across 4 samples (4.2×). At production scale (65 samples, 7M variants), residual would
  still hit 60-100 GB. Rejected — insufficient guarantee.

### Decision 3: polars-native derived columns

**Choice:** Use polars `when/then/otherwise` chain for `variant_type` and string operations
for `ti_tv`.

```python
df = df.with_columns([
    pl.when(pl.col("ALT").str.len_chars() > pl.col("REF").str.len_chars())
    .then(pl.lit("INS"))
    .when(pl.col("REF").str.len_chars() > pl.col("ALT").str.len_chars())
    .then(pl.lit("DEL"))
    .when((pl.col("REF").str.len_chars() == 1) & (pl.col("ALT").str.len_chars() == 1))
    .then(pl.lit("SNV"))
    .otherwise(pl.lit("MNV"))
    .alias("variant_type"),
])
```

These are fallback implementations — the Rust column-oriented parser computes both columns
natively. The polars expressions serve as the Python fallback path (cyvcf2 parser).

**Benchmark:** polars-native is 6× faster than `.to_list()` + Python list comprehension
(20ms vs 120ms for 500K rows).

### Decision 4: malloc_trim after each sample

**Choice:** Call `ctypes.CDLL("libc.so.6").malloc_trim(0)` after `write_parquet + del df + gc.collect()`.

With process isolation this is a belt-and-suspenders measure — the process exits immediately
after, so it has no practical effect on RSS. It's included because:
- It costs <1ms
- It's useful if a user runs with `--sample-workers 1` (no multiprocessing)
- It helps during the transition period (Fix A deployed before Fix B)

## Risks / Trade-offs

- **[Risk] `multiprocessing.Manager().Semaphore` creates a manager process that could crash**
  → **Mitigation:** The manager process is lightweight (only manages the semaphore). If it
  crashes, the pool workers will hang on `acquire()`. Add a timeout: `sem.acquire(timeout=300)`.
  If timeout expires, process the sample anyway (worse case: two large samples run concurrently,
  still within 200 GB after Fix A).

- **[Risk] `spawn` context means all worker args must be picklable**
  → **Mitigation:** All args are simple types (dict, str, int, bool). The `row` dict from
  `manifest.to_dicts()` is picklable. The `large_sem` from Manager is a proxy object designed
  for cross-process use.

- **[Risk] `maxtasksperchild=1` means no worker reuse — import overhead per sample**
  → **Mitigation:** Acceptable. 2-4s overhead per sample vs minutes of processing. For
  very small samples (<100K variants, <10s processing), the overhead dominates. But the
  pipeline is designed for large samples where minutes of processing dwarf 2-4s startup.

- **[Risk] Rust column-oriented parser might produce different output than row-oriented**
  → **Mitigation:** The cyvcf2 fallback path remains unchanged. Tests verify column-oriented
  Rust output matches row-oriented Rust output (same values, same types). During transition,
  keep `parse_rescue` as a wrapper that validates against `parse_rescue_columns`.

- **[Trade-off] Process isolation prevents memory reuse across samples**
  Each process starts fresh — no shared polars string cache, no shared Arrow buffers.
  This is intentional (the whole point is to prevent cross-sample memory accumulation).
  The per-sample memory is well within limits after Fix A reduces peak by ~40 GB.

## Migration Plan

1. Deploy column-oriented rescue parsing (Fix A + C) first — backward compatible, faster
2. Deploy process isolation (Fix B + D) — requires CLI flag changes
3. Keep `--sample-workers` flag; internally switches from ThreadPoolExecutor to multiprocessing.Pool when >1
4. Add `--process-mode` flag (`thread`/`spawn`) for gradual rollout
5. No data format changes — parquet files, CSVs, and charts are identical
6. Rollback: revert to `--process-mode thread` if spawn overhead is problematic

## Open Questions

- Should `maxtasksperchild` be configurable? Default 1 (max safety), but users may want
  2-4 for faster small-sample processing. → **Decision: make it configurable via
  `--max-tasks-per-child` flag, default 1.**
