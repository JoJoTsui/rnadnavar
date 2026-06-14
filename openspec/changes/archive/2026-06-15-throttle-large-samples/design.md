## Context

Samples vary from 250K to 7M variants. All overheads scale linearly with variant count. A 7M-variant sample creates 28M Python strings from `to_list()`, 7M Python tuples from `set(zip())`, and a 2 GB Rust HashSet clone. With `--sample-workers 4`, four large samples compound these overheads beyond cgroup limits.

## All Overheads Found (per 7M-variant sample)

| # | Location | Overhead | Size |
|---|---|---|---|
| 1 | cli.py:95-98 | `to_list()` ×4 cols | 1.25 GB Python objects |
| 2 | cli.py:99 | `set(zip(...))` | 2 GB Python tuples |
| 3 | cli.py:86 | `rescue_df` 191 cols | 10 GB (unavoidable) |
| 4 | caller_parser:350 | `caller_data` dict (6 callers) | 4 GB Python lists |
| 5 | caller.rs:155 | `target_positions.clone()` | 2 GB Rust HashSet |
| 6 | caller.rs:183 | `chrom.clone()` per match | ~50 MB strings |
| 7 | cli.py:205 | `gc.collect()` only after write | Slow GC, memory retained |
| 8 | cli.py:86 | rescue_df has unused columns | ~5 GB could be dropped |

## Goals / Non-Goals

**Goals:**
- Auto-throttle: only 1 large (>2M) sample at a time
- Eliminate `to_list()` + `set(zip())` Python overhead (items 1, 2)
- Avoid Rust HashSet clone (item 5)
- Add gc.collect() at key free points (item 7)

**Non-Goals:**
- Reducing rescue_df column count (item 3, 8 — separate optimization)
- Avoiding caller_data dict entirely (already column-oriented, item 4 acceptable with throttling)

## Decisions

### 1. Auto-throttle with Semaphore after rescue parse

**Decision:** Parse rescue VCF first (fast, memory-efficient in Rust). Then check variant count. If >2M, acquire semaphore before caller parsing (the memory-heavy phase).

```python
# In _process_one:
rescue_df = parse_fn(rescue_path)
n = len(rescue_df)
is_large = n > 2_000_000
if is_large:
    large_sem.acquire()
    print(f"  [{sample_id}] Large sample ({n} variants) — processing exclusively")
try:
    # ... caller parsing, join, write parquet, free
finally:
    if is_large:
        large_sem.release()
```

### 2. Avoid to_list() + set(zip()) — pass rescue parquet to Rust

**Decision:** Write the 4 target columns to a temp parquet file. Add a Rust function that reads the parquet and builds the target HashSet internally. Zero Python objects.

```python
# OLD (3.25 GB Python objects):
chroms = df["CHROM"].to_list()  # 7M strings
poss = df["POS"].to_list()      # 7M ints
refs = df["REF"].to_list()      # 7M strings
alts = df["ALT"].to_list()      # 7M strings
targets = set(zip(chroms, poss, refs, alts))  # 7M tuples

# NEW (0 Python objects):
tmp = f"/dev/shm/{sample_id}_targets.parquet"
df.select(["CHROM", "POS", "REF", "ALT"]).write_parquet(tmp)
# Rust reads tmp, builds HashSet internally
```

**Rust side:** New function `read_targets_from_parquet(path) -> HashSet<(String, i64, String, String)>`. Uses the `parquet` + `arrow` crates (available via pyo3-arrow). Or simpler: use polars in Rust (but that requires adding the polars crate).

**Even simpler approach:** Add a `#[pyfunction]` that accepts a Python list-of-lists and builds the HashSet in Rust with `py.detach()`. But that still requires Python lists...

**Pragmatic approach for this change:** Instead of full parquet round-trip, use chunked Arrow access:

```python
# Use polars' to_arrow() for zero-copy column access
chroms_arr = df["CHROM"].to_arrow()  # PyArrow array, zero-copy
poss_arr = df["POS"].to_arrow()
refs_arr = df["REF"].to_arrow()
alts_arr = df["ALT"].to_arrow()
# Pass Arrow arrays to Rust via pyo3-arrow
```

But pyo3-arrow integration requires API knowledge. Let me use a simpler interim solution:

**Interim: chunked iteration with immediate set insertion, then free:**
```python
targets = set()
for chunk_start in range(0, n, 100_000):
    chunk_end = min(chunk_start + 100_000, n)
    for chrom, pos, ref, alt in df[chunk_start:chunk_end, ["CHROM","POS","REF","ALT"]].iter_rows():
        targets.add((chrom, pos, ref, alt))
```

This creates 100K tuples at a time instead of 7M all at once. Peak Python tuple overhead: ~30 MB per chunk vs 2 GB all at once.

### 3. Avoid Rust HashSet clone

**Decision:** Change `parse_caller_vcf` to accept `&HashSet` (borrow) instead of taking ownership. Already done — the function takes `target_positions: &HashSet<...>`. But line 155 creates a clone:

```rust
let mut remaining = target_positions.clone();  // 2 GB clone
```

**Fix:** Remove the clone. Remove elements from a separate `HashSet` or use a visited counter approach. Since we match by position and the VCF is sorted, we could use a different approach: iterate the target set in order, advance through the VCF. But that's a bigger refactor.

**Simpler fix:** Instead of cloning the full set, clone only the keys into a `HashSet` of integer hashes or use a bitset. But the simplest: just remove the clone and pass a mutable reference.

Actually, `target_positions` is `&HashSet` — we can't remove from it. We need a mutable set. Two options:
a) Pass ownership: `target_positions: HashSet<...>` (consume, no clone needed)
b) Use a separate removal set: `let mut found = HashSet::new(); found.insert(key);`

Option (a) is simplest: change the function signature to take `HashSet` (owned), not `&HashSet`. The caller moves the set in. No clone.

### 4. Add gc.collect() at key points

**Decision:** Call `gc.collect()` after the `del` statements in `process_single_sample`:

```python
del chroms, poss, refs, alts  # or equivalent
del target_positions
del rescue_df, caller_data
gc.collect()
```

## Implementation Plan

1. Auto-throttle (cli.py `_process_one`) — Semaphore + variant count check
2. Chunked target set building (cli.py) — replaces to_list() + set(zip())
3. Remove Rust HashSet clone (caller.rs) — take owned HashSet
4. Add gc.collect() at free points (cli.py)
5. Tests
