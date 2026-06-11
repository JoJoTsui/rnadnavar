## Why

Samples with >2M variants (4071: 7M, 4072: 3.5M, 4069: 2.7M) cause memory spikes of 20-30 GB each when processed concurrently with `--sample-workers 4`. The `to_list()` calls create Python string/int objects for every variant (28M strings for a 7M-variant sample), and the `set(zip(...))` creates 7M Python tuples. With 4 concurrent large samples, peak memory exceeds cgroup limits. Auto-throttling large samples and avoiding Python object materialization fixes this.

## What Changes

- **Auto-throttle large samples**: Use a `threading.Semaphore` to limit concurrent processing of samples with >2M positions to 1 at a time. Small samples process normally with full parallelism.
- **Avoid `to_list()` for target positions**: Instead of materializing Python lists (`chroms = df["CHROM"].to_list()`), iterate the DataFrame columns in chunks to build the Rust HashSet. This eliminates 28M Python strings per 7M-variant sample.
- **Tests**: Add memory tests verifying auto-throttle behavior and no to_list regression.

## Capabilities

### New Capabilities
- `large-sample-throttling`: Automatic concurrency limiting when sample variant count exceeds threshold

### Modified Capabilities
- None — implementation change only, no spec-level behavior change

## Impact

- **Python code**: `cli.py` (auto-throttle + chunked target position building)
- **Rust code**: None (Rust functions unchanged — just called differently)
- **Tests**: `test_seq2neo_stats.py` (throttle + memory tests)
- **No new dependencies**
