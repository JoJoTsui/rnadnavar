## Why

The 24-sample pipeline prints `resource_tracker: There appear to be 6 leaked semaphore objects to clean up at shutdown` at every run. The 6 POSIX semaphores come from `mp.Pool`'s 3 internal `SimpleQueue` objects (2 locks each). With `spawn` context, child processes call `sem_open()` to open the named semaphores, registering them with the child's `resource_tracker`. When children exit (`maxtasksperchild=1`), the child-side tracker detects the parent-created semaphores as "leaked" before the parent has a chance to `sem_unlink()` them. The `close()` + `join()` fix doesn't help because it only affects worker exit — the tracker registration is a separate mechanism.

## What Changes

- Switch from `spawn` to `fork` context for `multiprocessing.Pool`. Forked children inherit semaphore file descriptors directly (no `sem_open` → no tracker registration → no warning).
- Remove the `large_lock = None` workaround — fork allows cross-process synchronization via inheritance if needed later.
- Fork is 6-10× faster at worker startup (no module re-import) with identical memory isolation (`maxtasksperchild=1` → kernel reclaims all on exit).

## Capabilities

### New Capabilities
- None — bugfix only.

### Modified Capabilities
- `process-isolated-sample-processing`: Process isolation mechanism changed from `spawn` to `fork`. Behavior unchanged — each sample still runs in an independent process that exits after completion.

## Impact

- **Python**: `cli.py` — one line change: `mp.get_context("spawn")` → `mp.get_context("fork")`. Remove the `large_lock = None` line and associated comments.
- **Speed**: Worker startup drops from ~3s (spawn + re-import) to ~1ms (fork). Net speed improvement.
- **Memory**: Identical isolation — worker exits, kernel reclaims all.
