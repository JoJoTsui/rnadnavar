## Context

The `fix-memory-retention-column-oriented-and-process-isolation` change introduced `multiprocessing.Pool` with `spawn` context for process-isolated sample processing. The `spawn` context creates fresh Python interpreters per worker, which causes POSIX named semaphores (`SemLock`) from the Pool's internal queues to be registered with each child's `resource_tracker`. When children exit after processing one sample (`maxtasksperchild=1`), the child-side tracker warns about 6 semaphores that the parent created but hasn't `sem_unlink`'d yet.

## Goals / Non-Goals

**Goals:**
- Eliminate the `resource_tracker: There appear to be 6 leaked semaphore objects` warning
- Maintain identical memory isolation (worker process exits → kernel reclaims all)
- Maintain or improve worker startup speed

**Non-Goals:**
- Changing the worker function or sample processing logic
- Adding cross-process synchronization (no longer needed with per-process isolation)
- Changing any other multiprocessing usage in the codebase

## Decisions

### Decision: Switch from `spawn` to `fork` context

**Choice:** `mp.get_context("fork")` instead of `mp.get_context("spawn")`.

**Rationale:**

With `fork`:
- Child inherits parent's memory via COW (copy-on-write)
- Semaphore file descriptors are inherited directly — no `sem_open()` call needed
- No `resource_tracker` registration in child processes
- Worker startup is ~1ms (fork syscall) vs ~3s (spawn + import all modules)
- `maxtasksperchild=1` ensures child exits after one sample → kernel reclaims all memory

The original concern about `fork` (COW memory pressure from CPython reference counting) is irrelevant with `maxtasksperchild=1`: the child exits immediately after processing one sample. Even if CPython dirties pages during processing, those pages are reclaimed by the kernel when the process exits. The parent's memory is unaffected because the parent never forks while holding large DataFrames (it forks workers at pool creation time, before any samples are processed).

**Alternatives considered:**
- **Keep spawn + suppress warning**: `warnings.filterwarnings('ignore', 'resource_tracker')`. Rejected — hides real leaks.
- **Keep spawn + `pool.terminate()` after `close()`+`join()`**: Tried, doesn't fix the child-side tracker registration.
- **Keep spawn + null internal queues + GC**: Accessing `pool._inqueue` etc. is fragile (private API).
- **Manual process management with `subprocess`**: More code, same result. Rejected — overengineered.

## Risks / Trade-offs

- **[Risk] `fork()` in multi-threaded parent**: If the parent has active threads at fork time, only the calling thread survives in the child. Any locks held by other threads are permanently locked in the child. → **Mitigation**: The pool is created at startup before any sample processing. The only threads at that point are polars' internal rayon pool (idle at startup) and the main thread. polars handles fork safety via `os.register_at_fork()`. No known issues.

- **[Trade-off] Fork briefly duplicates parent RSS**: The kernel reserves swap for the parent's RSS at fork time (~2-3 GB for the Python process with all modules loaded). → **Acceptable**: With a 200 GB cgroup, 4 workers × 3 GB reserved = 12 GB. The COW pages are never actually allocated unless dirtied, and actual RSS drops rapidly as workers run their own processing.
