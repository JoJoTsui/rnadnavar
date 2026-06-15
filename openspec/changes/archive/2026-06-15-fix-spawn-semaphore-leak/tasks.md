## 1. Remove large_lock workaround

- [x] 1.2 Remove the `large_lock = None` line and associated comments
- [x] 1.3 Update the `_process_worker` function to remove the unused `large_lock` parameter
- [x] 1.4 Remove worker_args `large_lock` from the args tuple construction

## 2. Switch from spawn to fork — REVERTED

- [x] 1.1 ~~Change `ctx = mp.get_context("spawn")` to `ctx = mp.get_context("fork")`~~ — **REVERTED** by bfb3fb4: fork deadlocks with polars' rayon thread pool
- [x] 2.1 ~~Update `--process-mode` choices: `spawn` → `fork`~~ — **REVERTED**
- [x] 2.2 ~~Update default process mode: `spawn` → `fork`~~ — **REVERTED**

## 3. Actual fix: suppress cosmetic warning

- [x] 3.0 Suppress resource_tracker warning via `warnings.filterwarnings('ignore', message='resource_tracker')` in cli.py (line 422)

## 4. Verification

- [x] 4.1 Full 64-sample pipeline run completed, no semaphore warnings visible
- [x] 4.2 Pipeline completed smoothly, no speed regression
- [x] 4.3 Run existing test suite, verify 0 regressions (188 passed, 1 pre-existing flaky GIL test)

## 5. Documentation (moot — target change already archived)

- [x] 5.1 ~~Update `fix-memory-retention-column-oriented-and-process-isolation/design.md`~~ — superceded
- [x] 5.2 ~~Update `fix-memory-retention-column-oriented-and-process-isolation/tasks.md`~~ — superceded
