## 1. Switch from spawn to fork

- [x] 1.1 Change `ctx = mp.get_context("spawn")` to `ctx = mp.get_context("fork")` in cli.py
- [x] 1.2 Remove the `large_lock = None` line and associated comments (no longer needed — fork inherits semaphores cleanly)
- [x] 1.3 Update the `_process_worker` function to remove the unused `large_lock` parameter
- [x] 1.4 Remove worker_args `large_lock` from the args tuple construction

## 2. Update CLI flags

- [x] 2.1 Update `--process-mode` choices and help text: `spawn` → `fork`
- [x] 2.2 Update default process mode: `spawn` → `fork` when `--sample-workers > 1`

## 3. Update documentation

- [ ] 3.1 Update `fix-memory-retention-column-oriented-and-process-isolation/design.md` Decision 2 to reflect the fork→spawn correction
- [ ] 3.2 Update `fix-memory-retention-column-oriented-and-process-isolation/tasks.md` with these changes

## 4. Tests

- [x] 4.1 Fork pool test: 8 results, no resource_tracker warning ✓
- [ ] 4.2 Run existing test suite, verify 0 regressions

## 5. Verification

- [ ] 5.1 Run 24-sample pipeline, confirm "6 leaked semaphore objects" warning is gone
- [ ] 5.2 Verify processing time is not slower (fork startup is faster, so wall time should decrease)
