## 1. Auto-throttle large samples (cli.py)

- [x] 1.1 Add `threading.Semaphore(1)` for large sample mutual exclusion
- [x] 1.2 In `_process_one`, check variant count after rescue parse, acquire semaphore if >2M
- [x] 1.3 Use try/finally to ensure semaphore release even on exception
- [x] 1.4 Add log message when large sample is throttled

## 2. Chunked target position building (cli.py)

- [x] 2.1 Replace `to_list()` × 4 + `set(zip())` with chunked iteration
- [x] 2.2 Use 100K chunk size, free each chunk after processing
- [x] 2.3 Remove `del chroms, poss, refs, alts` (no longer needed)

## 3. Remove Rust HashSet clone (caller.rs)

- [x] 3.1 Change `parse_caller_vcf` to accept owned `HashSet` instead of `&HashSet`
- [x] 3.2 Remove `let mut remaining = target_positions.clone()`
- [x] 3.3 Update lib.rs pyfunction to pass owned HashSet

## 4. Add gc.collect() at free points (cli.py)

- [x] 4.1 Add `gc.collect()` after `del target_positions` in process_single_sample
- [x] 4.2 Add `gc.collect()` after `del rescue_df, caller_data` in process_single_sample
- [x] 4.3 Add `import gc` at top of cli.py

## 5. Tests

- [x] 5.1 Add `test_chunked_target_building` — verify chunked iteration vs set(zip)
- [x] 5.2 Add `test_rust_no_hashset_clone` — verify function takes owned HashSet
- [x] 5.3 Add `test_gc_collect_in_process_sample` — code review check
- [x] 5.4 Run full test suite, ensure 0 regressions
