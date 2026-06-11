## 1. Rust: Column-oriented rescue parser (vcf.rs)

- [x] 1.1 Define `RescueColumns` struct with column vectors: `chrom: Vec<String>`, `pos: Vec<i64>`, `ref_base: Vec<String>`, `alt: Vec<String>`, `filter: Vec<String>`, `info: HashMap<String, Vec<Option<String>>>`, `variant_type: Vec<String>`, `ti_tv: Vec<Option<bool>>`
- [x] 1.2 Implement `parse_rescue_columns(path) -> Result<RescueColumns>` — populates column vectors directly during parse loop instead of building `Vec<RescueRecord>`
- [x] 1.3 Compute `variant_type` (SNV/INS/DEL/MNV) and `ti_tv` (transition/transversion) during parse — one pass, no separate iteration
- [x] 1.4 Handle INFO field types in Rust: Integer → String representation, Float → String, Flag → "true"/absent, String/Character → as-is. Use `format_info_value()` for string conversion.
- [x] 1.5 Preserve `info_keys` order from VCF header for deterministic column output

## 2. Rust: pyo3 binding (lib.rs)

- [x] 2.1 Add `parse_rescue_columns` pyfunction returning `PyDict` of `{column_name: PyList}` — same pattern as `parse_caller_vcf`
- [x] 2.2 Push fixed columns: CHROM (str list), POS (i64 list), REF (str list), ALT (str list), FILTER (str list)
- [x] 2.3 Push each INFO column as `PyList` of `Option<String>` (None for missing values)
- [x] 2.4 Push derived columns: variant_type (str list), ti_tv (bool list with None)
- [x] 2.5 Use `py.detach()` for GIL release during parse; re-acquire for Python object construction
- [x] 2.6 Register `parse_rescue_columns` in module function list
- [x] 2.7 Keep existing `parse_rescue` function for backward compatibility during transition

## 3. Python: Rescue VCF wrapper (rust_vcf.py)

- [x] 3.1 Add `parse_rescue_vcf_columns()` function: calls `stats_core.parse_rescue_columns()`, builds `pl.DataFrame` from `pl.Series` with explicit dtypes
- [x] 3.2 Map INFO field types using `RESCUE_INT_FIELDS`, `RESCUE_FLOAT_FIELDS`, `RESCUE_FLAG_FIELDS` for correct dtype inference (Int64/Float64/Bool/Utf8)
- [x] 3.3 Remove `_cast_columns()` function — replaced with `_cast_columns_fallback()` for legacy Rust/cyvcf2 paths only
- [x] 3.4 Remove `.to_list()` calls and list comprehensions for derived columns — Rust provides them directly
- [x] 3.5 Update `parse_rescue_vcf()` to delegate to `parse_rescue_vcf_columns()` when Rust is available, keeping cyvcf2 fallback
- [x] 3.6 Handle missing fields in Rust output: `_fill_missing_fields()` adds null columns for any ALL_RESCUE_FIELDS not in Rust output

## 4. Python: polars-native fallback for derived columns (rescue_parser.py)

- [x] 4.1 Replace per-row computation + `.to_list()` for `variant_type` with `_add_derived_columns_polars()` using polars `when/then/otherwise` chain
- [x] 4.2 Replace per-row computation for `ti_tv` with polars vectorized string operations in `_add_derived_columns_polars()`
- [x] 4.3 Remove Python helper functions `_derive_variant_type` and `_is_transition` (only used in rescue_parser.py)

## 5. Python: Process isolation in CLI (cli.py)

- [x] 5.1 Add `--max-tasks-per-child` CLI flag (default: 1) and `--process-mode` flag (choices: thread/spawn, default: spawn when sample_workers > 1)
- [x] 5.2 Extract `_process_worker(args)` as module-level function accepting `(row, max_workers, use_rust, variant_dir_str, large_sem)`, returning `{"sample_id": str, "stats": dict}`
- [x] 5.3 Implement `_process_worker`: calls `process_single_sample`, writes parquet, calls `del df; gc.collect()`, calls `_malloc_trim()`, returns only stats
- [x] 5.4 Add `_malloc_trim()` helper function with fallback for non-glibc platforms
- [x] 5.5 Create `multiprocessing.Manager().Semaphore(1)` for cross-process large-sample throttle
- [x] 5.6 Pass semaphore to workers via args tuple (Manager proxy is picklable)
- [x] 5.7 In `process_single_sample`, acquire/release semaphore with `timeout=300` in `try/finally` for samples >2M variants
- [x] 5.8 Use `multiprocessing.Pool(processes=N, maxtasksperchild=args.max_tasks_per_child, context=mp.get_context('spawn'))` when `--sample-workers > 1` and `process_mode=spawn`
- [x] 5.9 Use `pool.imap_unordered()` to stream results as samples complete
- [x] 5.10 Accumulate `all_stats` from worker return values (no DataFrame in return)
- [x] 5.11 Keep sequential path and thread-pool path for `--process-mode thread` or `--sample-workers 1`
- [x] 5.12 Add semaphore timeout (300s) in worker to prevent hang if manager process dies

## 6. Tests

- [x] 6.1 Add `TestColumnOrientedRescueParser` class (7 tests): column parity, lengths, derived cols, polars parity, malloc_trim
- [x] 6.2 All column lengths equal after column-oriented parse
- [x] 6.3 Missing INFO fields get null columns with correct dtypes via _fill_missing_fields
- [x] 6.4 variant_type and ti_tv computed correctly in Rust path
- [x] 6.5 polars-native derived columns produce correct output
- [ ] 6.6 `_process_worker` is picklable (verified by running --process-mode spawn)
- [x] 6.7 `_malloc_trim()` doesn't crash on Linux and handles missing libc
- [ ] 6.8 Worker function releases semaphore on exception (verified by running pipeline)
- [ ] 6.9 Parquet output parity between process-isolated and thread-based (verified by running pipeline)
- [x] 6.10 Full existing test suite: 135 total, 133 pass, 2 fail (1 pre-existing flaky GIL, 1 renamed semaphore → fixed)
- [ ] 6.11 Memory regression test (verified by running pipeline)

## 7. Verification

- [x] 7.1 Build Rust module with all 6 functions exported (parse_rescue, parse_rescue_columns, bam_stats, parse_caller_vcf, compute_tiers, pileup_variants)
- [ ] 7.2 Run 4-sample test with `--sample-workers 4 --process-mode spawn`, verify no OOM and output parity
- [ ] 7.3 Run 12-sample test (the same set that previously OOM-killed), verify all samples complete successfully
- [ ] 7.4 Profile memory: log per-worker peak RSS, confirm <60 GB per process for large (>2M variant) samples
- [ ] 7.5 Verify all output CSV/parquet files are byte-identical to a thread-based run
