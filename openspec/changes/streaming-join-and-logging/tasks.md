## 1. Column pruning (rescue_parser.py + cli.py)

- [x] 1.1 Define `SLIM_COLS` in cli.py — the 12 columns needed during processing
- [x] 1.2 Split rescue_df into rescue_slim + rescue_output after parsing
- [x] 1.3 Free full rescue_df immediately after split
- [x] 1.4 hstack rescue_output back before parquet write
- [x] 1.5 Remove 6 dropped columns from rescue_parser.py field list

## 2. Streaming join (cli.py + caller_parser.py)

- [x] 2.1 Extract `_join_one_caller` function from `join_caller_columns` — joins a single caller
- [x] 2.2 Replace `parse_all_callers` + `join_caller_columns` with streaming loop in process_single_sample
- [x] 2.3 After each caller join: `del cols; gc.collect()`
- [x] 2.4 Remove caller_data accumulation dict

## 3. Memory logging (cli.py)

- [x] 3.1 Add `_mem()` function using `/proc/self/status` (no psutil dep)
- [x] 3.2 Place `_mem()` at every step boundary in process_single_sample
- [x] 3.3 Place `_mem()` in `_process_one` after parquet write + free
- [x] 3.4 Use `flush=True` on all print statements

## 4. Log level separation (cli.py + caller_parser.py)

- [x] 4.1 Add `--verbose` CLI flag
- [x] 4.2 Default: show only sample-level progress + [MEM] lines
- [x] 4.3 `--verbose`: show per-caller scanning/found messages
- [x] 4.4 Remove duplicate "[N/M] SID - processing..." lines (redundant with progress)

## 5. Tests

- [x] 5.1 Add test_column_pruning — verify slim+output columns cover all original columns
- [x] 5.2 Add test_streaming_join — verify single-caller join produces same result as batch join
- [x] 5.3 Add test_mem_logging — verify _mem() exists and is called
- [x] 5.4 Add test_dropped_columns — verify 6 columns not in rescue output
- [x] 5.5 Run full test suite, ensure 0 regressions
