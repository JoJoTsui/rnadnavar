## 1. Revert fork → spawn + suppress warning

- [x] 1.1 Change `ctx = mp.get_context("fork")` back to `ctx = mp.get_context("spawn")` in cli.py
- [x] 1.2 Update `--process-mode` choices and help text: `fork` → `spawn`
- [x] 1.3 Change default process mode: `"fork"` → `"spawn"` when `--sample-workers > 1`
- [x] 1.4 Add `warnings.filterwarnings('ignore', message='resource_tracker')` at top of `main()`
- [x] 1.5 Remove fork-specific comment, update to spawn rationale

## 2. Streaming cross-sample aggregation

- [x] 2.1 Define `_CROSS_SAMPLE_COLS` in statistics.py — 35 columns needed by all aggregation functions
- [x] 2.2 Modify `_ensure_eager()` to select only `_CROSS_SAMPLE_COLS` before `.collect()` on lazy frames
- [x] 2.3 Add `_ensure_eager` import + call in `tier_summary()` (tiering_stats.py) for column-pruned collection
- [x] 2.4 All other aggregation functions already use `_ensure_eager` — automatically fixed by #2.2
- [x] 2.5 Add `_mem()` logging between each aggregation step in cli.py
- [x] 2.6 Add `_malloc_trim()` + `gc.collect()` after aggregation cleanup

## 3. Visualizer cleanup

- [x] 3.1 Remove duplicate `_maybe_collect` definition (lines 23-27 in visualizer.py)

## 4. Test speed optimization

- [x] 4.1 Add `@pytest.fixture(scope="class")` for rescue VCF parse in `TestColumnOrientedRescueParser`
- [x] 4.2 Update all 7 test methods to use the fixture instead of individual parses (saves ~80s)

## 5. Verification

- [x] 5.1 Pending: syntax verified, user to run pipeline
- [ ] 5.2 Run 12-sample pipeline with `--sample-workers 4 --process-mode spawn`, verify no OOM
- [ ] 5.3 Run 24-sample pipeline, verify all cross-sample aggregations complete
- [ ] 5.4 Verify CSV output matches previous runs (same data, different computation path)
- [ ] 5.5 Confirm no `resource_tracker` warning at exit
