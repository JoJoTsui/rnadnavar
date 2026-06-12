## Context

The seq2neo variant statistics pipeline processes 65 samples, each with ~200K variants × ~200 columns (~600 MB per sample DataFrame). The current architecture accumulates all DataFrames in `samples_data` dict, concatenates to `combined_df`, then runs aggregations and visualization on the combined DataFrame. Peak memory is ~80 GB, causing OOM kills at 12+ samples. The fix streams each sample to disk and uses lazy polars scans for cross-sample operations.

## Goals / Non-Goals

**Goals:**
- Peak memory < 2 GB for 65 samples (down from ~80 GB)
- All output CSV files and dashboard HTML identical to current output
- Same or better performance (lazy scans may be faster than eager concat)
- Backward compatible: existing single-sample tests pass unchanged

**Non-Goals:**
- Changing the output format or file structure
- Adding new dependencies
- Process-level parallelism for cross-sample aggregation (polars handles parallelism internally)

## Decisions

### 1. Per-sample parquet files (not partitioned, not Hive)

**Decision:** Write `variant_details_{sample_id}.parquet` to a subdirectory, one file per sample.

**Rationale:** polars `scan_parquet("*.parquet")` natively handles schema unification across files with different column sets (some callers may have SB/FAD while others don't). Partitioned writes would require consistent schemas and add complexity. One file per sample is simple, debuggable, and enables partial re-processing.

**Alternative considered:** Write a single partitioned parquet dataset with `pl.scan_parquet()` and partition columns. Rejected because column sets vary across callers, and partitioned writes need consistent schema.

### 2. Lazy frame API pattern: accept both eager and lazy

**Decision:** Aggregation functions check `isinstance(df, pl.LazyFrame)` and call `.collect()` once internally. Return types unchanged.

**Rationale:** Avoids breaking all callers. Existing tests that pass eager DataFrames continue to work. The CLI path uses lazy frames from `scan_parquet()`, while unit tests can use eager DataFrames.

```python
def dataset_summary(df: pl.LazyFrame | pl.DataFrame) -> dict:
    if isinstance(df, pl.LazyFrame):
        df = df.collect()
    # ... rest unchanged
```

### 3. Pre-aggregate for visualization, not sample

**Decision:** For charts that need raw data (scatter plots: DNA vs RNA VAF/DP), sample 10K rows from the lazy frame. For charts that aggregate (bar charts, box plots), pre-compute aggregations before passing to chart functions.

**Rationale:** Sampling 10K rows from 13M is statistically sufficient for scatter plots. Pre-aggregating for bar charts means altair receives DataFrames with 10-50 rows instead of millions. This eliminates `.to_pandas()` conversion overhead.

**Charts needing raw data (sampled):**
- `plot_dna_vs_rna_vaf` — sample 10K rows
- `plot_dna_vs_rna_dp` — sample 10K rows
- `plot_ref_alt_dp_scatter` — sample 10K rows
- `plot_dna_vs_rna_per_caller` — sample per caller

**Charts needing aggregation (pre-computed):**
- All bar charts, box plots, pie charts, histograms — pre-group

### 4. join_caller_columns: Series construction, not row dicts

**Decision:** Build `pl.Series` directly from column vectors, then construct DataFrame from Series list. Each Series specifies its dtype explicitly, avoiding schema inference.

**Rationale:** `pl.Series("DP", values, dtype=pl.Int64)` creates the Arrow array directly from the Python list, then the list can be freed. This is ~2× more memory-efficient than `pl.DataFrame(list_of_dicts)` which creates intermediate Python dicts for every row.

### 5. Visualization: remove .to_pandas(), use polars-native altair

**Decision:** All `.to_pandas()` calls removed. altair natively supports polars DataFrames.

**Rationale:** polars DataFrames pass directly to altair without conversion. The one exception is the validation heatmap which uses pandas `.pivot()` — replaced with polars `.pivot()`.

## Risks / Trade-offs

- **[Risk] Parquet write I/O per sample**: Writing 65 parquet files adds disk I/O. → **Mitigation:** Parquet write is ~1-2s per sample (compressed). Total ~2 minutes for 65 samples, negligible vs ~10 min processing time.
- **[Risk] Lazy scan schema mismatch**: Different samples may have different columns (missing fields for some callers). → **Mitigation:** `pl.scan_parquet` with `how="diagonal_relaxed"` handles schema differences automatically.
- **[Risk] Chart output changes with pre-aggregation**: Pre-computed quartiles may differ slightly from altair's internal boxplot computation. → **Mitigation:** Verify chart PNG outputs match within visual tolerance using the test suite.
- **[Risk] Memory regression if someone reverts to eager**: → **Mitigation:** `TestMemoryEfficiency` in CI catches regression.

## Migration Plan

1. Fix correctness bugs first (low risk, independent)
2. Refactor cli.py to streaming (core change)
3. Update statistics.py for lazy frames
4. Update visualizer.py for pre-aggregated inputs
5. Add memory regression tests
6. Run on 4-sample test set, verify output parity
7. Run on 20-sample set, verify < 2 GB peak memory
8. Run on full 65-sample set
