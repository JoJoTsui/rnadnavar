## Context

The seq2neo stats pipeline currently has two intertwined problems in the cross-sample aggregation phase (`cli.py` lines 1148–1275):

1. **OOM at `.collect()`**: After building a lazy scan over all 66 per-sample parquets and applying unified + hard filters, line 1223 calls `combined_df.collect()` which materializes all 165 columns × ~7.9M rows into a single polars DataFrame. This consumes ~219 GB RSS and triggers the cgroup OOM killer (confirmed via dmesg: `total-vm:281GB anon-rss:219GB`).

2. **Hard filter drops variants**: The hard filter (`build_hard_filter_expr`) drops ~3.1M variants (28%) before confidence tiering. This conflates two concerns — statistics should observe and characterize, while filtering belongs in downstream dataset preparation.

The ironic part: `_ensure_eager` (in `statistics.py`) already implements column pruning to load only ~60 `_CROSS_SAMPLE_COLS` — but it never gets a chance because line 1223 already materialized everything.

## Goals / Non-Goals

**Goals:**
- Eliminate the OOM kill by never materializing all 165 columns into memory at once
- Keep `combined_df` as a true `LazyFrame` through the entire aggregation pipeline
- Convert hard filter from dropping variants to generating flag columns (observe, don't exclude)
- Externalize hard filter conditions as a Python config module for maintainability
- Add hard filter breakdown statistics and visualization
- Remove dead code (re-scan block, `--no-hard-filter` flag)

**Non-Goals:**
- Rewrite `compute_confidence_tier` / `compute_soft_flags` to operate lazily (eager is fine if column-pruned)
- Change per-sample parquet format (flag columns are added lazily at aggregation time)
- Modify the unified filter pipeline (NoConsensus exclusion etc. remains as-is)
- Change BAM statistics or visualization architecture beyond the hard filter chart

## Decisions

### Decision 1: Column-pruned eager collect for confidence tier computation

**Choice**: Collect only the ~15 columns needed by `compute_confidence_tier` + `compute_soft_flags`, compute eagerly, then join `confidence_tier` and `soft_flags` back to the lazy scan.

**Alternatives considered**:
- *Full lazy rewrite of confidence tier*: Would require rewriting the conditional `pl.when().then().otherwise()` logic as lazy window operations. Higher risk of correctness drift, harder to verify parity with current output.
- *Compute confidence tier per-sample before aggregation*: Would require storing confidence tiers in per-sample parquets. More invasive change, but better long-term. Deferred to future iteration.
- *Batch the collect into chunks of N samples*: Addresses OOM but adds complexity (dedup, merge). Lazy scan already handles cross-sample window functions (`n_recurrent_samples`).

**Rationale**: Pruned collect is the smallest safe change. The ~15 columns (flag booleans, tier columns, FILTER, N_SUPPORT_CALLERS, etc.) for 7.9M rows produce a DataFrame of ~2-3 GB, well within the process's memory limit. All downstream functions already use `_ensure_eager` for column pruning.

### Decision 2: Python config module over YAML

**Choice**: Store hard filter conditions in `bin/vcf_stats/seq2neo/hard_filter_config.py` as a Python data structure with lambda-based expression builders, matching the `bin/common/tier_config.py` pattern.

**Alternatives considered**:
- *YAML config*: More declarative but requires a parser, expression type dispatch, and string-to-lambda translation. Adds complexity without benefit for a config that changes with code.
- *CLI flags per condition*: Unwieldy for 8 conditions, hard to document, no single source of truth.

**Rationale**: The project already uses Python modules for configuration (`tier_config.py`, `vcf_config.py`). Lambdas are the most direct way to express polars column expressions. The config lives next to the code that uses it — no deployment synchronization needed.

### Decision 3: Flag columns computed lazily, joined to lazy scan

**Choice**: Add `flag_hard_*` boolean columns and `hard_filter_flags` / `n_hard_flags` summary columns via `.with_columns()` on the lazy scan. No variants are dropped.

**Rationale**:
- Flag columns travel with the variant data through the lazy pipeline — downstream functions can reference them without re-joining
- `hard_filter_flags` (comma-joined string) and `n_hard_flags` (integer count) enable simple filtering in dataset preparation: `n_hard_flags == 0` for clean variants
- The hard filter breakdown stats use lazy `select(pl.len())` after filtering on each flag — tiny memory cost

### Decision 4: Remove the re-scan block (lines 1251–1275)

**Choice**: Delete the block that re-scanned from parquet, re-applied all filters, recomputed `n_recurrent_samples`, re-applied hard filter, and wrapped `combined_eager.lazy()`.

**Rationale**: This block was a workaround for the fact that `combined_eager` (the fully materialized DataFrame) polluted the lazy scan. With pruned collect + join-back, the main lazy scan stays clean throughout. The re-scan duplicated:
- Parquet I/O (re-reading all files)
- Unified filter application (already done)
- `n_recurrent_samples` window computation (already done)
- Hard filter application (no longer drops)

### Decision 5: Remove `--no-hard-filter` CLI flag

**Choice**: Delete the `--no-hard-filter` argument and the conditional block that gates hard filter application.

**Rationale**: The hard filter no longer drops variants — it only adds flag columns. There is zero performance cost to always computing them, and no reason a user would want to skip flag annotation. If flags are unwanted, they can be ignored downstream.

## Data Flow (New)

```
pl.scan_parquet("*_variants.parquet")        # lazy, zero memory
  .filter(sample_id in current_sample_ids)    # lazy
  .with_columns(partition, n_recurrent_samples)  # lazy
  .filter(unified_filter_expr)                # lazy, NoConsensus etc.
  .with_columns(hard_filter_flag_exprs)       # lazy, 8 flag cols + summary cols
  │
  ├── hard filter stats:
  │     select(pl.len()).filter(flag_hard_X).collect() × 8 → tiny TSV + chart
  │
  ├── confidence tier (pruned):
  │     select(15 cols).collect() → compute_confidence → compute_soft_flags
  │     → join confidence_tier + soft_flags back to lazy scan
  │     → confidence_tier_summary.tsv (small group_by)
  │
  └── downstream stats (all lazy):
        _ensure_eager prunes to ~60 cols per function
        disease_summary, tier_summary, caller_overlap_matrix, etc.
```

## Memory Model

| Stage | Columns | Rows | Est. Memory |
|---|---|---|---|
| Lazy scan + filters | 0 (no materialization) | — | <100 MB (plans only) |
| Hard filter stats | 1 (pl.len()) | 8× scalar | <1 MB |
| Confidence tier collect | ~15 | 7.9M | ~2-3 GB |
| Confidence join-back | +2 (confidence_tier, soft_flags) | — | lazy |
| _ensure_eager per stats fn | ~60 | aggregated (group_by) | <1 GB per fn |

Peak memory: ~3 GB (confidence tier pruned collect), down from 219 GB.

## Risks / Trade-offs

- **[Risk] Schema mismatch if flag columns already exist in parquet** → Mitigation: Use `pl.col()` reference with fallback; flag columns are computed lazily from existing columns, not written to parquet.
- **[Risk] Hard filter conditions reference columns not in all parquets** → Mitigation: `build_hard_filter_flag_exprs()` checks column availability per condition (same pattern as current `build_hard_filter_expr`).
- **[Trade-off] Hard filter breakdown stats add 8 extra `.collect()` calls** → Acceptable: each collects only `pl.len()` (single scalar), total <1ms per condition.
- **[Trade-off] Confidence tier + soft flags still computed eagerly** → Acceptable: ~2-3 GB for 7.9M variants × 15 columns. Future iteration can push this into per-sample parquet generation.
