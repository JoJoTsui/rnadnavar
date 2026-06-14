## MODIFIED Requirements

### Requirement: Cross-sample aggregation functions accept lazy frames
Cross-sample aggregation functions (`dataset_summary`, `disease_summary`, `set_summary`, `sample_tier_summary`) SHALL accept both `pl.DataFrame` (eager) and `pl.LazyFrame` inputs. When a lazy frame is provided, the function SHALL materialize it once via `.collect()` before performing multiple filter/group-by operations.

#### Scenario: dataset_summary with lazy frame
- **WHEN** `dataset_summary(lazy_frame)` is called with a `pl.LazyFrame` from `pl.scan_parquet()`
- **THEN** it SHALL call `.collect()` once and produce identical output to `dataset_summary(eager_frame)`

#### Scenario: Backward compatible with eager frame
- **WHEN** `dataset_summary(eager_df)` is called with a `pl.DataFrame`
- **THEN** it SHALL work unchanged (no `.collect()` call needed)
