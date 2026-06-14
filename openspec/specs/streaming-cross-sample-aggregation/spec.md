# streaming-cross-sample-aggregation

## Purpose

Cross-sample aggregation functions use column-projected lazy polars queries that only read needed columns from per-sample parquet files, avoiding full dataset materialization.

## Requirements

### Requirement: Cross-sample aggregation uses column-projected lazy queries
All cross-sample aggregation functions SHALL use polars lazy queries with explicit column selection (`.select(needed_cols)`) pushed into the parquet scan before `.collect()`. No function SHALL call `.collect()` on the full lazy frame without prior column projection and aggregation.

#### Scenario: Aggregation loads only required columns
- **WHEN** `disease_summary` is called with a lazy frame over 24 parquet files
- **THEN** only the columns needed for disease-level aggregation SHALL be read from disk
- **AND** the `.collect()` SHALL return a small aggregated DataFrame (tens of rows, not millions)

#### Scenario: Aggregation memory is bounded
- **WHEN** any cross-sample aggregation function processes 175M variants
- **THEN** peak RSS during that function SHALL NOT exceed 5 GB
- **AND** memory SHALL be freed before the next aggregation function runs

#### Scenario: Functions handle both eager and lazy frames
- **WHEN** a cross-sample aggregation function receives an already-eager DataFrame
- **THEN** it SHALL compute aggregations without error (no re-collect needed)
- **AND** it SHALL NOT duplicate the data in memory

### Requirement: Streaming aggregation query pattern
Each aggregation function SHALL follow the pattern: select needed columns → filter/transform → group_by → aggregate → collect. The `.collect()` SHALL be the final operation, after all aggregations are defined.

#### Scenario: disease_summary streaming pattern
- **WHEN** `disease_summary(df)` is called
- **THEN** it SHALL select only the columns needed for disease grouping and metric computation
- **AND** it SHALL group by disease and disease_normalized
- **AND** it SHALL call `.collect()` after all aggregations are defined

#### Scenario: tier_summary streaming pattern
- **WHEN** `tier_summary(df)` is called
- **THEN** it SHALL select only the columns needed for tier-level aggregation
- **AND** it SHALL call `.collect()` after all aggregations are defined
