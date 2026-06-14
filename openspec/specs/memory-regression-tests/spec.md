# memory-regression-tests

## Purpose

Automated tests that prevent memory-expensive patterns from being reintroduced into the seq2neo stats pipeline. Covers DataFrame lifecycle, lazy/eager parity, peak memory bounds, and chart function input patterns.

## Requirements

### Requirement: DataFrame memory freed after per-sample processing
The system SHALL verify that DataFrames are freed from memory after each sample's processing completes and the parquet file is written.

#### Scenario: DataFrame reference deleted
- **WHEN** a sample's variant DataFrame is written to parquet
- **THEN** the `del df` operation SHALL succeed and a subsequent `gc.collect()` SHALL reclaim the memory

### Requirement: Lazy scan parity with eager concat
The system SHALL verify that lazy `pl.scan_parquet()` aggregation produces identical results to eager `pl.concat()` + `pl.DataFrame` aggregation for all cross-sample statistics.

#### Scenario: dataset_summary parity
- **WHEN** 2 samples are processed in both streaming and eager mode
- **THEN** `dataset_summary`, `disease_summary`, `tier_summary`, `sample_tier_summary` outputs SHALL be identical

### Requirement: Peak memory under threshold for multi-sample run
The system SHALL verify that peak memory usage does not exceed 4 GB when processing 20+ samples.

#### Scenario: Memory stays under limit
- **WHEN** 20 samples are processed in streaming mode
- **THEN** peak memory SHALL be less than 4 GB (measured via `tracemalloc` or OS-level monitoring)

### Requirement: Chart functions accept pre-aggregated DataFrames
The system SHALL verify that visualization functions produce valid altair charts when given pre-aggregated summary DataFrames instead of raw variant-level data.

#### Scenario: Chart with aggregated data
- **WHEN** `plot_vc_distribution` receives a DataFrame with columns (set_number, VC, count)
- **THEN** it SHALL produce the same chart as when receiving the raw variant-level DataFrame
