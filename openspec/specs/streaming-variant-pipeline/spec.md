# streaming-variant-pipeline

## Purpose

Streaming architecture for the seq2neo variant statistics pipeline — per-sample parquet output, lazy cross-sample aggregation via polars scan_parquet, and memory-safe sample statistics accumulation.

## Requirements

### Requirement: Per-sample variant data streamed to parquet files
The system SHALL write each sample's variant details DataFrame to a per-sample parquet file immediately after processing and free the DataFrame memory before processing the next sample. At no point SHALL all samples' DataFrames be held simultaneously in memory.

#### Scenario: Single sample processed and freed
- **WHEN** a sample's variant DataFrame is computed
- **THEN** the DataFrame SHALL be written to `variant_details/<sample_id>_variants.parquet` and the Python reference SHALL be deleted (del) and garbage-collected before the next sample begins processing

#### Scenario: 65 samples processed without OOM
- **WHEN** 65 samples are processed sequentially or with sample-workers
- **THEN** peak memory SHALL be less than 4 GB

### Requirement: Cross-sample aggregation via lazy parquet scan
The system SHALL use `pl.scan_parquet()` to lazily load per-sample parquet files for cross-sample aggregation, avoiding materialization of the full combined dataset in memory.

#### Scenario: Lazy scan aggregates correctly
- **WHEN** per-sample parquet files exist in the variant_details directory
- **THEN** `pl.scan_parquet("variant_details/*.parquet")` SHALL produce identical aggregation results to loading all DataFrames eagerly and concatenating

### Requirement: Sample statistics accumulate without variant data
The system SHALL accumulate only per-sample summary statistics (`sample_summary` dicts) across samples, NOT the raw variant DataFrames. Only after all samples are processed SHALL cross-sample aggregation begin.

#### Scenario: all_stats list is memory-safe
- **WHEN** 65 samples are processed
- **THEN** the `all_stats` list SHALL contain 65 dicts of ~1 KB each, totaling less than 1 MB
