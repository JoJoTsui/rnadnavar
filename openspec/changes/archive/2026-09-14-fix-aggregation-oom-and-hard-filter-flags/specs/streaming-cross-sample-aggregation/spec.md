# streaming-cross-sample-aggregation

## Purpose

Cross-sample aggregation functions use column-projected lazy polars queries that only read needed columns from per-sample parquet files, avoiding full dataset materialization.

## ADDED Requirements

### Requirement: Confidence tier computation uses column-pruned collect
The confidence tier and soft flags computation SHALL collect only the columns required by `compute_confidence_tier` and `compute_soft_flags` (~15 columns), rather than all columns in the lazy scan. After computing confidence tiers and soft flags eagerly on the pruned DataFrame, the results SHALL be joined back to the main lazy scan via a left join on (sample_id, CHROM, POS).

#### Scenario: Pruned collect for confidence tier
- **WHEN** confidence tiers need to be computed from a lazy scan over 66 parquet files with 7.9M variants and 165 columns
- **THEN** only the columns listed in `_CONFIDENCE_TIER_COLS` SHALL be materialized by the `.collect()` call
- **AND** peak RSS during the collect SHALL NOT exceed 5 GB

#### Scenario: Join-back preserves lazy scan
- **WHEN** confidence tier and soft flags are joined back to the main lazy scan
- **THEN** the resulting `combined_df` SHALL remain a `pl.LazyFrame`
- **AND** downstream aggregation functions SHALL be able to scan from parquet without re-materializing

## REMOVED Requirements

### Requirement: Re-scan after confidence tier computation
**Reason**: The re-scan block (previously at cli.py lines 1251–1275) re-read all parquet files, re-applied filters, and wrapped the in-memory DataFrame as `.lazy()`. This was a workaround for the full `.collect()` polluting the lazy scan. With column-pruned collect + join-back, the lazy scan stays clean and the re-scan is unnecessary.
**Migration**: Delete the re-scan block. The join-back mechanism replaces its function with less I/O and lower memory.
