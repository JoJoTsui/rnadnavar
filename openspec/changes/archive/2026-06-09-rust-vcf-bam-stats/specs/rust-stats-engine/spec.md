## ADDED Requirements

### Requirement: Rust aggregation engine
The system SHALL provide a Rust-based statistics aggregation engine that computes per-variant metrics (VAF = AD_ALT/DP, DNA/RNA mean DP/REF_DP/ALT_DP/VAF across callers) and returns results as Arrow arrays for polars ingestion. The engine SHALL be at least 2× faster than equivalent polars expressions for per-variant computation.

#### Scenario: Per-variant VAF computation
- **WHEN** caller AD_ALT and DP columns are provided
- **THEN** per-caller VAF columns are computed as AD_ALT/DP with null handling for DP=0

#### Scenario: Mean computation across callers
- **WHEN** per-caller metrics exist for DNA and RNA callers
- **THEN** DNA_VAF_mean, RNA_VAF_mean, DNA_DP_mean, RNA_DP_mean, and REF/ALT DP means are computed

### Requirement: Arrow interop with polars
The system SHALL return all computed columns as Arrow record batches that polars can ingest without copy (zero-copy where possible via Arrow C Data Interface).

#### Scenario: Arrow-to-polars handoff
- **WHEN** Rust returns an Arrow RecordBatch
- **THEN** polars creates a DataFrame from the batch without data copy
