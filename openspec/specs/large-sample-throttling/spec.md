# large-sample-throttling

## Purpose

Auto-throttling of large samples (>2M positions) to prevent concurrent memory spikes, with chunked target position construction and explicit garbage collection.

## Requirements

### Requirement: Large samples are auto-throttled to prevent concurrent memory spikes
The system SHALL limit concurrent processing of samples with more than 2,000,000 variant positions to one at a time. Samples below the threshold SHALL process with full parallelism.

#### Scenario: Large sample acquires exclusive semaphore
- **WHEN** a sample's rescue VCF contains >2,000,000 positions
- **THEN** a semaphore SHALL be acquired before caller VCF parsing begins and released after the sample completes (including on error)

#### Scenario: Small samples bypass throttling
- **WHEN** a sample's rescue VCF contains ≤2,000,000 positions
- **THEN** the semaphore SHALL NOT be acquired and the sample SHALL process concurrently with other samples

### Requirement: Target positions built without materializing full Python lists
The system SHALL NOT call `to_list()` on rescue DataFrame columns to build target positions. Target positions SHALL be built using chunked iteration (≤100,000 rows per chunk) to limit Python object creation.

#### Scenario: Chunked target position construction
- **WHEN** a rescue VCF with 7,000,000 positions is processed
- **THEN** target positions SHALL be built in chunks of ≤100,000 rows, with each chunk's Python objects freed before the next chunk is processed

### Requirement: Rust caller parser avoids HashSet clone
The Rust `parse_caller_vcf` function SHALL accept an owned `HashSet` of target positions (not a reference), eliminating the 2 GB clone for large samples.

#### Scenario: No clone of target positions
- **WHEN** `parse_caller_vcf` is called with 7,000,000 target positions
- **THEN** no clone of the HashSet SHALL be created; the function SHALL consume the owned set directly

### Requirement: Explicit garbage collection at memory free points
The system SHALL call `gc.collect()` after deleting large Python objects in `process_single_sample` to ensure memory is returned to the OS promptly.

#### Scenario: GC after freeing target positions and DataFrames
- **WHEN** `del target_positions` and `del rescue_df, caller_data` are executed
- **THEN** `gc.collect()` SHALL be called immediately after to reclaim memory
