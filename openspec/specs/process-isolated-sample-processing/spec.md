# process-isolated-sample-processing

## Purpose

Each sample is processed in an independent Python process using `multiprocessing.Pool` with `maxtasksperchild=1` and `context='spawn'`, ensuring the OS reclaims all memory after each worker exits. Includes cross-process large-sample throttling and `malloc_trim` calls.

## Requirements

### Requirement: Each sample processed in an independent Python process
The pipeline SHALL process each sample in a separate Python process using
`multiprocessing.Pool` with `maxtasksperchild=1` and `context='spawn'`. When a worker
process completes its sample, it SHALL exit, and the operating system SHALL reclaim
all memory allocated by that process.

#### Scenario: Memory reclaimed after sample completes
- **WHEN** a worker process finishes processing a sample and exits
- **THEN** the process RSS attributable to that sample SHALL return to the operating system
- **AND** subsequent samples processed by new workers SHALL start with a clean memory slate

#### Scenario: Concurrent samples in separate processes
- **WHEN** 4 worker processes are processing 4 samples concurrently
- **THEN** each process SHALL have its own independent memory space
- **AND** memory freed in one process SHALL NOT affect memory available to other processes

#### Scenario: Sequential mode fallback
- **WHEN** `--sample-workers 1` is specified
- **THEN** samples SHALL be processed sequentially in the main process
- **AND** no subprocess SHALL be spawned

#### Scenario: resource_tracker warning suppressed
- **WHEN** the pipeline starts
- **THEN** `warnings.filterwarnings('ignore', message='resource_tracker')` SHALL be called
- **AND** no "leaked semaphore objects" warning SHALL appear at shutdown

### Requirement: Cross-process large-sample throttle
Samples with more than 2 million variants SHALL be processed exclusively (one at a time)
using a cross-process semaphore from `multiprocessing.Manager`. Other samples MAY be
processed concurrently.

#### Scenario: Large sample acquires exclusive access
- **WHEN** a sample has >2M positions in the rescue VCF
- **THEN** the worker SHALL acquire the cross-process semaphore before proceeding with the streaming join
- **AND** no other large sample SHALL process concurrently while the semaphore is held

#### Scenario: Large sample releases semaphore on completion
- **WHEN** a large sample completes (successfully or with error)
- **THEN** the worker SHALL release the cross-process semaphore in a `finally` block

#### Scenario: Small samples process concurrently with large sample
- **WHEN** a large sample holds the exclusive semaphore
- **THEN** small samples (<2M positions) SHALL continue processing concurrently without acquiring the semaphore

### Requirement: Worker function is picklable and importable
The worker function passed to `multiprocessing.Pool` SHALL be defined at module level (not nested) and SHALL accept only picklable arguments. It SHALL NOT reference any thread-local or process-unsafe state.

#### Scenario: Worker receives manifest row dict
- **WHEN** a worker process is started
- **THEN** it SHALL receive a plain Python dict containing the sample manifest row
- **AND** all values in the dict SHALL be picklable types (str, int, float, bool, None)

#### Scenario: Worker returns only small result
- **WHEN** a worker process completes
- **THEN** it SHALL return only a dict with `sample_id` (str) and `stats` (dict of scalars)
- **AND** the returned dict SHALL NOT contain any DataFrame, Series, or file handle

### Requirement: Output format unchanged
The per-sample parquet files and aggregated statistics produced by process-isolated processing SHALL be identical to those produced by thread-based processing (given the same inputs and the same column-oriented rescue parser).

#### Scenario: Parquet file parity
- **WHEN** the same sample is processed via process isolation and via thread-based mode
- **THEN** the resulting `{sample_id}_variants.parquet` files SHALL be identical

#### Scenario: Stats CSV parity
- **WHEN** the same manifest is processed via both modes
- **THEN** all aggregate CSVs (sample_summary, tier_summary, dataset_summary, etc.) SHALL be identical

### Requirement: malloc_trim called after each sample
After writing the per-sample parquet file and freeing the DataFrame, the worker SHALL call `malloc_trim(0)` via `ctypes` to release freed glibc pages back to the kernel.

#### Scenario: malloc_trim after successful sample
- **WHEN** a sample is successfully processed and written to parquet
- **THEN** `malloc_trim(0)` SHALL be called after `del df` and `gc.collect()`

#### Scenario: malloc_trim skipped if unavailable
- **WHEN** `ctypes.CDLL("libc.so.6")` raises an exception (non-glibc platform)
- **THEN** the worker SHALL silently skip the `malloc_trim` call and continue
