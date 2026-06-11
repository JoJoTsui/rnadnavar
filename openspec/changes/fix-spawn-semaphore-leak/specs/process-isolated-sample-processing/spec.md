## MODIFIED Requirements

### Requirement: Each sample processed in an independent Python process
The pipeline SHALL process each sample in a separate Python process using
`multiprocessing.Pool` with `maxtasksperchild=1` and `context='fork'`. When a
worker process completes its sample, it SHALL exit, and the operating system
SHALL reclaim all memory allocated by that process. No `resource_tracker`
semaphore warnings SHALL be emitted at shutdown.

#### Scenario: Memory reclaimed after sample completes
- **WHEN** a worker process finishes processing a sample and exits
- **THEN** the process RSS attributable to that sample SHALL return to the operating system
- **AND** subsequent samples processed by new workers SHALL start with a clean memory slate

#### Scenario: Concurrent samples in separate processes
- **WHEN** 4 worker processes are processing 4 samples concurrently
- **THEN** each process SHALL have its own independent memory space
- **AND** memory freed in one process SHALL NOT affect memory available to other processes

#### Scenario: No semaphore leak warning
- **WHEN** the pipeline completes all samples and the Pool is properly closed
- **THEN** the `resource_tracker` SHALL NOT emit "leaked semaphore objects" warnings

#### Scenario: Sequential mode fallback
- **WHEN** `--sample-workers 1` is specified
- **THEN** samples SHALL be processed sequentially in the main process
- **AND** no subprocess SHALL be spawned

#### Scenario: Fork context is used
- **WHEN** `--sample-workers > 1` and `--process-mode` is not explicitly set to `thread`
- **THEN** the pool SHALL be created with `multiprocessing.get_context("fork")`
- **AND** worker processes SHALL inherit the parent's module state via fork
