## MODIFIED Requirements

### Requirement: Each sample processed in an independent Python process
The pipeline SHALL process each sample in a separate Python process using
`multiprocessing.Pool` with `maxtasksperchild=1` and `context='spawn'`. When a
worker process completes its sample, it SHALL exit, and the operating system
SHALL reclaim all memory allocated by that process. The `resource_tracker`
warning about leaked semaphore objects SHALL be suppressed via `warnings.filterwarnings`.

#### Scenario: Memory reclaimed after sample completes
- **WHEN** a worker process finishes processing a sample and exits
- **THEN** the process RSS attributable to that sample SHALL return to the operating system
- **AND** subsequent samples processed by new workers SHALL start with a clean memory slate

#### Scenario: Spawn context is used
- **WHEN** `--sample-workers > 1` and `--process-mode` is not explicitly set to `thread`
- **THEN** the pool SHALL be created with `multiprocessing.get_context("spawn")`
- **AND** worker processes SHALL be fresh Python interpreters (fork+exec)

#### Scenario: resource_tracker warning suppressed
- **WHEN** the pipeline starts
- **THEN** `warnings.filterwarnings('ignore', message='resource_tracker')` SHALL be called
- **AND** no "leaked semaphore objects" warning SHALL appear at shutdown
