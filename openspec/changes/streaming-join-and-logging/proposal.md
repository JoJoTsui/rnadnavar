## Why

Samples with 7M variants still OOM-kill at 8-10 samples despite auto-throttling. Analysis shows three compounding memory issues: (1) the rescue DataFrame carries 63 output-only columns through processing (~10 GB dead weight), (2) all 6 callers' column data accumulates in memory before joining (~5.5 GB), and (3) polars join operations create new DataFrames without freeing old ones between iterations (~10-15 GB stale DataFrames). Additionally, there's no memory instrumentation to identify which step causes the kill.

## What Changes

- **Streaming join**: Join each caller immediately after parsing, free the caller's column data, then parse the next caller. Eliminates the 6-caller accumulation.
- **Column pruning**: Split rescue_df into a slim processing DataFrame (12 columns) and an output-only DataFrame (63 columns). Only the slim df participates in the join; output columns are hstack'd back before writing parquet.
- **Drop 6 unnecessary columns**: Remove N_CALLERS, RaVeX_FILTER, N_DNA_CALLERS, N_RNA_CALLERS, DNA_SUPPORT, RNA_SUPPORT from the rescue parser.
- **Memory logging**: Add `_mem()` instrumentation at every step with immediate flush. The last printed line before a kill identifies the failing step.
- **Log level separation**: Critical progress shown on terminal; caller-level detail hidden behind `--verbose`.

## Capabilities

### New Capabilities
- `streaming-caller-join`: Join callers one-at-a-time, freeing intermediates between iterations
- `memory-instrumentation`: RSS logging at every processing step with immediate flush

### Modified Capabilities
- None — implementation change only

## Impact

- **Python code**: `cli.py` (streaming join + column pruning + memory logging), `caller_parser.py` (per-caller join function), `rescue_parser.py` (drop 6 fields)
- **Rust code**: None
- **Tests**: Memory logging tests, streaming join tests, column pruning tests
