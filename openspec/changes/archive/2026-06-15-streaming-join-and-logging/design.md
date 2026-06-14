## Context

The pipeline processes 65 samples, some with 7M variants. Each sample's rescue VCF has 191 columns (75 INFO fields + 105 caller columns + derived + meta). During processing, only 12 columns are actively used; the remaining 63 output-only columns consume 10 GB of dead weight for a 7M-variant sample. Additionally, all 6 callers' column data (Python lists) accumulates in a dict before joining, peaking at 5.5 GB. The polars join loop creates new DataFrames without GC between iterations, leaving 10-15 GB of stale DataFrames in Rust memory.

## Goals / Non-Goals

**Goals:**
- Reduce per-sample peak memory from ~28 GB to ~17 GB
- Pinpoint which step causes OOM kills via memory logging
- Clean up terminal output: critical progress only, detail on `--verbose`

**Non-Goals:**
- Changing Rust functions (no rebuild needed)
- Reducing the final joined DataFrame size (12 GB is the data, unavoidable)

## Decisions

### 1. Streaming join — one caller at a time

**Decision:** Replace `parse_all_callers` + `join_caller_columns` with a streaming loop in `process_single_sample`:

```python
df = rescue_slim
for caller_name, cfg in CALLER_CONFIGS.items():
    cols = _parse_one_caller(...)           # returns column-oriented dict
    df = _join_one_caller(df, cols, caller_name)  # join + rename + free
    del cols
    gc.collect()
```

Each caller's 3.5 GB of Python lists is freed immediately after joining. The largest single caller (DNA_deepsomatic at 6.5M rows) sets the peak, not the sum of all 6.

### 2. Column pruning

**Decision:** After rescue parsing, split into slim (12 cols) and output (63 cols):

```python
SLIM_COLS = ["CHROM", "POS", "REF", "ALT", "FILTER", "FILTERS_NORMALIZED",
             "GNOMAD_AF", "COSMIC_CNT", "REDI_EVIDENCE",
             "N_DNA_CALLERS_SUPPORT", "N_RNA_CALLERS_SUPPORT"]

rescue_slim = rescue_df.select(SLIM_COLS)
rescue_output = rescue_df.select([c for c in rescue_df.columns if c not in SLIM_COLS])
del rescue_df
```

After all processing: `df = df.hstack(rescue_output)` before writing parquet.

### 3. Drop 6 columns permanently

**Decision:** Remove from `RESCUE_STRING_FIELDS` in `rescue_parser.py`:
- N_CALLERS
- RaVeX_FILTER
- N_DNA_CALLERS
- N_RNA_CALLERS
- DNA_SUPPORT  
- RNA_SUPPORT

These are redundant with other fields or not used.

### 4. Memory logging with immediate flush

**Decision:** Add `_mem(step_name)` function using `psutil`:

```python
import psutil

def _mem(msg):
    rss = psutil.Process().memory_info().rss / (1024**3)
    print(f"  [MEM {rss:.1f}GB] {msg}", flush=True)
```

Called at every step boundary. The last printed `[MEM X.XGB]` before a kill identifies the failing step.

Log points:
- `_mem(f"start [{sample_id}]")`
- `_mem("after rescue parse")`
- `_mem("after column split")`  
- `_mem("after caller N: {caller_name}")` (for each of 6 callers)
- `_mem("after all callers joined")`
- `_mem("after stats")`
- `_mem("after parquet write — free")`

### 5. Log level separation

**Decision:** Default terminal output is concise. `--verbose` shows caller-level detail.

```
Default:                                   --verbose:
  [1/65] PRJNA298376_3812  OK  (0.5 GB)      [DNA_mutect2] scanning ...mutect2.dec.norm.vcf.gz
  [2/65] PRJNA298376_3942  OK  (1.2 GB)      [DNA_mutect2] found 5576 variants (Rust)
  [MEM 2.1GB] after rescue parse             [RNA_mutect2] scanning ...mutect2.dec.norm.vcf.gz
  [MEM 2.3GB] after column split             [RNA_mutect2] found 10220 variants (Rust)
  ...                                         ...
  [MEM 8.4GB] after join #4: DNA_strelka     [MEM 2.1GB] after rescue parse
  [MEM 9.1GB] after join #5: RNA_strelka     [MEM 2.3GB] after column split
  [MEM 9.1GB] after stats                    ...
  [MEM 0.2GB] after free — done
```

## Memory Budget (7M-variant sample, after fixes)

```
Step                              Memory
──────────────────────────────────────────
after rescue parse                 10.7 GB  (191 cols)
after column split                  9.6 GB  (slim 0.7 + output 8.9 idle)
caller 1 parsed (DeepSomatic)      13.1 GB  (+3.5 caller cols)
after caller 1 join+free           10.2 GB  (df 1.6 + output 8.9, cols freed)
caller 2 parsed (DNA_strelka)      11.2 GB  (+1.0)
after caller 2 join+free           11.1 GB  (df 2.2 + output 8.9)
caller 3-5 (smaller)               11-12 GB
after caller 6 join+free           15.1 GB  (df 6.2 + output 8.9)
after stats                        15.7 GB
hstack output + write + free        0.2 GB  (only all_stats remains)
──────────────────────────────────────────
Peak:                              ~16 GB   (was ~28 GB)
```

## Risks / Trade-offs

- **[Risk] `psutil` not available**: → **Mitigation:** Fall back to `/proc/self/status` VmRSS, no dependency needed.
- **[Risk] Streaming join is slower**: Sequential caller parse+join takes longer than parallel parse. → **Mitigation:** Rust callers are fast (4s total), join is the bottleneck either way. Acceptable.
- **[Risk] Column list out of sync with rescue parser**: New INFO fields added to rescue_parser.py but not to SLIM_COLS would be lost. → **Mitigation:** Test verifies all columns present in output parquet.
