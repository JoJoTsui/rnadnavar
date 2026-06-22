## Context

The `expand-stats-visualizations` proposal (P2) added 4 expanded BAM metrics to `bam_stats.py` using pysam functions that scan the entire BAM file. Each function (`_compute_duplication_rate`, `_compute_properly_paired_pct`, `_compute_insert_size_stddev`, `_compute_coverage_bins`) opens the BAM independently and iterates all reads. With 183 BAM files and 8 workers, this takes ~23 hours. The Rust `stats_core` backend already performs a single-pass BAM scan in ~30 seconds for the 6 original metrics — these 4 new metrics can be computed in the same single pass with zero additional I/O.

The `stats_core` crate uses `noodles-bam` 0.90 which exposes all needed SAM flags: `is_duplicate()` (0x400), `is_properly_segmented()` (0x2), and `template_length()` for insert size. All three are already accessed in the existing loop — they just need counters.

## Goals / Non-Goals

**Goals:**
- Compute duplication rate, properly paired pct, and insert size stddev in the existing Rust BAM scan (single pass, zero extra I/O)
- Compute coverage bins via a new Rust function using indexed BAM queries per BED region
- Remove pysam post-processing from the Rust path in `compute_bam_stats()`
- Add timeout safety net to `bam_future.result()`
- Add `ensure_bam_stats_columns()` to `--resume` path
- Keep pysam fallback functions for environments without Rust, with `max_reads` cap

**Non-Goals:**
- Changing the output column schema — same 8 column names, same semantics
- Changing the `_CROSS_SAMPLE_COLS` or variant-level parquet — BAM stats are sample-level
- Optimizing the pysam fallback path beyond adding max_reads
- Adding coverage bins to WGS mode (requires reference genome FASTA for chromosome lengths)

## Decisions

### Decision 1: 3 counters in existing BamStats struct, not a separate function

**Chosen**: Add `duplicate_count`, `proper_pair_count`, `insert_size_sum`, `insert_size_sum_sq` fields to the existing `BamStats` struct and increment them in the existing `for result in reader.records()` loop.

**Alternatives considered**:
- Separate Rust function per metric → Rejected: each would open and scan the BAM independently, same problem as pysam
- Compute in Python from the Rust record iterator → Rejected: noodles-bam records aren't exposed to Python

**Rationale**: The existing loop already checks `flags.is_unmapped()`, `flags.is_properly_segmented()`, `record.mapping_quality()`, and `record.template_length()`. Adding 3 boolean flag checks per record is ~nanoseconds of overhead. The insert size stddev uses Welford's online algorithm (sum and sum-of-squares) to avoid storing all values.

### Decision 2: coverage_bins as a separate Rust function

**Chosen**: New `coverage_bins(bam_path, bed_regions)` function using noodles-bam indexed queries per BED region, counting per-base depth.

**Alternatives considered**:
- Add to the existing pileup infrastructure → Rejected: pileup is position-specific, coverage binning needs region-wide per-base depth
- Compute in the same single-pass scan → Rejected: per-base depth tracking across the entire genome (or BED regions) requires storing a coverage array, which doesn't fit in the existing stat-accumulation pattern

**Rationale**: Coverage binning is inherently a separate pass — it counts per-base depth, not per-read statistics. But it uses indexed queries per BED region, not full-BAM scanning, so it's O(BED size) not O(BAM size).

### Decision 3: Rust implementation + pysam fallback coexist

**Chosen**: Implement in Rust, keep pysam functions for fallback, add `max_reads` cap to pysam fallbacks.

**Alternatives considered**:
- Remove pysam functions entirely → Rejected: need fallback for environments without `stats_core` compiled
- Make pysam the default and Rust optional → Rejected: defeats the purpose

**Rationale**: The Rust path is the production path. The pysam path is a development/fallback path. Adding `max_reads=10_000_000` ensures the fallback never takes more than ~30 seconds per BAM file.

### Decision 4: Timeout on future.result()

**Chosen**: `bam_future.result(timeout=3600)` — 1 hour timeout.

**Rationale**: Even with Rust, a very large BAM could take minutes. 1 hour is conservative for 183 BAM files with 8 workers. If BAM stats exceed 1 hour, something is wrong and it's better to fail fast with a warning than hang indefinitely.

## Risks / Trade-offs

- **[Rust build]** — `maturin develop --release` requires a Rust toolchain. If the build fails, the pysam fallback is still available. → **Mitigation**: The pysam fallback functions remain intact.
- **[coverage_bins not available for WGS]** — The new function requires BED regions. For WGS mode, coverage bins remain null. → **Mitigation**: This matches the existing behavior — coverage bins were always BED-only.
- **[Insert size stddev precision]** — Welford's online algorithm can accumulate floating-point error for large datasets. → **Mitigation**: Use f64 throughout; error is < 0.01 for realistic BAM sizes.
