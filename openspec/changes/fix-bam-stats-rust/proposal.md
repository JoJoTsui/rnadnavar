## Why

The BAM statistics subsystem has 4 critical bugs that produce wrong or missing metrics: (1) the Rust backend computes expanded metrics (duplication_rate_pct, properly_paired_pct, insert_size_stddev) but the Python wrapper silently drops them — 7 of 20 columns are always null; (2) RNA `mean_mapq ≈ 0.2` for all 60+ RT samples because the Rust code excludes MAPQ=255 (STAR's "not available" sentinel) from the numerator but keeps it in the denominator; (3) `coverage_bins` has a BED 0-based ↔ noodles 1-based off-by-one error plus silently skips regions starting at position 0; (4) the Rust tier fallback misassigns 598K RNAedit variants to C3D1 instead of C7D1, making the C7 tier completely absent from all output. The pysam fallback is dead code per the user's direction (Rust-only).

## What Changes

- **Fix expanded metrics drop**: Read `duplication_rate_pct`, `properly_paired_pct`, `insert_size_stddev` from the Rust backend's return dict in `_compute_bam_stats_rust` at `bam_stats.py:474-481`. Currently the Python wrapper builds the return dict with only 6 base metrics and never reads the 3 expanded fields that Rust already computes.
- **Fix MAPQ=255 handling in Rust**: Exclude MAPQ=255 reads from BOTH the numerator AND denominator of `mean_mapq` in `stats_core/src/bam.rs:129-133`. STAR sets MAPQ=255 for ~99% of uniquely-mapped RNA reads; noodles represents 255 as `None`, so the current code adds 0 to `mq_sum` but still increments `mapped`. Fix: only increment `mapped` when `mapping_quality()` is `Some`, and use a separate `mapq_count` for the denominator.
- **Fix coverage_bins off-by-one**: Convert BED 0-based coordinates to noodles 1-based by adding `+1` to both start and end positions in `stats_core/src/bam.rs:247-257`. Fix the start=0 skip by handling `start=0` as a valid BED interval (convert to 1-based `start=1`). Replace the bare `except Exception` at `bam_stats.py:537-553` with targeted error handling and logging.
- **Fix Rust tier fallback for RNAedit**: In `stats_core/src/tier.rs:199-202`, remove the raw-caller-count fallback when concordant=(0,0). Instead, return C7 (0 concordant callers) as the Python path does. The fallback to raw `N_DNA/RNA_CALLERS_SUPPORT` was intended for cases where FILTERS_NORMALIZED is empty, but it misassigns RNAedit variants (which have 0 concordant callers by design — RNAedit is assigned by REDIportal annotation, not by callers).
- **Remove pysam fallback code**: Per the user's direction, the pysam backend is unused. Remove all pysam-specific computation paths in `bam_stats.py` (`_compute_duplication_rate`, `_compute_properly_paired_pct`, `_compute_insert_size_stddev`, `_compute_coverage_bins`, the pysam branch of `compute_bam_stats`). Keep only the Rust path. Add a startup check that verifies `stats_core` loads.

## Capabilities

### New Capabilities

- `rust-bam-stats-correction`: Corrected Rust BAM statistics computation — MAPQ=255 exclusion, expanded metrics pass-through, coverage_bins coordinate conversion, and tier fallback fix.

### Modified Capabilities

- `bam-statistics`: Expanded metrics (duplication_rate_pct, properly_paired_pct, insert_size_stddev) SHALL be populated from the Rust backend instead of being null-filled. MAPQ=255 reads SHALL be excluded from mean_mapq computation. Coverage_bins SHALL use correct 1-based coordinates.
- `variant-tiering-stats`: Rust tier engine SHALL return C7 for (0,0) concordant callers without falling back to raw caller counts, matching the Python path behavior.

## Impact

- `bin/vcf_stats/seq2neo/bam_stats.py` — ~100 lines: read 3 expanded metrics from Rust dict (~10 lines), remove pysam fallback paths (~70 lines), replace bare except with targeted logging (~10 lines), remove pysam-specific compute functions (~10 lines)
- `bin/vcf_stats/seq2neo/rust_bam.py` — ~5 lines: remove pysam-related fallback logic
- `bin/vcf_stats/seq2neo/stats_core/src/bam.rs` — ~30 lines: fix MAPQ=255 handling (~10 lines), fix coverage_bins coordinate conversion (~15 lines), remove pysam-inconsistent metric definitions (~5 lines)
- `bin/vcf_stats/seq2neo/stats_core/src/tier.rs` — ~10 lines: remove raw-caller-count fallback at line 199-202, return C7 for (0,0) concordant
- `bin/vcf_stats/seq2neo/stats_core/src/lib.rs` — ~5 lines: update if tier function signature changes
- **BREAKING**: `mean_mapq` for RNA samples will change from ~0.2 to the correct value (mean MAPQ over reads with explicit MAPQ). `bam_stats.tsv` will have populated expanded metric columns. Tier distribution will show C7D1 with ~599K RNAedit variants.
- Requires Rust recompilation of `stats_core.so`
