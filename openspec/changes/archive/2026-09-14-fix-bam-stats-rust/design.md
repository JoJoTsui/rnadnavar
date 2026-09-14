## Context

The BAM statistics subsystem (`bam_stats.py` 724 lines, `rust_bam.py` 175 lines, `stats_core/src/bam.rs`, `stats_core/src/tier.rs`) computes per-sample alignment metrics and per-variant tier assignments via a compiled Rust extension (`stats_core.so`). Four bugs produce wrong or missing data:

1. **Expanded metrics dropped**: `bam_stats.py:474-481` (`_compute_bam_stats_rust`) builds the return dict with only 6 base metrics (total_reads, mapped_reads, mapping_rate_pct, mean_coverage, mean_insert_size, mean_mapq). The Rust backend (`bam.rs:200-217`) computes and returns 3 additional metrics (duplication_rate_pct, properly_paired_pct, insert_size_stddev), but the Python wrapper never reads them. `bam_stats.py:532-534` then does `result["duplication_rate_pct"] = result.get("duplication_rate_pct")` which returns `None`. Result: 7 of 20 columns in `bam_stats.tsv` are always null.

2. **MAPQ=255 artifact**: `bam.rs:131-133` — `if let Some(mq) = record.mapping_quality()` — noodles represents MAPQ=255 (SAM-spec "not available") as `None`. STAR sets 255 for ~99% of uniquely-mapped RNA reads. These reads contribute 0 to `mq_sum` but are still counted in `mapped` (incremented at `bam.rs:129`). `mean_mapq = mq_sum / mapped` ≈ 0.2 for ALL 60+ RT samples. DNA samples (BWA, real MAPQ 0-60) show 54-59. The pysam fallback would give ~254 (includes 255 literally) — also wrong.

3. **coverage_bins off-by-one**: `bam.rs:247-257` — BED is 0-based half-open `[start, end)`; noodles `Position` is 1-based inclusive. The code passes `Position(start)` and `Position(end)` with no `+1` conversion, shifting coverage one base left. Also, `NonZero::new(0)` returns `None` → regions starting at 0 are skipped via `continue`, but `total_bases` was already incremented → `cov_*_pct` deflated. Combined with bare `except Exception` at `bam_stats.py:537-553`, all `cov_*` columns are null with no diagnosability.

4. **Tier fallback misassigns RNAedit**: `tier.rs:199-202` — When concordant callers = (0,0), Rust falls back to raw `N_DNA/RNA_CALLERS_SUPPORT`. RNAedit variants have 0 concordant callers (RNAedit is assigned by REDIportal, not by any caller voting "RNAedit"). So 598,343 RNAedit variants are reassigned to C3D1 (≥2 raw RNA callers) instead of C7D1 (0 concordant + database). The C7 tier is completely absent from all output. The Python path (`tiering_stats.py:199`) correctly returns C7 for (0,0) concordant — it only falls back on exception, which never fires.

The pysam fallback in `bam_stats.py` is dead code — the user confirmed only the Rust backend is used. Removing it eliminates ~100 lines of maintenance burden and the backend-divergence risk (pysam and Rust disagree on MAPQ, insert-size capping, stddev divisor, duplication rate sampling, properly_paired definition, and coverage scan behavior).

## Goals / Non-Goals

**Goals:**
- Populate the 3 expanded metrics that Rust already computes
- Fix MAPQ=255 to produce correct RNA mean_mapq values
- Fix coverage_bins coordinate conversion and start=0 handling
- Fix Rust tier fallback to return C7 for (0,0) concordant, matching Python
- Remove pysam fallback code (Rust-only)

**Non-Goals:**
- Fixing visualization bugs (that belongs in `fix-broken-visualizations`)
- Fixing multi-allelic classification (that belongs in `fix-multiallelic-logic`)
- Adding new BAM metrics beyond the 3 already computed by Rust
- Changing the BAM file discovery logic (`_locate_bam_file`)
- Rebuilding `stats_core.so` — this proposal specifies the Rust source changes; the build step is an implementation task

## Decisions

### Decision 1: Read all Rust-returned metrics in the Python wrapper

**Chosen**: Add `duplication_rate_pct`, `properly_paired_pct`, `insert_size_stddev` to the return dict in `_compute_bam_stats_rust` by reading from `raw["..."]`. These are already computed by `bam.rs:200-217` and returned through `lib.rs:121-123`. One line per metric.

**Alternatives considered**:
- Re-compute in Python → Rejected: duplicates Rust logic, slower, introduces inconsistency
- Add a generic "copy all keys" loop → Rejected: less explicit, could accidentally pass through debug keys

### Decision 2: Exclude MAPQ=255 from both numerator and denominator

**Chosen**: In `bam.rs:129-133`, only increment `mapped` when `mapping_quality()` is `Some`. Use a separate `mapq_count` variable for the denominator:
```rust
if let Some(mq) = record.mapping_quality() {
    mq_sum += u8::from(mq) as f64;
    mapq_count += 1;
}
```
Then `mean_mapq = mq_sum / mapq_count` (not `mq_sum / mapped`). Reads with MAPQ=255 are excluded from both. `mapped` still counts all non-unmapped reads for `mapping_rate_pct`.

**Alternatives considered**:
- Include 255 in the numerator as 255 → Rejected: inflates mean to ~254, meaningless
- Set MAPQ=255 to 0 → Rejected: deflates mean, same problem as current but less extreme
- Document as "mean MAPQ over reads with explicit MAPQ" → This IS what the chosen approach implements

### Decision 3: BED 0-based → noodles 1-based with +1 conversion

**Chosen**: In `bam.rs:247-257`, convert `start` and `end` to 1-based by adding 1:
```rust
let pos_start = (start + 1) as usize;  // 0-based BED → 1-based noodles
let pos_end = (end) as usize;          // BED end is exclusive → 1-based inclusive = end
```
For `start=0`: `pos_start = 1` (valid NonZero). No more `NonZero::new(0)` skip. The depth array indexing at `bam.rs:284` also needs adjustment: `base_pos = pos + offset - (start + 1) as i64`.

**Alternatives considered**:
- Use 0-based noodles API → Rejected: noodles doesn't have a 0-based API; Position is inherently 1-based
- Subtract 1 from start instead → Rejected: doesn't fix the start=0 case

### Decision 4: Remove the raw-caller-count tier fallback

**Chosen**: In `tier.rs:199-202`, replace the fallback with a direct C7 return:
```rust
if dna_count == 0 && rna_count == 0 {
    // No concordant callers — return C7 tier
    // (RNAedit, unclassified, or annotation-only variants)
    return TierResult { caller_tier: "C7", database_tier: ..., ... };
}
```
Do NOT fall back to raw `N_DNA/RNA_CALLERS_SUPPORT`. The Python path already does this correctly.

**Alternatives considered**:
- Keep fallback but exclude RNAedit → Rejected: requires knowing the FILTER category in Rust, which is available but adds complexity; the Python path doesn't need the fallback at all
- Move all tier computation to Python → Rejected: Rust is faster for the common case; only the (0,0) fallback is wrong

### Decision 5: Remove pysam fallback entirely

**Chosen**: Remove all pysam-specific functions in `bam_stats.py`: `_compute_duplication_rate` (206-216), `_compute_properly_paired_pct` (240-256), `_compute_insert_size_stddev` (280-292), `_compute_coverage_bins` (322-343), and the pysam branch of `compute_bam_stats` (387-439). Remove the `HAS_RUST_BAM` check at `bam_stats.py:134` — Rust is required, not optional. Add a startup check: `if not hasattr(stats_core, 'bam_stats'): raise RuntimeError("stats_core.bam_stats not found — Rust extension required")`.

**Alternatives considered**:
- Keep as documented fallback → Rejected: user confirmed Rust-only; the fallback has 6 backend-divergence issues that would all need fixing
- Keep but never call → Rejected: dead code with maintenance burden

## Risks / Trade-offs

- **[Rust recompilation required]**: Changes to `bam.rs` and `tier.rs` require rebuilding `stats_core.so`. The build process must be documented or automated. → **Mitigation**: Add a build script or Makefile target for `stats_core.so`; document the Rust toolchain version.
- **[MAPQ value change]**: RNA `mean_mapq` will change from ~0.2 to the correct value (likely 20-40 for STAR multi-mappers, or null if all reads are MAPQ=255). Users comparing old and new runs will see different values. → **Mitigation**: Document the change; the old value was wrong.
- **[Tier distribution change]**: 598K RNAedit variants will move from C3D1 to C7D1. All downstream stats, charts, and summaries will show a new C7D1 tier. → **Mitigation**: This is the correct behavior; C7D1 was supposed to exist per `tier_config.py:206`.
- **[Coverage values change]**: Fixed coordinate conversion will produce slightly different `cov_*_pct` values. → **Mitigation**: The old values were wrong (shifted by one base); the new values are correct.
- **[Hard dependency on Rust]**: Removing pysam fallback means the pipeline cannot run without `stats_core.so`. → **Mitigation**: Add a clear startup error; the Rust extension is already required for VCF parsing and pileup.
