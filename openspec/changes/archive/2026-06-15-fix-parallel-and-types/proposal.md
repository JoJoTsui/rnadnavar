## Why

The Rust implementation (`stats_core`) was built iteratively against pre-1.0 noodles APIs. It works but has gaps: the VCF parser returns all-strings (breaking downstream statistics), the BAM module only does whole-genome stats (no per-position pileup), the VCF parser double-reads files (1.4× speedup vs expected 5-10×), and there's no Arrow interop (Python dicts instead of zero-copy). This proposal covers a thorough review and fix of all Rust-side issues plus the Python integration layer.

## What Changes

### 1. Rust VCF Parser — Fix + Optimize
- **Fix type casting**: `_cast_columns()` in `rust_vcf.py` to match `rescue_parser.py` output schema (Int64/Float64/Boolean/String columns). This fixes 6 downstream type errors across `statistics.py`, `tiering_stats.py`, `visualizer.py`, `cli.py`.
- **Optimize single-pass**: Currently reads the file twice (header dedup + records). Pass the bgzf reader through or write deduplicated header as a temp file to avoid double I/O.
- **Handle edge cases**: Missing INFO fields, empty VCFs, multi-alt alleles (only first ALT currently extracted), records with no POS.
- **Verify output parity**: Run a comparison script that parses the same VCF with both Rust and Python, then diffs every column.

### 2. Rust BAM Module — Audit + Document Gaps
- **Whole-genome stats**: `mean_coverage` is always 0 (needs reference genome lengths from BAM header). Fix by reading `bam.header().reference_sequences()`.
- **Missing features**: Per-position pileup with strand bias + base quality blocked by noodles-sam version conflict. Document the dependency chain and version requirements. The Python pysam fallback (`rust_bam.py::_pileup_pysam`) handles this correctly and is the recommended path until noodles stabilizes.
- **BAM index validation**: Warn on stale `.bai` files (index older than BAM). Currently pysam emits a C-level warning; the Rust reader should surface this more cleanly.

### 3. Python Wrapper Layer — Complete Integration
- **`rust_vcf.py`**: Add type casting, verify fallback works (stats_core absent → cyvcf2), add logging for which parser is active.
- **`rust_bam.py`**: `pileup_variants()` pysam implementation works but is single-threaded. Add ThreadPoolExecutor for per-position parallelism. Remove the `stats_core` pileup path until Rust implementation exists.
- **CLI integration**: `process_single_sample()` currently imports directly from `rescue_parser`. Change to import from `rust_vcf` (already done). Add `--parser` flag to switch between rust/python for debugging.

### 4. Build System — Make It Reproducible
- **Environment**: Document `PYO3_PYTHON` and `CONDA_PREFIX` requirements in a Makefile or script.
- **Version pinning**: Current Cargo.toml pins work. Add a `cargo update --dry-run` check to detect breakage.
- **CI readiness**: Add a `make build-rust` target that handles the maturin develop flow.

### 5. Performance — Benchmark + Profile
- **Current baseline**: VCF 1.4× (13K rec/s vs 9K), BAM 1.4M reads/s. Profile the VCF parser to identify the bottleneck (double-read, Python dict conversion, or noodles-vcf parsing itself).
- **Expected after optimization**: VCF 5-10× with single-pass + Arrow interop.
- **Memory**: 5.8M variants × ~100 columns = ~3 GB. Verify no leaks in the Rust parser.

## Capabilities

### Modified Capabilities
- `rust-vcf-parsing`: Fix type casting, optimize single-pass, edge case handling, output parity
- `rust-bam-pileup`: Fix coverage computation, document gaps, improve BAM index handling
- `variant-visualization`: Fix N_SUPPORT_CALLERS ordinal sorting after type fix
- `variant-tiering-stats`: Fix type comparisons after column casting

## Impact

- `bin/vcf_stats/seq2neo/rust_vcf.py` — add `_cast_columns()`, logging
- `bin/vcf_stats/seq2neo/rust_bam.py` — add ThreadPoolExecutor for pileup, document gaps
- `bin/vcf_stats/seq2neo/stats_core/src/vcf.rs` — optimize single-pass, fix edge cases
- `bin/vcf_stats/seq2neo/stats_core/src/bam.rs` — fix coverage, improve index handling
- `bin/vcf_stats/seq2neo/stats_core/src/lib.rs` — add coverage function export
- `bin/vcf_stats/seq2neo/cli.py` — add `--parser` flag
- `bin/vcf_stats/seq2neo/statistics.py` — defensive type casts where helpful
- `bin/vcf_stats/seq2neo/visualizer.py` — fix ordinal scale after type fix
