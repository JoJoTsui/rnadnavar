# stats_core — Rust backend for seq2neo statistics

Rust implementation of VCF parsing, BAM statistics, variant tiering, and BAM
pileup for the rnadnavar seq2neo pipeline. Exposed to Python via PyO3.

## Requirements

| Tool | Minimum Version | Notes |
|------|----------------|-------|
| Rust (cargo, rustc) | 1.70+ | https://rustup.rs |
| maturin | 1.0+ | `pip install maturin` |
| Python | 3.10+ | Target environment must be active |

## Quick Start

```bash
# From the stats_core directory:
./build_rust.sh

# Or directly with maturin:
maturin develop --release
```

The built module is installed into the active Python environment as `stats_core`.

## Build Modes

| Mode | Command | Compile Time | Runtime Speed | Use Case |
|------|---------|-------------|---------------|----------|
| **release** | `./build_rust.sh` or `maturin develop --release` | Slow (~2-5 min) | Fast (optimized) | Production, benchmarking, full pipeline runs |
| **debug** | `./build_rust.sh debug` or `maturin develop` | Fast (~30-60s) | Slow (unoptimized) | Development, debugging, quick iteration |

**Recommendation:** Use **release** for any pipeline run with >1 sample. The
runtime difference is 5-10× for VCF parsing and BAM scanning. Use **debug**
only when iterating on Rust code and recompiling frequently.

## Project Structure

```
stats_core/
├── Cargo.toml          # Rust dependencies (noodles, pyo3)
├── build_rust.sh       # Build script wrapper
├── README.md           # This file
└── src/
    ├── lib.rs          # PyO3 module exports
    ├── vcf.rs          # VCF parser (rescue VCF → column-oriented dicts)
    ├── bam.rs          # BAM statistics (reads, coverage, insert size)
    ├── pileup.rs       # Per-position BAM pileup for variant validation
    ├── caller.rs       # Caller VCF parsing (Mutect2, Strelka2, SAGE, DeepSomatic)
    └── tier.rs         # Variant tiering (evidence-based classification)
```

## Key Dependencies

| Crate | Version | Purpose |
|-------|---------|---------|
| pyo3 | 0.28 | Python bindings |
| noodles-vcf | 0.88 | VCF parsing |
| noodles-bam | 0.90 | BAM record iteration (re-exports noodles-sam types) |
| noodles-bgzf | 0.36 | BGZF block reading |
| noodles-csi | 0.42 | BAM index (BAI) reading |
| noodles-core | 0.20 | Core types (Position, Region) |
| pyo3-arrow | 0.17 | Arrow/Polars interop (future: zero-copy) |

## Version Compatibility

All noodles crates are pinned to the latest available versions on crates.io
(as of 2026-06). There are no known version conflicts — `noodles-bam 0.90`,
`noodles-vcf 0.88`, `noodles-bgzf 0.36`, `noodles-csi 0.42`,
and `noodles-core 0.20` work together. `noodles-sam` is not a direct
dependency — its alignment types (CIGAR ops, flags) are re-exported by
`noodles-bam`.

### Rust vs Python

The Rust module handles all BAM operations:
- **Whole-genome/WES statistics** (reads, coverage, insert size, mapping quality)
  via `bam.rs` — exposed as `stats_core.bam_stats()` and `stats_core.bam_stats_bed()`
- **Per-position pileup** (DP, REF_DP, ALT_DP, strand bias F1R2/F2R1, mean_BQ,
  mean_MQ) via `pileup.rs` — exposed as `stats_core.pileup_variants()` and
  `stats_core.pileup_variants_multi()`

There is no Python fallback for BAM pileup. If `stats_core` is not built, the
pipeline raises `ImportError` with instructions to run `./build_rust.sh`.

## Troubleshooting

### `error: can't find crate for 'stats_core'`
The build script must be run from the `stats_core/` directory (where `Cargo.toml` lives).

### `PYO3_PYTHON` not set
maturin auto-detects the active Python. If you see errors, explicitly set it:
```bash
export PYO3_PYTHON=$(which python)
maturin develop --release
```

### Linker errors about missing Python symbols
Ensure the Python environment has shared libraries:
```bash
# Conda/mamba:
conda install -c conda-forge python

# Or set:
export LD_LIBRARY_PATH=$(dirname $(which python))/../lib:$LD_LIBRARY_PATH
```

### Module not found after build
maturin installs into the active Python's site-packages. Verify:
```bash
python -c "import stats_core; print(dir(stats_core))"
```
