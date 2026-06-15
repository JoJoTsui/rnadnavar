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
| noodles-bam | 0.90 | BAM record iteration |
| noodles-sam | 0.85 | SAM/BAM alignment types (CIGAR ops, flags) |
| noodles-bgzf | 0.36 | BGZF block reading |
| noodles-csi | 0.42 | BAM index (BAI) reading |
| noodles-core | 0.20 | Core types (Position, Region) |
| pyo3-arrow | 0.17 | Arrow/Polars interop (future: zero-copy) |

## Version Compatibility Notes

### noodles-sam (0.85) and per-position pileup

The `noodles-sam` crate provides alignment record types including CIGAR
operations (`Kind::consumes_reference()`) and flags. However, the current
noodles ecosystem has a version conflict that blocks implementing per-position
pileup (strand bias, base quality) entirely in Rust:

- `noodles-bam 0.90` depends on `noodles-sam ~0.85`
- `noodles-csi 0.42` depends on `noodles-core 0.20`
- Full per-position pileup (querying individual bases, strand, quality) would
  require `noodles-sam >= 0.90` for the complete alignment record API

**Current state:** The Rust BAM module does whole-genome/WES statistics (reads,
coverage, insert size, mapping quality) via `noodles-bam 0.90`. The per-position
pileup with strand bias and base quality is handled by the Python pysam fallback
in `rust_bam.py::_pileup_pysam()`. This is the recommended path until noodles
stabilizes its SAM/BAM type hierarchy.

**Upgrade path:** When noodles releases aligned versions (estimated noodles 0.120+),
the pileup can be migrated to Rust for a 5-10× speedup. Track:
- https://crates.io/crates/noodles-sam
- https://crates.io/crates/noodles-bam

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
