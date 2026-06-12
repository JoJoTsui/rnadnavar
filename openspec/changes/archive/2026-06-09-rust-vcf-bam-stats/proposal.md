## Why

The current Python implementation of VCF/BAM parsing and statistics is a draft with known bugs, missing analytical dimensions, and will not scale to the full dataset (3.3M+ variants across 65 samples, each with 3 BAM files). Variant-wise BAM pileup with strand/quality analysis takes ~9 hours in pysam for 65 samples. A Rust core using noodles (VCF/BAM) + maturin (Python binding) solves performance and correctness together, while a unified 4-level statistics hierarchy (dataset, per-tier, per-sample, per-sample×per-tier) ensures all analyses are consistently available at every granularity. Python handles orchestration (manifest loading, argument parsing), output generation (altair visualization, CSV/Parquet), and tiering integration — Rust handles only the compute-intensive VCF/BAM I/O.

## Technology Stack

- **Python**: orchestration, manifest processing, tiering engine, altair visualization, output writing — managed by `uv`
- **Rust** (`stats_core` crate via maturin): VCF parsing (noodles-vcf 0.88), BAM pileup (noodles-bam 0.90), stats engine — exposed as a Python extension module
- **Key dependency versions**: pyo3 0.28.3, pyo3-arrow 0.17.0, noodles 0.111, noodles-vcf 0.88, noodles-bam 0.90

## What Changes

### Rust Core — `stats_core` (maturin + noodles)
- **BREAKING**: Replace `rescue_parser.py` cyvcf2 parsing with Rust `noodles-vcf` reader — faster, handles edge cases correctly
- **BREAKING**: Replace `caller_parser.py` cyvcf2 parsing with Rust VCF reader — single-pass extraction of FORMAT fields for all 6 callers at target positions
- **BREAKING**: Replace `bam_stats.py` pysam parsing with Rust `noodles-bam` reader — per-position pileup with strand bias and base quality
- New Rust stats engine for per-variant metric computation (VAF, DNA/RNA means, REF/ALT DP)

### Python Layer (orchestration + visualization)
- Manifest loading, argument parsing, output writing stay in Python
- Tiering engine integration (`tiering_stats.py`) stays in Python
- All visualization stays in Python using altair + vl-convert-python
- `uv` manages the entire Python environment (packages, venv, maturin builds)

### Statistics Hierarchy (4 Levels) — computed in Python/polars from Rust-produced DataFrames
- **Level 1 — Dataset-wide**: All variants across all samples
- **Level 2 — Per-tier (CxDy)**: All samples grouped by C1D1..C7D0
- **Level 3 — Per-sample**: Each sample individually, plus BAM whole-genome stats
- **Level 4 — Per-sample × Per-tier**: Each sample broken by CxDy tier

### BAM Statistics (Rust)
- 3 BAM types per sample: DNA normal (DN), DNA tumor (DT), RNA tumor (RT)
- Whole-genome stats: total/mapped reads, mapping rate, coverage, insert size, MAPQ
- Variant-wise pileup: DP, REF/ALT depth, strand bias (F1R2/F2R1 per allele), base quality, mapping quality
- Two modes: `--pileup-mode all` (default) and `--pileup-mode filtered` (exclude NoConsensus)

### Visualization Updates (Python/altair)
- Existing charts updated for 4-level data
- New BAM charts: per-sample metrics bars, BAM vs caller DP/VAF scatter, coverage distribution, validation heatmap
- Dashboard organized in 4 sections: Overview, Tier Analysis, BAM & Validation, Per-Sample

## Capabilities

### New Capabilities
- `rust-vcf-parsing`: Rust-based VCF parsing for rescue and caller VCFs via noodles-vcf 0.88
- `rust-bam-pileup`: Rust-based BAM pileup at variant positions via noodles-bam 0.90
- `rust-stats-engine`: Rust per-variant metrics computation (VAF, means)
- `four-level-statistics`: Unified statistics at dataset, per-tier, per-sample, and per-sample×per-tier
- `bam-validation`: Cross-validation of BAM pileup metrics against caller VCF FORMAT fields
- `bam-visualization`: BAM statistics and validation charts (altair)

### Modified Capabilities
- `variant-visualization`: Updated dashboard with 25+ charts in 4 sections
- `bam-statistics`: Upgraded to 3 BAM types (DN/DT/RT), whole-genome + variant-wise, Rust-accelerated

## Impact

- **New**: `bin/vcf_stats/seq2neo/stats_core/` — Rust crate (~500 LOC)
- **New**: `bin/vcf_stats/seq2neo/rust_vcf.py` — Python wrapper for Rust VCF parser
- **New**: `bin/vcf_stats/seq2neo/rust_bam.py` — Python wrapper for Rust BAM pileup
- **Modified**: `bin/vcf_stats/seq2neo/statistics.py` — 4-level aggregation, `sample_tier_summary()`
- **Modified**: `bin/vcf_stats/seq2neo/visualizer.py` — 10+ new altair charts, dashboard sections
- **Modified**: `bin/vcf_stats/seq2neo/cli.py` — wire Rust core, CLI flags, 4-level outputs
- **Modified**: `bin/vcf_stats/seq2neo/bam_stats.py` — 3 BAM types (DN/DT/RT)
- **Modified**: `pyproject.toml` — maturin build dependency, pysam dev dependency
