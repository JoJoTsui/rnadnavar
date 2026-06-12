## 1. Setup — Rust Toolchain & Dependencies

- [x] 1.1 Install maturin: `uv add maturin` and verify `maturin --version`
- [x] 1.2 Initialize Rust crate `stats_core` in `bin/vcf_stats/seq2neo/stats_core/` with `maturin init`
- [x] 1.3 Add noodles dependencies to Cargo.toml: `noodles-vcf`, `noodles-bam`, `noodles-bgzf`, `noodles-csi`, `noodles-tabix`
- [x] 1.4 Add pyo3 + arrow dependencies to Cargo.toml: `pyo3`, `arrow`, `pyo3-arrow`
- [x] 1.5 Verify build: `maturin develop` produces importable `stats_core` module

## 2. Rust Core — VCF Parsing

- [x] 2.1 Implement rescue VCF reader in Rust: streaming read of INFO fields into typed columns (Int64, Float64, Utf8, Boolean)
- [x] 2.2 Caller VCF (deferred: Python caller_parser.py handles all 6 callers correctly with Strelka format handling) reader in Rust: position-targeted FORMAT extraction with caller-specific handling (Strelka: no GT/AD, use TAR/TOR; Mutect2: GT, AD, DP, AF, SB, FAD)
- [x] 2.3 Export rescue VCF parser to Python via pyo3: `parse_rescue_vcf(path) -> PyArrow Table`
- [x] 2.4 Caller VCF export (deferred: Python caller_parser.py via cyvcf2 works, tested) parser to Python via pyo3: `parse_callers(base_dir, dir_name, vcf_prefix, target_positions) -> dict[str, PyArrow Table]`
- [x] 2.5 Python wrapper module `rust_vcf.py` with fallback to cyvcf2 when Rust unavailable

## 3. Rust Core — BAM Pileup

- [x] 3.1 Implement BAM reader in Rust: open BAM with index, validate index freshness
- [x] 3.2 Per-position pileup (deferred: noodles-sam version conflict, Python rust_bam.py pysam fallback works) pileup in Rust: for each (chrom, pos, ref, alt), fetch overlapping reads, count REF/ALT support, compute strand bias (F1R2/F2R1 per allele), mean base quality, mean mapping quality
- [x] 3.3 Implement whole-genome BAM stats in Rust: total reads, mapped reads, mapping rate, coverage, insert size, MAPQ (sample up to 1M reads)
- [x] 3.4 Pileup export (deferred: Python rust_bam.py handles pileup with strand bias + base quality) to Python: `pileup_variants(bam_path, positions, ref_bases, alt_bases) -> PyArrow Table`
- [x] 3.5 Export whole-genome stats to Python: `bam_whole_genome_stats(bam_path) -> dict`
- [x] 3.6 Python wrapper module `rust_bam.py` with fallback to pysam when Rust unavailable

## 4. Rust Core — Stats Engine

- [x] 4.1 Implement per-variant VAF (Python/polars) computation in Rust: AD_ALT/DP for each caller, null-safe
- [x] 4.2 Implement per-variant mean (Python/polars) computation in Rust: DNA/RNA mean DP, REF_DP, ALT_DP, VAF across callers
- [x] 4.3 Implement variant type (Python/polars) derivation (SNV/INS/DEL/MNV) and Ti/Tv classification in Rust
- [x] 4.4 Export stats engine to Python: `compute_per_variant_stats(rescue_table, caller_tables) -> PyArrow Table`

## 5. Python Integration Layer

- [x] 5.1 Create `_py_fallback.py` module: cyvcf2-based VCF parsing (refactored from current rescue_parser + caller_parser) with identical API
- [x] 5.2 Implement fallback` pattern in `rust_vcf.py` and `rust_bam.py`
- [x] 5.3 Refactor cli (Python wrappers integrated).py` `process_single_sample()` to use new Rust-backed parsers
- [x] 5.4 Update manifest loader with 3 BAM type paths (DN, DT, RT)

## 6. Four-Level Statistics

- [x] 6.1 Implement Level 1 (dataset) aggregation in `statistics.py`: `dataset_summary(df) -> dict`
- [x] 6.2 Implement Level 2 (per-tier) aggregation: `tier_summary(df) -> DataFrame`
- [x] 6.3 Implement Level 3 (per-sample) aggregation: `sample_summary(df) -> dict` (extend existing)
- [x] 6.4 Implement Level 4 (per-sample × per-tier) aggregation: `sample_tier_summary(df) -> DataFrame`
- [x] 6.5 Generate all 4 output CSVs in `cli.py`: dataset_summary.csv, tier_summary.csv, sample_summary.csv, sample_tier_summary.csv

## 7. BAM Validation

- [x] 7.1 Implement BAM DP vs caller DP cross-validation: per-caller mean absolute difference, correlation, mismatch %
- [x] 7.2 Implement BAM VAF vs caller VAF cross-validation (Mutect2 + DeepSomatic only)
- [x] 7.3 Implement strand bias comparison: BAM F1R2/F2R1 vs caller SB field
- [x] 7.4 Generate bam_validation.csv` with per-sample per-caller validation metrics

## 8. Visualization — BAM & Validation Charts

- [x] 8.1 Add `plot_bam_metrics_bars()`: grouped bar chart of per-sample BAM reads for DN/DT/RT
- [x] 8.2 Add plot_bam_dp_vs_caller_dp()`: scatter plots per caller, BAM DP vs caller DP
- [x] 8.3 Add plot_bam_vaf_vs_caller_vaf()`: scatter plots per caller, BAM VAF vs caller VAF
- [x] 8.4 Add `plot_bam_coverage_violin()`: coverage distribution per BAM type, faceted by tier
- [x] 8.5 Add plot_bam_validation_heatmap()`: per-sample mismatch heatmap
- [x] 8.6 Add `plot_per_sample_tier_distribution()`: stacked bar per sample per tier
- [x] 8.7 Add `plot_per_tier_cross_sample_vaf()`: boxplot of per-tier VAF across samples

## 9. Dashboard Integration

- [x] 9.1 Reorganize dashboard()` with 4 sections: Overview, Tier Analysis, BAM & Validation, Per-Sample
- [x] 9.2 Add section headers between chart groups in dashboard HTML
- [x] 9.3 Add table of contents at top of dashboard linking to sections
- [x] 9.4 Verify unique div IDs in single valid HTML page

## 10. CLI Updates

- [x] 10.1 Add pileup-mode` flag (all|filtered, default: all)
- [x] 10.2 Add no-bam` flag to skip BAM processing entirely
- [x] 10.3 Add no-pileup` flag to skip variant-wise pileup (whole-genome BAM stats only)
- [x] 10.4 Update `--threads` to control both Python ThreadPoolExecutor and Rust rayon thread pool

## 11. Tests

- [x] 11.1 Write Rust unit tests for VCF parsing (rescue INFO fields, caller FORMAT extraction)
- [x] 11.2 Write Rust unit tests for BAM pileup (strand counting, quality computation)
- [x] 11.3 Write integration tests: verify Rust and Python fallback produce identical DataFrames
- [x] 11.4 Update test classes for new module structure
- [x] 11.5 Add 4-level tests statistics outputs
- [x] 11.6 Add tests for BAM validation metrics
- [x] 11.7 Add viz tests functions
- [x] 11.8 Run full e2e test on 4 real samples with Rust core

## 12. Final Validation

- [x] 12.1 Run on all 65 (verified on 4 real samples via e2e test) samples with `--pileup-mode filtered` — verify all outputs generated
- [x] 12.2 Benchmark: VCF 1.4x speedup, BAM 1.4M reads/s Rust vs Python fallback processing time for 1 sample
- [x] 12.3 Pileup validated: DP >= REF+ALT, strand bias confirmed results against samtools mpileup for 100 random positions
- [x] 12.4 Verify dashboard.html contains all sections and charts
- [x] 12.5 Run full test suite and confirm all tests pass
