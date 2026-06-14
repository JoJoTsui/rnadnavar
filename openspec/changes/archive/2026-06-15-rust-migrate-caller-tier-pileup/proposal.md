## Why

The seq2neo variant statistics pipeline has 4 remaining Python/GIL-bound bottlenecks after fixing rescue VCF and BAM stats: caller VCF parsing (17s/sample), tiering computation (13s/sample), BAM pileup validation (no Rust path), and a slow Python for-loop join (3s/sample). Together these account for ~60% of per-sample processing time. Moving them to Rust with proper GIL release, using normalized caller VCFs already produced by the pipeline, will cut per-sample time from ~50s to ~5s.

## What Changes

- **Caller VCF parsing moved to Rust**: Replace cyvcf2-based `parse_single_caller()` with a Rust noodles-vcf implementation that reads bcftools-normalized caller VCFs (`*.dec.norm.vcf.gz`). Match variants on (CHROM, POS, REF, ALT) against the rescue VCF — not just (CHROM, POS). Extract all FORMAT fields per caller type. Release GIL during parsing.
- **join_caller_columns optimized with polars**: Replace Python for-loops building pl.Series with a polars `.join()` on all 4 coordinate columns (CHROM, POS, REF, ALT). This fixes a latent bug where multiallelic sites could match the wrong allele.
- **Tiering computation moved to Rust**: Implement the CxDy tiering logic (concordant caller counting, database evidence checking, tier assignment from `tier_config.py`) as a Rust `#[pyfunction]` with full parity to the Python `TieringEngine`. Release GIL during computation.
- **BAM pileup moved to Rust**: Implement per-position BAM pileup using noodles-bam indexed queries (`reader.query()` with BAI/CSI index). Extract DP, REF_DP, ALT_DP, F1R2/F2R1 strand counts, mean BQ, mean MQ at each variant position. Release GIL.
- **CALLER_CONFIGS updated**: Point to normalized caller VCF paths instead of raw filtered VCFs.

## Capabilities

### New Capabilities
- `caller-vcf-rust-parser`: Rust-based parsing of 6 caller VCFs (Mutect2, DeepSomatic, Strelka × DNA/RNA) from bcftools-normalized files, extracting FORMAT fields matched on (CHROM, POS, REF, ALT) coordinates
- `bam-pileup-rust`: Rust-based per-position BAM pileup using noodles-bam indexed queries, extracting depth, allele counts, strand bias, and quality metrics

### Modified Capabilities
- `variant-tiering-stats`: Tiering computation backend changed from Python TieringEngine to Rust implementation with identical output — requirement unchanged, implementation replaced
- `multi-level-aggregation`: join_caller_columns now joins on (CHROM, POS, REF, ALT) instead of (CHROM, POS) — fixes multiallelic site matching; caller VCF paths updated to normalized files

## Impact

- **Rust crate**: 3 new source files (`caller.rs`, `tier.rs`, `pileup.rs`) + modified `lib.rs`
- **Python code**: `caller_parser.py` (wire Rust + polars join), `tiering_stats.py` (wire Rust), `rust_bam.py` (wire Rust pileup), `manifest_loader.py` (normalized VCF paths)
- **Tests**: `test_seq2neo_stats.py` — new classes for Rust caller parsing, Rust tiering, Rust pileup; extend TestGILRelease; update TestCallerParser for 4-column join
- **Dependencies**: Existing noodles-vcf 0.88, noodles-bam 0.90, noodles-csi 0.42 — no new deps needed
- **Data**: Uses normalized VCFs at `normalized/<caller>/*.dec.norm.vcf.gz` and `vcf_realignment/normalized/<caller>/*.dec.norm.vcf.gz` — already produced by pipeline, read-only
