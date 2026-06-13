## Why

BAM pileup (variant-wise depth/strand/quality validation from alignment files) is enabled by default but has three critical issues that make it unusable in production:

1. **Performance**: Per-position BAI query loops over 1.4M positions × 3 BAM types = 4.2M individual BAM seeks. A single large sample takes ~9 hours to complete, and the process appears hung (no progress output during GIL-released Rust execution).

2. **Data correctness**: Pileup results are computed per (CHROM, POS, REF, ALT) but joined back to the variant DataFrame on only (CHROM, POS). At 3,859 multiallelic sites (up to 16 alleles each), pileup data is cross-joined to wrong alleles — 129,430 rows affected in a single sample.

3. **Dead functionality**: `--pileup-mode filtered` flag (exclude NoConsensus variants) is defined but never wired, making it impossible to focus pileup on high-confidence calls.

## What Changes

### Fix 1 — Windowed BAM queries in Rust (performance)
- Group positions into ~1Mb genomic windows per chromosome
- One BAI query per window instead of per position (~3K queries vs 1.4M)
- Same pileup logic within each window, matching reads to positions via HashMap
- Estimated: ~30 seconds per BAM type instead of ~3 hours

### Fix 2 — 4-column join key (correctness)
- Add REF and ALT columns to pileup output (`rust_bam.py`)
- Join pileup results on `["CHROM", "POS", "REF", "ALT"]` instead of `["CHROM", "POS"]` (`cli.py`)
- Eliminates cross-join at multiallelic sites

### Fix 3 — Wire `--pileup-mode` flag
- When `--pileup-mode filtered`, filter positions to exclude `FILTER == "NoConsensus"` before pileup
- Adds progress logging before and after each BAM type's pileup

### Fix 4 — Tests
- 7 new test methods covering: 4-column output, multiallelic join, windowed query parity, filtered mode, large-batch performance, integration, join key verification

## Impact

- **Rust**: `pileup.rs` — rewrite per-position loop to windowed queries (~80 lines changed)
- **Python**: `rust_bam.py` — add REF/ALT to output; `cli.py` — 4-column join, wire pileup_mode, add logging
- **Tests**: `test_seq2neo_stats.py` — 7 new test methods
- **Performance**: 9 hours → ~2 minutes per sample for pileup
