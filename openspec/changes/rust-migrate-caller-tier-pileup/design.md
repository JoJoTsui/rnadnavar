## Context

The seq2neo statistics pipeline parses 6 caller VCFs per sample (Mutect2, DeepSomatic, Strelka × DNA/RNA), joins their FORMAT fields with the rescue VCF, computes CxDy tiers, and performs BAM pileup validation. Currently all 4 steps run in Python with cyvcf2/pysam, holding the GIL. After fixing `parse_rescue` and `bam_stats` GIL issues, these remain the bottlenecks: caller parsing (17s), joining (3s), tiering (13s), pileup (many minutes for 200K+ positions). The pipeline already produces bcftools-normalized caller VCFs (`*.dec.norm.vcf.gz`) with CSI indexes that are currently unused.

## Goals / Non-Goals

**Goals:**
- Move caller VCF parsing to Rust (noodles-vcf) with GIL release, reading normalized VCFs
- Replace join_caller_columns Python loop with polars `.join()` on (CHROM, POS, REF, ALT)
- Move tiering computation to Rust with identical output to Python TieringEngine
- Move BAM pileup to Rust (noodles-bam indexed queries) with GIL release
- Update CALLER_CONFIGS to use normalized VCF paths
- Tests that verify Rust output matches Python output exactly

**Non-Goals:**
- Changing the tiering algorithm — must produce identical results to Python
- Adding new FORMAT fields beyond what caller_parser.py already extracts
- Supporting non-normalized caller VCFs — normalized paths are the new default
- Pileup for whole-genome positions (only variant positions from rescue VCF)
- Replacing cyvcf2/pysam imports completely — fallback paths preserved

## Decisions

### 1. Use normalized caller VCFs (*.dec.norm.vcf.gz)

**Decision:** Switch CALLER_CONFIGS from raw filtered VCFs to bcftools-normalized VCFs.

**Rationale:** Normalized VCFs have consistent REF/ALT representation matching the rescue VCF. This enables exact (CHROM, POS, REF, ALT) matching instead of just (CHROM, POS), fixing a latent bug at multiallelic sites. The files are already produced by the pipeline — no additional processing needed.

**Paths:**
```
DNA: normalized/{caller}/{prefix}DT_vs_{prefix}DN/*.dec.norm.vcf.gz
RNA: vcf_realignment/normalized/{caller}/{prefix}RT_realign_vs_{prefix}DN/*.dec.norm.vcf.gz
```

### 2. 4-column matching for caller lookup

**Decision:** Match caller variants to rescue variants on (CHROM, POS, REF, ALT), not just (CHROM, POS). The polars join uses all 4 columns.

**Rationale:** At multiallelic sites, a caller VCF may have multiple records at the same position with different REF/ALT. Using only (CHROM, POS) can match the wrong allele. Normalized VCFs guarantee consistent REF/ALT representation, making exact matching safe.

**Alternative considered:** Use only (CHROM, POS) and rely on caller VCFs being decomposed. Rejected because it's fragile in edge cases.

### 3. Tiering in Rust with hardcoded rules

**Decision:** Implement tier rules from `tier_config.py` directly in Rust rather than calling Python from Rust.

**Rationale:** The tier rules are static data (7 caller tiers × 2 database tiers = 14 final tiers) with simple boolean conditions. Hardcoding in Rust avoids pyo3 overhead and GIL roundtrips. The quality scores, tier colors, and display names stay in Python (not performance-critical).

**Parity requirement:** Must exactly reproduce `count_concordant_callers()` and `compute_database_tier()` logic. FILTERS_NORMALIZED string parsing must match the Python regex. Database thresholds must match: gnomAD_AF > 0.001, COSMIC_CNT > 0, REDI_EVIDENCE truthy.

### 4. BAM pileup via noodles-bam indexed queries

**Decision:** Use `bam::io::Reader::query()` with CSI index for per-position BAM pileup.

**Rationale:** noodles-bam 0.90 supports indexed queries via `reader.query(header, index, region)`. The CSI index exists alongside normalized VCFs. Per-position pileup is inherently position-by-position — indexed queries avoid full BAM scans.

**Alternative considered:** Sequential BAM scan with position filter. Rejected because full scan of 200M+ reads per BAM is wasteful when we need data at specific positions.

**Alternative considered:** rust-htslib. Rejected because noodles-bam is already in the dependency tree and avoids C library linking complexity.

### 5. Column vectors for Rust↔Python data transfer

**Decision:** Use `.to_list()` (Python → Rust) and return dicts from Rust (Rust → Python). Avoid pyo3-arrow for now.

**Rationale:** For 200K variants, column vectors are ~few MB — acceptable memory overhead. pyo3-arrow would avoid copies but adds complexity and potential version conflicts. The simpler approach is sufficient for this scale.

### 6. Rust module structure

**Decision:** 3 new source files in stats_core: `caller.rs`, `tier.rs`, `pileup.rs`. Each exports a single `#[pyfunction]` with `py.detach()`.

**Rationale:** Follows the existing pattern from `vcf.rs` and `bam.rs`. One file per concern. GIL release via `py.detach()` matches the pattern recently fixed in `lib.rs`.

## Risks / Trade-offs

- **[Risk] Tiering output mismatch:** Rust implementation produces different tier assignments than Python for edge cases → **Mitigation:** Run both implementations on 1000 random variants and assert exact match before switching.
- **[Risk] Normalized VCF missing for some samples:** Not all samples may have normalized VCFs → **Mitigation:** Fall back to original caller VCF paths if normalized file not found. Log warning.
- **[Risk] noodles-sam version conflict with pileup:** Previous attempt was blocked → **Mitigation:** Verify `reader.query()` compiles and runs on a real BAM before implementing full pileup. If blocked, defer pileup to separate change.
- **[Risk] 4-column join increases complexity for callers without REF/ALT:** Strelka reports TAR/TIR/TOR but may represent alleles differently → **Mitigation:** Normalized VCF guarantees consistent representation. If matching fails, the variant gets null-filled caller columns (graceful degradation).

## Migration Plan

1. Add Rust functions to stats_core, build, install
2. Add Python wrappers with Rust-primary + Python-fallback pattern
3. Update CALLER_CONFIGS with normalized paths
4. Update join_caller_columns to polars join on 4 columns
5. Run test suite, verify output parity
6. Run on 4-sample test set, compare with previous output

## Open Questions

- None — all decisions resolved during exploration phase
