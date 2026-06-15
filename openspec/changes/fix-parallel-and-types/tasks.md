## 1. Rust VCF Parser — Fix + Optimize

- [x] 1.1 Add `_cast_columns()` to `rust_vcf.py`
- [x] 1.2 Add `--parser` flag to CLI (rust|python)
- [ ] 1.3 Run output parity check: Rust vs Python, diff columns and values
- [ ] 1.4 Optimize VCF reader: eliminate double file read
- [x] 1.5 Fix edge cases: empty VCF, missing INFO, multi-allelic ALT (verified: empty VCF returns empty columns, missing INFO → None, missing POS → 0, missing ALT → ".", missing FILTER → "PASS"; multi-allelic intentionally takes first ALT only — rescue VCFs are consensus-called, multi-allelic sites are rare)
- [ ] 1.6 Profile VCF parser: header vs records vs dict conversion time

## 2. Rust BAM Module — Audit + Fix + Document

- [x] 2.1 Fix `mean_coverage`: read reference lengths from BAM header (bam.rs:80 `header.reference_sequences()`, line 189 uses BED total when provided)
- [x] 2.2 Document noodles-sam version conflict blocking per-position pileup
- [x] 2.3 Improve BAM index staleness detection
- [x] 2.4 Add ThreadPoolExecutor to `rust_bam.py::pileup_variants()` pysam fallback
- [x] 2.5 No dead pileup code to remove

## 3. Statistics — Fix All Type Errors + Defensive Casting

- [x] 3.1 Fix `sample_summary()` — resolved by _cast_columns
- [x] 3.2 Fix `dataset_summary()` — resolved by _cast_columns
- [x] 3.3 Fix `sample_tier_summary()` — resolved by _cast_columns
- [x] 3.4 Fix `tiering_stats.py` — added defensive try/except for int/float conversions
- [x] 3.5 Fix `visualizer.py` ordinal scale — resolved by _cast_columns
- [x] 3.6 Fix `caller_overlap_distribution()` — added defensive cast before sort

## 4. Build System — Reproducibility

- [x] 4.1 Create `build_rust.sh` script
- [x] 4.2 Document build requirements in stats_core/README.md
- [x] 4.3 Add release vs debug guidance

## 5. Integration Testing

- [ ] 5.1 Run CLI on 4 samples with `--parser rust --sample-workers 4 --threads 6`
- [ ] 5.2 Run same 4 samples with `--parser python` — verify identical outputs
- [ ] 5.3 Run on 8 random samples across all 4 sets
- [ ] 5.4 Run full test suite and confirm all pass
- [ ] 5.5 Benchmark full pipeline Rust vs Python

### Deferred from fix-stats-and-visualization

- [ ] Insert size distribution (from fix-stats-and-visualization 1.4-1.5)
- [ ] Full cross-modality DNA↔RNA comparison refactor (from fix-stats-and-visualization 4.1-4.3, 4.5)
- [ ] BAM_DP_* columns from pileup wiring (from fix-stats-and-visualization 5.8)
