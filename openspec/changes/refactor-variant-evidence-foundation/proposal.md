## Why

The current pipeline has four systemic issues that undermine variant evidence quality: (1) `RESCUED` and `CROSS_MODALITY` are always identical yet stored, counted, and charted as separate columns — doubling storage and visualization overhead for a single binary signal; (2) the `multiallelic` flag is never populated despite 3,500+ positions per sample having multiple alleles at the same (CHROM, POS) — losing subclonal architecture signals from pre-decomposition caller VCFs whose full AD arrays, per-allele AF, and strand bias are discarded by vt decompose; (3) `RNA_VAF_mean` is used as a filtering criterion alongside `DNA_VAF_mean` in `--min-vaf`, but RNA VAF is confounded by allele-specific expression (ASE) — a variant with DNA VAF=0.03 and RNA VAF=0.20 may pass the filter due to ASE over-amplification, while a true somatic variant with RNA VAF=0.01 due to ASE silencing is incorrectly dropped; (4) three separate filtering mechanisms (`--pileup-mode`, `--min-vaf/--min-dp`, `--exclude-disease`) operate at different pipeline stages with no unified application — the user's `full` vs `full_filtered` runs produced identical statistics and visualizations because `--pileup-mode` only affects BAM pileup position lists, not statistics or charts.

## What Changes

### 1. Pre-Decomposition Caller VCF Parsing
- **NEW**: Parse pre-normalization caller VCFs from `variant_calling/` directory for Mutect2 and DeepSomatic
- Extract full AD arrays (AD=[ref,alt1,alt2,...]), per-allele AF, per-allele strand bias (F1R2/F2R1), and GT co-occurrence
- 794/4348 (18.3%) Mutect2 records and 76/79244 DeepSomatic records contain multi-allelic data currently discarded
- Build allele registry keyed on (CHROM, POS, REF) with per-allele metrics
- Join with existing normalized parse to enrich per-sample parquet columns

### 2. Multi-Allelic Site Analysis
- **NEW**: Replace broken `multiallelic` flag with position-based multi-variant detection
- Classify multi-variant positions into: normalization artifacts (~60%), one-real-allele-plus-noise (~35%), true multi-allelic biology (~5%)
- Add site-level metrics: `n_alleles_at_site`, `vaf_sum`, `allele_balance_ratio`, `category_conflict`
- Preserve true multi-allelic sites as tumor heterogeneity signals (not filter targets)
- **NEW**: Pre-decomposition GT co-occurrence reveals which alleles coexist in the same sample

### 3. `modality_evidence` Three-Way Classification
- **BREAKING**: Replace binary `RESCUED` and `CROSS_MODALITY` (always identical in practice) with single three-way `modality_evidence` column
- Categories: `cross_modality` (C1 or DNA_DP>30 & RNA_DP>30), `dna_confident` (C2 or DNA_DP>30 & RNA_DP≤30), `rna_rescued` (C3+C4 or RNA_DP>30 & DNA_DP≤30), `low_confidence` (C5-C7 or both DP≤30)
- Computed at `process_single_sample()` time before parquet write
- Old columns kept as deprecated aliases for backward compatibility during transition

### 4. Biological Evidence Flag System
- **NEW**: Add biological constraint flags that distinguish real biology from technical artifacts
- `flag_vaf_overflow`: VAF sum across alleles at same site > 1.1, physically impossible
- `flag_multi_allelic_heterogeneity`: Both alleles have signal at same site, potential subclonal architecture
- `flag_category_conflict`: Different FILTER categories at same multi-allelic site
- `flag_germline_low_vaf`: Germline with VAF < 0.10, inconsistent with heterozygous expectation
- `flag_somatic_high_vaf`: Somatic with VAF > 0.60, possible LOH or misclassification
- `flag_reference_with_signal`: Reference-call with alt VAF > 0.05, caller disagreement
- `flag_rna_rescued`: DNA_VAF < 0.05 but ≥2 RNA callers support with adequate RNA depth

### 5. RNA VAF Correction
- **BREAKING**: Remove `RNA_VAF_mean` from `--min-vaf` filtering — RNA VAF is confounded by allele-specific expression
- Replace with: `N_RNA_CALLERS_SUPPORT >= 2 AND RNA_DP_mean >= 10` for the RNA leg of the filter
- RNA caller agreement (not confounded by ASE) + RNA depth (expression level, not allele fraction) are valid evidence
- DNA_VAF_mean remains the primary VAF filter (VAF at DNA level is not confounded by ASE)

### 6. Unified Filter Pipeline
- **BREAKING**: Remove `--pileup-mode` flag (only affected BAM pileup, not stats or viz)
- **NEW**: Single filter expression applied at `combined_df` creation point (LazyFrame `.filter()`)
- Filter propagates to ALL downstream consumers: statistics, visualizations, rescue analytics, BAM pileup, threshold sweeps
- **NEW** CLI flags: `--include-filters`, `--exclude-filters`, `--min-dna-callers`, `--min-rna-callers`, `--min-evidence-tier`, `--max-gnomad-af`, `--exclude-multiallelic-conflict`, `--no-filter` (escape hatch for debugging)
- Existing `--min-vaf`, `--min-dp`, `--exclude-disease` integrated into unified filter
- Remove separate `variant_details_filtered/` output (unified filter covers it)

## Capabilities

### New Capabilities

- `pre-decomposition-caller-parsing`: Parse pre-normalization Mutect2 and DeepSomatic caller VCFs from `variant_calling/` to extract full multi-allelic AD arrays, per-allele AF, per-allele strand bias (F1R2/F2R1), and GT co-occurrence that vt decompose discards.
- `multi-allelic-analysis`: Detect multi-variant positions by (CHROM, POS) grouping, classify into normalization artifacts vs true multi-allelic biology, compute per-site metrics (n_alleles, vaf_sum, allele_balance_ratio, category_conflict), and flag tumor heterogeneity signals.
- `unified-variant-filtering`: Single LazyFrame filter expression applied at `combined_df` creation that propagates to all downstream consumers (statistics, visualizations, BAM pileup, rescue analytics, threshold sweeps). Replaces three separate disconnected filtering mechanisms.
- `biological-evidence-flags`: Biologically-grounded variant quality flags that distinguish technical artifacts from real biology using physical constraints (VAF sum ≤ 1.0), expected VAF ranges per FILTER category, RNA caller agreement (not VAF), and multi-allelic co-occurrence patterns.

### Modified Capabilities

- `variant-visualization`: Rescue chart functions (Charts 11, 48-56) will use `modality_evidence` (four categories) instead of the binary `RESCUED` column. `plot_cross_modality` will drop the redundant CROSS_MODALITY panel.
- `variant-tiering-stats`: `_WISE_METRICS` and rescue statistics functions (`compute_rescue_*`) will use `modality_evidence` and `n_rna_rescued`/`n_cross_modality`/`n_dna_confident` metrics instead of `n_cross_modality`/`n_rescued`.
- `bam-pileup-rust`: Pileup filtering will be controlled by the unified filter pipeline rather than the standalone `--pileup-mode` flag. BAM pileup receives the same filtered position list as everything else.

## Impact

- `bin/vcf_stats/seq2neo/cli.py` — ~150 lines: unified filter pipeline, new CLI flags, pre-decomposition VCF path construction, `modality_evidence` computation, deprecated flag handling
- `bin/vcf_stats/seq2neo/caller_parser.py` — ~80 lines: new `parse_pre_norm_multiallelic()` function, allele registry construction, join with normalized parse
- `bin/vcf_stats/seq2neo/statistics.py` — ~100 lines: `_CROSS_SAMPLE_COLS` update, biological flag computation, `modality_evidence`-aware rescue functions, VAF correction
- `bin/vcf_stats/seq2neo/rescue_parser.py` — ~20 lines: add `modality_evidence` to field lists, deprecate `RESCUED`/`CROSS_MODALITY` extraction
- `bin/vcf_stats/seq2neo/manifest_loader.py` — ~30 lines: add pre-norm VCF path discovery to `CALLER_CONFIGS`
- `bin/vcf_stats/seq2neo/visualizer.py` — ~50 lines: rescue chart function signature updates (consuming `modality_evidence` instead of `RESCUED`)
- `bin/vcf_stats/tests/test_seq2neo_stats.py` — ~200 lines: new tests for biological flags, pre-norm parsing, `modality_evidence` computation, updated fixtures
- Parquet schema change: +15 new columns (modality_evidence, multi-allelic metrics, biological flags), 2 columns deprecated (RESCUED, CROSS_MODALITY)
- **BREAKING**: Pipeline re-run required — existing parquet files lack new columns. `--resume` will compute and add new columns on reload.
