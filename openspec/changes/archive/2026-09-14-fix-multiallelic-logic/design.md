## Context

The archived `2026-06-22-refactor-variant-evidence-foundation` proposal specified a multi-allelic analysis system with pre-decomposition VCF parsing, per-allele classification, site-level metrics, biological flags, and modality_evidence classification. The implementation skeleton exists (functions, columns, CLI flags, plot registration) but the core logic diverges from the spec in 8 critical ways:

1. **Classification**: Uses `allele_balance_ratio < 0.03` / `vaf_sum < 0.15` / catch-all instead of per-allele REF/ALT-length, VAF/DP, and ALT-base tests. The length-difference criterion for normalization_artifact is entirely missing. Result: 0.03% normalization_artifact, 73% noise, 27% true_multi_allelic — inverted from the expected 60%/35%/5%.
2. **allele_balance_ratio**: Uses `alt_dp_2nd_max / alt_dp_max` instead of `max_alt_dp / total_alt_dp`. For ALT_DP=[15,5], code gives 0.333 instead of 0.75.
3. **flag_multi_allelic_heterogeneity**: Includes `normalization_artifact` in the flag; spec says only `true_multi_allelic`.
4. **flag_category_conflict**: Plain copy of `category_conflict` without the "both categories biological" guard; no `category_conflict_resolution` column.
5. **--exclude-multiallelic-conflict**: Checks `flag_vaf_overflow` instead of `flag_category_conflict`.
6. **Pre-norm data dead code**: `parse_pre_norm_multiallelic` extracts AD arrays, AF, F1R2/F2R1, GT — but `pre_norm_*` columns are not in `_CROSS_SAMPLE_COLS` and no downstream function references them.
7. **modality_evidence_caller**: Uses `>=1 DNA AND >=1 RNA` for `cross_modality` instead of C1-C7 tier-based mapping (C1=cross_modality, C2=dna_confident, C3+C4=rna_rescued, C5-C7=low_confidence).
8. **Visualizations**: `allele_balance` plot shows per-caller VAF max/min (not per-allele); `category_conflict` plot uses cross-sample group_by on a per-sample flag producing nonsensical single-category "conflicts".

Additionally, cyvcf2 is used as a fallback in `caller_parser.py` and `rescue_parser.py` but the Rust parser is the only backend in use. The cyvcf2 fallback adds maintenance burden and was found to have a silent join-failure bug (empty REF/ALT). Per the user's direction, cyvcf2 should be removed from the seq2neo stats subsystem (it remains in the main workflow's `bin/run_consensus_vcf.py`, `bin/filter_vcf.py`, etc.).

## Goals / Non-Goals

**Goals:**
- Fix multi-allelic classification to match the spec's per-allele criteria
- Fix allele_balance_ratio formula to `max_alt_dp / total_alt_dp`
- Fix 3 biological flags to match spec (heterogeneity, category_conflict, add strand_bias)
- Wire pre-decomposition data into downstream metrics (vaf_sum, allele_balance, gt_cooccurrence, strand_balance)
- Unify modality_evidence_caller with C1-C7 tier-based mapping
- Fix 2 multi-allelic visualizations
- Remove cyvcf2 fallback from seq2neo stats

**Non-Goals:**
- Fixing the Rust tier fallback for RNAedit (that belongs in `fix-bam-stats-rust`)
- Fixing BAM stats metrics (MAPQ, coverage_bins, expanded metrics) — that belongs in `fix-bam-stats-rust`
- Adding new filtering logic (that belongs in `add-low-confidence-variant-filtering`)
- Changing the upstream Nextflow pipeline or VCF normalization
- Removing cyvcf2 from the main workflow (`bin/run_consensus_vcf.py`, `bin/filter_vcf.py`, etc.) — only the seq2neo stats subsystem

## Decisions

### Decision 1: Classification uses per-allele criteria, not aggregate thresholds

**Chosen**: Implement the spec's three-tier classification:
- `normalization_artifact`: Different REF/ALT string lengths at the same position OR one allele has VAF≈0 (≤0.001) AND DP≈0 (≤1)
- `noise`: Exactly one allele has VAF≥0.01 AND DP≥5; all others have VAF<0.01 OR DP<5
- `true_multi_allelic`: ≥2 alleles with VAF≥0.01 AND DP≥5 AND different ALT base sequences

The `when()`-chain ordering: check normalization_artifact first (length + zero-signal), then noise (single-signal), then true_multi_allelic (multi-signal).

**Alternatives considered**:
- Keep aggregate thresholds but tune them → Rejected: the fundamental problem is aggregate vs per-allele; no threshold tuning fixes the missing length criterion
- Use a machine learning classifier → Rejected: over-engineered for a 3-class biological rule

### Decision 2: allele_balance_ratio = max_alt_dp / total_alt_dp

**Chosen**: Change `statistics.py:444` from `alt_dp_2nd_max / alt_dp_max` to `alt_dp_max / total_alt_dp`. This naturally handles 3+ alleles (the denominator includes all alleles, not just top-2). Use pre-norm `ad_alts` data when available; fall back to `DNA_ALT_DP_mean` from decomposed records.

**Alternatives considered**:
- Use VAF-based ratio (`max_vaf / sum_vaf`) → Rejected: spec's primary definition is ALT_DP-based; VAF is a secondary definition in one section of the spec

### Decision 3: modality_evidence_caller uses C1-C7 tier mapping

**Chosen**: Map directly from `final_tier`:
- C1D0, C1D1 → `cross_modality`
- C2D0, C2D1 → `dna_confident`
- C3D0, C3D1, C4D0, C4D1 → `rna_rescued`
- C5D0, C5D1, C6D0, C6D1, C7D0, C7D1 → `low_confidence`

This is a simple lookup, not a condition chain. The `modality_evidence_dp` column remains depth-based (DNA_DP > 30 AND RNA_DP > 30).

**Alternatives considered**:
- Keep raw caller counts but change thresholds → Rejected: tiers already encode the concordant caller counts; re-deriving from raw counts duplicates logic and can diverge from tier assignment (especially with the Rust fallback bug)

### Decision 4: Pre-norm data added to _CROSS_SAMPLE_COLS selectively

**Chosen**: Add `pre_norm_DNA_mutect2_n_alleles`, `pre_norm_DNA_mutect2_af`, `pre_norm_DNA_deepsomatic_af`, and the combined `gt_cooccurrence`, `strand_balance` to `_CROSS_SAMPLE_COLS`. Do NOT add raw `ad_ref`/`ad_alts` arrays (they are list-typed and not suitable for cross-sample aggregation). Instead, compute derived scalar metrics (`gt_cooccurrence`, `strand_balance`) at per-sample time and add those to the cross-sample schema.

**Alternatives considered**:
- Add all pre_norm_* columns → Rejected: list-typed columns break cross-sample scans and aggregation
- Don't use pre-norm at all → Rejected: the proposal's central motivation was recovering multi-allelic data lost to vt decompose

### Decision 5: Remove cyvcf2 from seq2neo stats only

**Chosen**: Remove `from cyvcf2 import VCF` and all fallback code paths in `caller_parser.py:23,296-322,488` and `rescue_parser.py:3,12,146+`. Remove the `--vcf-parser` CLI flag at `cli.py:802`. The Rust parser (`rust_vcf.py`, `stats_core.so`) becomes the only backend. If the Rust parser is unavailable, the pipeline errors with a clear message instead of silently falling back to broken cyvcf2 code.

**Alternatives considered**:
- Keep cyvcf2 as a documented fallback → Rejected: the fallback has a silent join-failure bug (empty REF/ALT) and the user confirmed Rust is the only backend
- Fix the cyvcf2 fallback → Rejected: maintaining two backends doubles the testing surface; Rust is faster and already handles all cases

## Risks / Trade-offs

- **[Parquet schema change]**: Adding `gt_cooccurrence`, `strand_balance`, `flag_extreme_strand_bias`, `category_conflict_resolution` columns changes the parquet schema. `--resume` will need to re-process samples lacking these columns. → **Mitigation**: `--resume` checks for `gt_cooccurrence` column presence; if absent, forces re-processing of that sample.
- **[modality_evidence BREAKING]**: Variants previously classified as `cross_modality` (1 DNA + 1 RNA) will become `low_confidence` or `rna_rescued`. Downstream stats and charts will show different distributions. → **Mitigation**: Document the change in SUMMARY.md; the old classification was wrong (1+1 is not "cross-modality strong").
- **[Pre-norm VCF availability]**: Pre-decomposition VCFs may not exist for all samples. → **Mitigation**: Pre-norm parsing remains best-effort; null-fill derived columns when pre-norm data is absent.
- **[Rust parser dependency]**: Removing the cyvcf2 fallback means the pipeline fails hard if `stats_core.so` is missing or incompatible. → **Mitigation**: Add a clear startup check that verifies the Rust extension loads; fail with an actionable error message.
- **[Classification distribution shift]**: Fixing the classifier will change the distribution from 73% noise / 27% true_multi_allelic to closer to 60% normalization_artifact / 35% noise / 5% true_multi_allelic. Flag counts will change significantly. → **Mitigation**: This is the intended behavior; the current distribution is wrong.
