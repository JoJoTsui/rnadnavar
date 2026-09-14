## Context

The seq2neo variant statistics pipeline processes 107.9M variants across 61 samples. After excluding NoConsensus (98.2M), the filtered dataset is 9.75M variants. However, 46.9% of these filtered variants (4.57M) are low-confidence:

- **3.31M Reference-with-signal** (FILTER=Reference but VAF>0.05) — the Reference label disagrees with the quantitative alt allele fraction. 89.8% are in C2D0 (≥2 DNA callers voted Reference) with mean ALT_DP=12.5, meaning the alt is real but the label is wrong.
- **1.53M low-evidence tiers** (C5/C6) — single-caller or no-caller support. C5D0 is 99.9% Artifact; C6D0 is 99.9% Artifact. Only 622 Somatic variants in C5/C6 (single-caller somatic, automatically suspect).
- **189K category conflicts** — different FILTER categories at the same (CHROM, POS), indicating caller disagreement.
- **35K Germline with VAF<0.10** — inconsistent with heterozygous expectation (~0.50). 98% have gnomAD AF (mean 0.204), suggesting common-variant false calls at low DNA VAF.
- **2K Somatic with VAF>0.60** — possible LOH, homozygous somatic, or misclassified germline.

The strongest under-utilized artifact discriminator is **strand bias**: `*_SB` (Mutect2 4-element F1R2/F2R1), `pre_norm_*_f1r2/f2r1` (per-allele), and `BAM_{DT,RT}_F1R2*/F2R1*` (per-position from pileup) are all in the parquet but drive ZERO flags or filters. Cross-sample recurrence is not computed at all — a "Somatic" variant found in 50/61 samples is almost certainly a mapping artifact, but no metric captures this.

The existing `build_unified_filter` (`cli.py:175`) handles `--include-filters`, `--exclude-filters`, `--min-vaf`, `--min-dp`, `--max-gnomad-af`, `--min-evidence-tier`, `--exclude-multiallelic-conflict` but does NOT use any biological flag, strand bias, recurrence, or paradoxical VAF logic.

## Goals / Non-Goals

**Goals:**
- Compute 6 new metrics (recurrence, strand bias, 4 new flags) that distinguish true from false variants
- Implement a 3-stage filtering strategy (hard/soft/confidence) that handles all three low-confidence categories
- Emit both filtered and unfiltered statistics so users can compare
- Add a real caller overlap matrix + support distribution
- Wire existing but unused strand bias and BAM per-position metrics into flags

**Non-Goals:**
- Fixing the Rust BAM stats bugs (MAPQ, coverage, expanded metrics) — that belongs in `fix-bam-stats-rust` (this proposal depends on it)
- Fixing multi-allelic classification — that belongs in `fix-multiallelic-logic` (this proposal depends on it)
- Machine learning-based variant classification — the 3-stage rule-based filter is sufficient for QC
- Tumor purity/ploidy estimation — requires external tools
- Changing the upstream variant calling pipeline

## Decisions

### Decision 1: Cross-sample recurrence via window function

**Chosen**: Compute `n_recurrent_samples` using polars `pl.len().over(['CHROM','POS','REF','ALT'])` on `combined_df` after all samples are scanned. This is a single window operation that adds a per-variant column without a separate group_by+join.

**Alternatives considered**:
- Separate group_by + join → Rejected: more memory, two-pass operation
- Compute at per-sample parquet write time → Rejected: can't know cross-sample recurrence at write time

### Decision 2: Strand bias via Fisher exact test

**Chosen**: Compute `max_strand_bias_fisher_p` as the minimum (most significant) Fisher exact p-value across all callers that have F1R2/F2R1 data. For each caller, construct a 2×2 contingency table: `[[F1R2_ref, F2R1_ref], [F1R2_alt, F2R1_alt]]` and compute Fisher exact p-value using `scipy.stats.fisher_exact`. The maximum across callers is taken because any single caller with extreme bias is sufficient to flag. When no F1R2/F2R1 data is available, the column is null.

**Alternatives considered**:
- Use only Mutect2 SB field → Rejected: misses DeepSomatic and pre-norm strand data
- Simple ratio threshold (F1R2/(F1R2+F2R1) < 0.10) → Rejected: doesn't account for read count; 2/48 and 20/480 have the same ratio but very different significance
- G-test instead of Fisher → Rejected: Fisher is exact for small counts, which is common for alt alleles

### Decision 3: Stage 0 hard filters use combined modality evidence

**Chosen**: Hard exclusion conditions use BOTH DNA and RNA evidence where applicable. For example, `(DNA_DP_mean < 5 AND RNA_DP_mean < 5)` excludes variants where neither modality has coverage, not just DNA-low variants. This prevents dropping RNA-rescued variants that have low DNA depth but adequate RNA depth.

**Alternatives considered**:
- DNA-only hard filters → Rejected: would drop legitimate RNA-rescued variants
- Separate DNA and RNA filter stages → Rejected: adds complexity; a single combined condition is clearer

### Decision 4: Confidence tiers map from final_tier + soft flag presence

**Chosen**: Confidence tier assignment is a simple lookup:
- HIGH: no soft flags AND final_tier ∈ {C1D0, C1D1, C2D1, C3D1, C4D1} AND (N_DNA_CALLERS_SUPPORT ≥ 1 OR N_RNA_CALLERS_SUPPORT ≥ 2)
- MEDIUM: no soft flags AND final_tier ∈ {C2D0, C3D0, C4D0}
- LOW: passes hard filters but has any soft flag
- DROPPED: fails hard filters (not in the output parquet, written to a separate `failed.tsv` with reason column)

**Alternatives considered**:
- Continuous confidence score → Rejected: less interpretable; users need discrete tiers for ML train/test splits
- Only HIGH/LOW (no MEDIUM) → Rejected: C2D0 (4.9M variants, mostly Reference) is not HIGH confidence but not LOW either

### Decision 5: Both filtered and unfiltered statistics emitted

**Chosen**: When the unified filter is active, compute and write both versions:
- `sample_summary.tsv` (unfiltered, 107.9M) + `sample_summary_filtered.tsv` (filtered, 9.75M)
- `set_summary.tsv` (unfiltered) + `set_summary_filtered.tsv` (filtered)
- `dataset_summary.tsv` (filtered, current behavior) + `dataset_summary_unfiltered.tsv` (unfiltered)
- `rescue_validation_summary.tsv` (currently unfiltered) — add a filtered version

The unfiltered versions are computed from `combined_df` before `.filter()`. The filtered versions are computed from the filtered `combined_df` after `.filter()`.

**Alternatives considered**:
- Only filtered → Rejected: user needs to see filtering impact; current behavior already has unfiltered sample/set summaries by accident (computed before filter)
- Only unfiltered → Rejected: the whole point of filtering is to produce clean downstream data

## Risks / Trade-offs

- **[Dependency on other proposals]**: `flag_low_rna_mapq` requires the MAPQ=255 fix from `fix-bam-stats-rust`. Stage 0 noise exclusion requires the corrected multi-allelic classification from `fix-multiallelic-logic`. → **Mitigation**: Implement in dependency order; `fix-bam-stats-rust` first, then `fix-multiallelic-logic`, then this proposal.
- **[Hard filter impact]**: Dropping ~1.45M variants (15%) may remove some true positives. The Reference-with-signal drop (VAF>0.10 AND ALT_DP≥3, ~900K) is the most aggressive. → **Mitigation**: Hard filter conditions are configurable via CLI flags; `--no-hard-filter` escape hatch preserves all variants with soft flags only.
- **[Strand bias computation cost]**: Fisher exact test on millions of variants adds compute time. → **Mitigation**: Only compute for variants with alt_reads ≥ 10 (skip low-alt variants where strand bias is not meaningful); use scipy's vectorized fisher_exact where possible.
- **[Recurrence computation memory]**: `pl.len().over([...])` on 9.75M variants requires a window operation. → **Mitigation**: polars handles this efficiently with streaming; the window key is (CHROM, POS, REF, ALT) which has high cardinality but is a standard group-by operation.
- **[Fisher exact dependency]**: Requires `scipy` which may not be installed. → **Mitigation**: Add scipy as a dependency or implement a pure-Python Fisher exact for small tables (2×2 is tractable).
