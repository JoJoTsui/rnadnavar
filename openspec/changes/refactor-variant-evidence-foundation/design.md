# Design: Refactor Variant Evidence Foundation

## Context

The seq2neo variant statistics pipeline (`bin/vcf_stats/seq2neo/`) processes rescue VCFs and individual caller VCFs through a Rust-accelerated pipeline: parse → join → tier → pileup → write per-sample parquet → scan all parquets → compute cross-sample statistics → generate Altair visualizations. The current design has four architectural gaps:

1. **Duplicate columns**: `RESCUED` and `CROSS_MODALITY` are always set to the same boolean (`in_dna and in_rna` in `bin/vcf_utils/tagging.py:143-148`), yet both are stored, counted, and charted.

2. **Lost multi-allelic data**: The pipeline only parses decomposed/normalized VCFs (`*.dec.norm.vcf.gz`), where vt decompose strips FORMAT to GT:DP, discarding full AD arrays, per-allele AF, and strand bias. Pre-decomposition VCFs exist at `variant_calling/{caller}/.../*.vcf.gz` with 794/4348 (18.3%) multi-allelic Mutect2 records containing full data.

3. **ASE-confound**: `cli.py:1414` uses `RNA_VAF_mean >= min_vaf` as a filtering criterion, but RNA VAF is confounded by allele-specific expression — a true somatic variant can appear at RNA VAF=0.01 due to expression silencing of the variant allele.

4. **Disconnected filters**: `--pileup-mode` (BAM pileup only), `--min-vaf/--min-dp` (filtered parquet only), and `--exclude-disease` (variant_details_filtered/ only) operate at three different pipeline stages with no unified application.

## Goals / Non-Goals

**Goals:**
- Parse pre-decomposition caller VCFs to recover full multi-allelic AD arrays, per-allele AF, and per-allele strand bias
- Replace binary RESCUED/CROSS_MODALITY with three-way modality_evidence computed at process_single_sample() time
- Add biological evidence flags that distinguish real biology from technical artifacts
- Remove RNA_VAF from filtering decisions; use RNA caller support + RNA depth instead
- Build a single unified filter pipeline at the combined_df injection point
- Maintain backward compatibility through deprecated columns and flag aliases

**Non-Goals:**
- Visualization changes (text marks, per-caller VAF+DP scatters, multi-allelic charts) — these belong in Proposal 2
- BAM alignment stats expansion — Proposal 2
- Changing the upstream Nextflow pipeline or VCF normalization — the pre-decomposition VCFs are read as-is
- Adding tumor purity/ploidy estimation — requires external tools and manifest schema changes
- Integrating PyClone or other subclonal reconstruction tools

## Decisions

### Decision 1: Pre-decomposition VCFs parsed in new optional function, not a separate pipeline

**Chosen**: Add `parse_pre_norm_multiallelic()` to `caller_parser.py`, called from `process_single_sample()` as an optional enrichment step that produces a polars DataFrame joined with the normalized parse.

**Alternatives considered**:
- Run a separate pipeline stage before the main pipeline → Rejected: adds I/O overhead, complicates orchestration, would need to write intermediate files
- Replace normalized VCF parsing entirely with pre-norm parsing → Rejected: pre-norm VCFs lack normalization (indels may be misaligned), don't exist for all samples, Strelka is always biallelic

**Rationale**: The enrichment pattern (parse extra data, join in-place) is consistent with the existing BAM pileup enrichment in `process_single_sample()`. The function returns a polars DataFrame keyed on (CHROM, POS, REF, ALT) for direct `.join()`.

### Decision 2: modality_evidence computed as two scheme columns

**Chosen**: Compute both `modality_evidence_caller` (based on N_DNA_CALLERS_SUPPORT / N_RNA_CALLERS_SUPPORT) and `modality_evidence_dp` (based on DNA_DP_mean / RNA_DP_mean > 30) as separate columns.

**Alternatives considered**:
- Single column using caller support only → Rejected: user requested both caller-based and DP-based definitions
- Single column that combines both → Rejected: conflates two independent evidence dimensions; users should choose which definition to filter on

**Rationale**: Two columns allow filtering on either definition independently. The caller-based scheme is more robust (3 orthogonal callers agree = strong evidence). The DP-based scheme is useful when caller support counts are unreliable (e.g., Strelka-only called variants).

### Decision 3: Unified filter at the LazyFrame level

**Chosen**: Apply `.filter(filter_expr)` on `combined_df` immediately after `pl.scan_parquet()`. The filter expression is a composition of AND-ed boolean conditions built from CLI flags. Because `combined_df` is a LazyFrame, polars pushes the filter into the parquet scan.

**Alternatives considered**:
- Apply filter at per-sample parquet write time → Rejected: requires re-running per-sample processing; doesn't work with --resume which reloads existing parquet files
- Apply filter in each downstream function individually → Rejected: error-prone, inconsistent, duplicates filter logic
- Apply filter on the eager DataFrame after `.collect()` → Rejected: defeats the purpose of lazy scanning — would materialize all columns before filtering

**Rationale**: Single injection point, zero overhead (predicate pushdown), automatic propagation to all downstream consumers. The `--no-filter` escape hatch skips the filter entirely.

### Decision 4: RNA_VAF removed from --min-vaf, replaced with caller support + depth

**Chosen**: The RNA leg of the VAF filter becomes `(N_RNA_CALLERS_SUPPORT >= 2) & (RNA_DP_mean >= 10)` instead of `(RNA_VAF_mean >= threshold)`.

**Alternatives considered**:
- Keep RNA_VAF with a warning → Rejected: doesn't fix the problem; users may ignore warnings
- Remove RNA entirely from the filter → Rejected: RNA caller evidence IS valid and useful; discarding it entirely loses valuable orthogonal confirmation
- Normalize RNA VAF by some factor to "correct" for ASE → Rejected: ASE is gene-specific and sample-specific; no universal correction factor exists

**Rationale**: RNA caller agreement (N_RNA_CALLERS_SUPPORT) tells us "is the variant real?" — this is not confounded by ASE. RNA depth tells us "do we have enough RNA data?" — also not confounded. Together they provide valid orthogonal evidence without the ASE confound.

### Decision 5: Biological flags are additive annotations, not filters

**Chosen**: Biological flags are computed as new DataFrame columns. The unified filter CAN reference them (e.g., `--exclude-multiallelic-conflict` checks `flag_category_conflict`), but by default they only annotate. Flagged variants survive filtering unless explicitly excluded.

**Alternatives considered**:
- Hard-filter all flagged variants → Rejected: many flags indicate "review this" not "this is wrong" (e.g., germline_low_vaf might be a real mosaic variant)
- Don't add flags at all, just filter → Rejected: loses the opportunity to annotate variants for downstream review

**Rationale**: Flags provide transparency — users can see WHY a variant was flagged and decide whether to filter. The `--no-filter` mode with flags set provides a complete picture of variant evidence quality.

## Key Data Flow

```
process_single_sample():
  1. Parse rescue VCF (normalized) → df
  2. Parse pre-norm caller VCFs → allele_registry (NEW)
  3. Join caller data to df (existing) 
  4. Join allele_registry to df on (CHROM, POS, REF, ALT) (NEW)
  5. Add sample metadata
  6. compute_all_per_variant(df) — existing VAF/DP means
  7. compute_tiers_for_dataframe(df) — existing tiering
  8. compute_modality_evidence(df) — NEW: modality_evidence_caller + modality_evidence_dp
  9. compute_biological_flags(df) — NEW: 8 flag columns
  10. compute_multi_allelic_metrics(df) — NEW: position-grouped metrics
  11. BAM pileup (existing, with positions from filtered df)
  12. Write per-sample parquet (existing, with new columns)

main():
  combined_df = pl.scan_parquet(variant_dir / "*_variants.parquet")
  filter_expr = build_unified_filter(args) — NEW
  if not args.no_filter:
      combined_df = combined_df.filter(filter_expr) — NEW: single injection point
  # ALL downstream functions receive filtered combined_df
  disease_summary(combined_df)
  compute_tier_summary(combined_df)
  dataset_summary(combined_df)
  # ... all 20+ statistics functions
  # ... all 40+ chart functions
```

## Column Changes to Per-Sample Parquet

**New columns (+15)**:
| Column | Type | Source |
|--------|------|--------|
| `modality_evidence_caller` | String | Computed from N_DNA/RNA_CALLERS_SUPPORT |
| `modality_evidence_dp` | String | Computed from DNA/RNA_DP_mean |
| `n_alleles_at_site` | UInt32 | (CHROM, POS) group_by |
| `vaf_sum` | Float64 | sum of DNA_VAF_mean per position |
| `allele_balance_ratio` | Float64 | max_alt_dp / total_alt_dp |
| `category_conflict` | Boolean | n_unique(FILTER) > 1 per position |
| `multiallelic_class` | String | "normalization_artifact" / "noise" / "true_multi_allelic" |
| `flag_vaf_overflow` | Boolean | vaf_sum > 1.1 |
| `flag_multi_allelic_heterogeneity` | Boolean | multiallelic_class == "true_multi_allelic" |
| `flag_category_conflict` | Boolean | category_conflict AND both categories are biological |
| `flag_germline_low_vaf` | Boolean | FILTER=Germline & VAF<0.10 & DP≥10 |
| `flag_somatic_high_vaf` | Boolean | FILTER=Somatic & VAF>0.60 |
| `flag_reference_with_signal` | Boolean | FILTER=Reference & VAF>0.05 |
| `flag_rna_rescued` | Boolean | DNA_VAF<0.05 & N_RNA≥2 & RNA_DP≥10 |
| `gt_cooccurrence` | String | From pre-norm GT (nullable) |

**Deprecated columns (retained, not in _CROSS_SAMPLE_COLS)**:
- `RESCUED` → maps to `modality_evidence_caller IN ("cross_modality", "rna_rescued")`
- `CROSS_MODALITY` → same mapping

**Added to _CROSS_SAMPLE_COLS**: `modality_evidence_caller`, `modality_evidence_dp`, all biological flag columns, `n_alleles_at_site`, `vaf_sum`, `multiallelic_class`

## CLI Flag Changes

**New flags**:
- `--include-filters FILTER...` — Only include variants with these FILTER values
- `--exclude-filters FILTER...` — Exclude variants with these FILTER values
- `--min-dna-callers N` — Minimum DNA callers required (default: 0, no filter)
- `--min-rna-callers N` — Minimum RNA callers required (default: 0, no filter)
- `--min-evidence-tier TIER` — Minimum modality_evidence tier
- `--max-gnomad-af FLOAT` — Maximum gnomAD allele frequency
- `--exclude-multiallelic-conflict` — Exclude variants at multi-allelic sites with category conflicts
- `--no-filter` — Disable all filtering (escape hatch)

**Deprecated flags (still accepted, emit warning)**:
- `--pileup-mode {all,filtered}` → Maps to `--exclude-filters NoConsensus` (for "filtered")
- The existing `--min-vaf` and `--min-dp` flags retain their names but their behavior changes for RNA (caller support + depth instead of VAF)

**Integrated flags (unchanged interface, now part of unified filter)**:
- `--min-vaf FLOAT` — Minimum VAF (DNA leg) or RNA caller support (RNA leg)
- `--min-dp INT` — Minimum DP (EITHER modality passes)
- `--exclude-disease DISEASE...` — Sample-level disease exclusion

## Risks / Trade-offs

- **[Pre-norm VCF availability]** — Pre-decomposition VCFs may not exist for all samples (e.g., if the pipeline was run with `--save_bam_mapped false` or output was cleaned up). → **Mitigation**: Pre-norm parsing is best-effort; if files are absent, the pipeline logs a warning and continues with null-filled multi-allelic columns.
- **[vt decompose format variability]** — Different versions of vt may write AD differently after decomposition. → **Mitigation**: The pre-norm parse is validated against OLD_MULTIALLELIC when available; mismatches are logged but don't halt the pipeline.
- **[Parquet schema drift]** — Adding 15 columns changes the parquet schema. `--resume` reloads parquet files; if they lack new columns, they must be recomputed. → **Mitigation**: `--resume` checks for the presence of `modality_evidence_caller` column; if absent, forces re-processing of that sample.
- **[Backward compatibility of RESCUED/CROSS_MODALITY]** — Downstream scripts may read these columns from the parquet directly. → **Mitigation**: Both columns are retained as deprecated aliases for one release cycle. They are excluded from `_CROSS_SAMPLE_COLS` to prevent accidental use in new code.
- **[Performance impact of pre-norm parsing]** — Parsing an additional VCF per caller adds I/O. Mutect2 pre-norm VCFs are small (~4,300 records per sample). → **Mitigation**: Pre-norm parsing is enabled only for Mutect2 and DeepSomatic callers; Strelka is skipped. The 18.3% multi-allelic records are a small fraction of the ~4,300-record VCF. I/O overhead is minimal compared to the existing caller VCF parsing (6 VCFs × tens of thousands of records).
- **[ASE correction may drop real RNA-rescued variants]** — If a variant has strong RNA caller support but low RNA_DP, it fails the new filter. → **Mitigation**: The RNA_DP≥10 threshold is a default, not hardcoded — future work could add `--min-rna-dp` flag for configurability.

## Migration Plan

1. **Phase 1 — Data model**: Add new columns in `process_single_sample()`, keep RESCUED/CROSS_MODALITY as deprecated aliases. Update `_CROSS_SAMPLE_COLS`. All existing tests should pass (new columns are additive).

2. **Phase 2 — Unified filter**: Implement `build_unified_filter()` in cli.py, apply at combined_df creation. Add new CLI flags. Deprecate `--pileup-mode`. Update CLI tests.

3. **Phase 3 — Statistics migration**: Update `_WISE_METRICS`, rescue statistics functions, and sample/dataset summary to use modality_evidence. Keep old metric names as deprecated aliases.

4. **Phase 4 — Verification**: Full pipeline re-run with `--resume` to verify new columns are computed. Run with `--no-filter` to verify parity with existing output. Run with various filter combinations to verify unified filter behavior.

Rollback: Remove the `.filter()` call for the unified filter (reverting to `--no-filter` behavior). The new columns are additive — old code that doesn't reference them continues to work. Deprecated columns are still present.

## Open Questions

1. **Should pre-norm VCF paths be added to the manifest CSV/parquet, or derived from `base_output_dir` + `dir_name`?** Derivation is simpler (no manifest regeneration needed) and consistent with how `CALLER_CONFIGS` already constructs normalized paths. But explicit paths in the manifest are more robust. **Current leaning**: Derive from `base_output_dir` (matching existing pattern for normalized paths).

2. **Should the `multiallelic` flag from the rescue VCF INFO be repaired (set correctly) or replaced by position-based detection?** The flag is currently always false. Position-based detection finds 3,500+ sites per sample. Pre-norm parsing provides the ground truth. **Current leaning**: Replace the broken flag with position-based detection + pre-norm validation. Don't attempt to fix the upstream VCF flag.

3. **DP threshold of 30 for modality_evidence_dp — is this the right value?** 30× is a standard high-quality coverage threshold in WES/WGS. But this may need to be sample-specific (some samples have lower mean coverage). **Current leaning**: Use 30 as default, make it configurable via `--modality-dp-threshold` in a follow-up.
