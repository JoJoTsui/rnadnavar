## Context

The seq2neo statistics pipeline processes rescue VCFs through Rust parsers, computes per-sample and cross-sample statistics, and generates ~35 altair charts. The first successful 64-sample run surfaced bugs and gaps that escaped testing with small sample sets. The current architecture is mature — all Rust paths are active, pysam fallbacks removed — so fixes are surgical: bug corrections, missing feature additions, and visualization quality improvements.

## Goals / Non-Goals

**Goals:**
- Fix 3 P0 bugs (BAM column swap, empty coverage violin, missing BAM integrity check)
- Add depth threshold sweep (total/REF/ALT DP) mirroring the VAF sweep pattern
- Add DP distribution plots (per-caller, BAM pileup, per-tier)
- Apply chromosome ordering to all CHROM-grouped charts consistently
- Fix VAF distribution plot styling (palette, axis clamping, threshold lines)
- Grid sample-wise charts by set_number for cross-set comparison
- Fix database enrichment chart (skip or annotate D=0 tiers)
- Add data leakage exclusion flags for downstream ML

**Non-Goals:**
- Refactoring the chart architecture (wise-based system stays as-is)
- Adding new chart types beyond DP distributions
- Changing the VCF parsing or tiering engine
- Modifying the Rust stats_core module
- Per-read insert size distribution (deferred to future proposal)

## Decisions

### D1: BAM integrity check — pre-processing, not per-BAM
**Choice:** Validate all BAMs upfront before any processing, by checking the BGZF EOF marker (magic bytes `1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00`) on the last 28 bytes of each file.
**Rationale:** Pre-processing check fails fast — broken files halt before hours of variant processing. The EOF check is O(1) per BAM (read last 28 bytes) and doesn't require parsing.
**Alternative:** Lazy check during pileup. Rejected — fails after expensive VCF parsing.

### D2: Depth threshold sweep — replicate VAF sweep pattern
**Choice:** Add `compute_dp_threshold_sweep()` in `statistics.py` that iterates per-caller DP columns (`{caller}_DP`) and computes retention % at each threshold. Thresholds: `[1, 2, 5, 10, 20, 50, 100, 200]` for total DP. For REF and ALT DP: `[0, 1, 2, 5, 10, 20, 50]`. Return a long-format DataFrame for plotting.
**Rationale:** The VAF sweep pattern is already proven, well-tested, and integrates with the existing chart system. Nothing new to invent.

### D3: DP distribution plots — reuse existing box+violin pattern
**Choice:** Call `plot_dp_distribution` (already exists for per-caller DP) and add `plot_bam_dp_distribution` (BAM pileup total/REF/ALT DP per BAM type) and `plot_per_tier_dp_boxplot` (per-tier DNA DP). All reuse the `_box_violin_overlay` helper.
**Rationale:** The `_box_violin_overlay` helper at `visualizer.py:163` already handles the combined chart generation. New functions are thin wrappers that melt columns and delegate.

### D4: Strelka VAF > 1 — clamp for viz, don't modify data
**Choice:** Add `vaf_display` column computed as `min(VAF, 1.0)` for visualization purposes. Add annotation text to Strelka charts. Do NOT modify the raw `{caller}_VAF` column.
**Rationale:** The caller VAF column is ground truth (Strelka uses tier-1 depth denominator — correct by design). Clamping is a visual convenience, not a data correction.

### D5: Chromosome ordering — add `_sort_chromosomes` calls, don't refactor
**Choice:** Add `_sort_chromosomes()` calls to the remaining chart functions that group by CHROM: `plot_vc_distribution`, `plot_dna_vs_rna_per_caller`, and any statistics aggregation with CHROM ordering. Do not change the helper signature.
**Rationale:** The helper (`visualizer.py:84`) already works correctly. The fix is just calling it consistently.

### D6: Data leakage exclusions — write filtered parquet, don't remove data
**Choice:** When `--exclude-disease`, `--min-vaf`, or `--min-dp` are specified, write additional filtered parquet files (`variant_details_filtered/`) alongside the full dataset. Charts and CSVs are generated from the full dataset only. Filtered files are for downstream ML consumption.
**Rationale:** Separating filtered output from analysis output prevents accidental data leakage while keeping the full pipeline output intact for review.

## Risks / Trade-offs

- **[Risk]** BAM integrity check adds ~1s per sample (reading last 28 bytes of 3 BAMs) → negligible vs hours of processing
- **[Risk]** DP threshold sweep adds ~40 lines to statistics.py → minimal complexity, well-tested pattern
- **[Risk]** Clamping Strelka VAF hides the denominator difference → mitigated by annotation text on the chart
- **[Risk]** Data leakage flags could be confusing (filtered vs full) → mitigated by writing to separate `variant_details_filtered/` directory
