## Why

The multi-allelic variant analysis system was partially implemented per the archived `2026-06-22-refactor-variant-evidence-foundation` proposal, but the core logic is substantially wrong: the classification uses crude aggregate thresholds instead of the spec's per-allele criteria (producing an inverted distribution — 73% noise, 27% "true multi-allelic" vs the expected 35%/5%); `allele_balance_ratio` uses the wrong formula (`2nd_max/max` instead of `max/total`); two of three multi-allelic flags diverge from spec; pre-decomposition VCF data is parsed but never consumed (dead code); and the `--exclude-multiallelic-conflict` CLI flag targets the wrong flag. Additionally, the `modality_evidence_caller` classification uses `>=1 DNA AND >=1 RNA` for `cross_modality` instead of the tier-based C1-C7 mapping specified in the proposal.

## What Changes

- **Fix multi-allelic classification criteria**: Replace the aggregate `allele_balance_ratio < 0.03` / `vaf_sum < 0.15` / catch-all logic with the spec's per-allele criteria: normalization_artifact (different REF/ALT lengths OR one allele VAF≈0/DP≈0), noise (one allele with signal, others VAF≈0/DP≈0), true_multi_allelic (≥2 alleles with VAF≥0.01 AND DP≥5 AND different ALT bases).
- **Fix `allele_balance_ratio` formula**: Change from `alt_dp_2nd_max / alt_dp_max` to `max_alt_dp / total_alt_dp` (matching spec and design.md).
- **Fix `flag_multi_allelic_heterogeneity`**: Exclude `normalization_artifact` from the flag — only `true_multi_allelic` sites should set this flag.
- **Fix `flag_category_conflict`**: Add the "both categories biological" guard — only flag when both/all FILTER categories have biological significance (not just NoConsensus/Artifact). Add `category_conflict_resolution` column with `keep_germline`/`drop_both` resolution.
- **Fix `--exclude-multiallelic-conflict`**: Change from checking `flag_vaf_overflow` to checking `flag_category_conflict` as the proposal specifies.
- **Wire pre-decomposition data into metrics**: Add `pre_norm_*` columns to `_CROSS_SAMPLE_COLS` so they survive the cross-sample scan. Use per-allele AF from pre-norm registry in `vaf_sum` and `allele_balance_ratio` when available. Compute `gt_cooccurrence` from pre-norm GT data. Compute per-allele `strand_balance` from F1R2/F2R1 and add `flag_extreme_strand_bias`.
- **Unify `modality_evidence_caller` with tier-based mapping**: Replace the current `>=1 DNA AND >=1 RNA` definition for `cross_modality` with the C1-C7 tier-based mapping: `cross_modality` = C1 (≥2 DNA + ≥2 RNA concordant), `dna_confident` = C2 (≥2 DNA), `rna_rescued` = C3+C4 (≥2 RNA), `low_confidence` = C5-C7.
- **Fix multi-allelic visualizations**: Fix `allele_balance` plot to show per-allele VAF balance (not per-caller VAF max/min). Fix `category_conflict` plot to use per-sample flag data (not cross-sample group_by that produces nonsensical single-category "conflicts").
- **Fix `OLD_MULTIALLELIC` validation**: Use the extracted `OLD_MULTIALLELIC` field to validate pre-norm parsing against the rescue VCF.
- **Remove cyvcf2 fallback from seq2neo stats**: Remove the cyvcf2 import and fallback path in `caller_parser.py` and `rescue_parser.py` — the Rust parser is the only backend. Remove the `--vcf-parser` CLI flag choice.

## Capabilities

### New Capabilities

- `multi-allelic-analysis`: Corrected multi-allelic site detection, classification (normalization_artifact/noise/true_multi_allelic with per-allele criteria), per-site metrics (allele_balance_ratio as max/total, vaf_sum with pre-norm enrichment, gt_cooccurrence, strand_balance), and pre-decomposition data consumption.
- `biological-evidence-flags`: Corrected biological flags — flag_multi_allelic_heterogeneity (true_multi_allelic only), flag_category_conflict (with biological-significance guard and resolution column), flag_extreme_strand_bias (new), and unified modality_evidence_caller (tier-based C1-C7 mapping).

### Modified Capabilities

- `variant-tiering-stats`: Unify modality_evidence_caller to use C1-C7 tier-based mapping instead of raw caller counts. Export FINAL_TIER_ORDER for visualizer consumption.

## Impact

- `bin/vcf_stats/seq2neo/statistics.py` — ~150 lines: fix `compute_multi_allelic_metrics` classification (~40 lines), fix `allele_balance_ratio` formula (~15 lines), fix `compute_biological_flags` (~30 lines), fix `compute_modality_evidence` (~25 lines), add strand_balance + gt_cooccurrence computation (~30 lines), add pre-norm column references (~10 lines)
- `bin/vcf_stats/seq2neo/caller_parser.py` — ~30 lines: remove cyvcf2 import and fallback path, ensure `parse_pre_norm_multiallelic` output columns are compatible with cross-sample scan
- `bin/vcf_stats/seq2neo/rescue_parser.py` — ~10 lines: remove cyvcf2 import and Python fallback
- `bin/vcf_stats/seq2neo/cli.py` — ~30 lines: add pre_norm columns to _CROSS_SAMPLE_COLS, fix `--exclude-multiallelic-conflict` to check `flag_category_conflict`, remove `--vcf-parser` flag, add `OLD_MULTIALLELIC` validation
- `bin/vcf_stats/seq2neo/visualizer.py` — ~40 lines: fix `plot_allele_balance_scatter` to use per-allele VAF, fix `plot_category_conflict_summary` to use per-sample data
- `bin/vcf_stats/seq2neo/tiering_stats.py` — ~15 lines: export FINAL_TIER_ORDER, update modality_evidence mapping
- Parquet schema: +3 new columns (gt_cooccurrence, strand_balance, flag_extreme_strand_bias, category_conflict_resolution), 1 column definition changed (allele_balance_ratio formula, modality_evidence_caller mapping)
- **BREAKING**: `modality_evidence_caller` semantics change — variants previously classified as `cross_modality` with only 1 DNA + 1 RNA caller will be reclassified as `low_confidence` or `rna_rescued`. Pipeline re-run required for correct modality stats.
