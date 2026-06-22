## Why

The refactor-variant-evidence-foundation change added new data columns (`modality_evidence_caller`, multi-allelic metrics, biological flags) that currently have no visual representation. The visualization code has accumulated technical debt: text marks on 21 bar chart functions lack overlap avoidance and auto-contrast coloring, `_add_bar_labels` is dead code, rescue charts reference the deprecated `RESCUED` column, `plot_dna_vs_rna_per_caller` only shows VAF (not DP), and the variant-category wise lacks per-caller DNA-vs-RNA comparison charts. Additionally, BAM alignment stats are limited to 6 metrics with no sample-wise faceted visualization.

## What Changes

### 1. Rescue Chart Migration to modality_evidence
- Update `plot_cross_modality` → `plot_modality_evidence` with 4-category stacked bar (cross_modality, dna_confident, rna_rescued, low_confidence)
- Update all 9 rescue chart functions to consume `modality_evidence_caller` as primary grouping column
- Add 4-category color scale for modality evidence categories
- Drop redundant CROSS_MODALITY panel (always identical to RESCUED)

### 2. DNA vs RNA Per-Caller VAF+DP Dual Scatters
- Extend `plot_dna_vs_rna_per_caller` to generate both VAF and DP scatter plots per caller
- Layout: 3 callers × 2 metrics = 6 panels (VAF row + DP row per caller, or side-by-side per metric)
- DP scatters clipped at 2000 (consistent with existing `plot_dna_vs_rna_dp`)
- Register in variant-category wise (already done in Proposal 1) and disease wise

### 3. Multi-Allelic Heterogeneity Visualization (NEW)
- `plot_multiallelic_classification`: Bar chart showing distribution of multiallelic_class categories per set
- `plot_allele_balance_scatter`: Scatter of major vs minor allele VAF at multi-allelic sites
- `plot_vaf_sum_histogram`: Histogram of vaf_sum across multi-allelic sites
- `plot_category_conflict_summary`: Bar chart of conflict types (Germline+NoConsensus, Somatic+NoConsensus, etc.)

### 4. Text Mark Optimization
- **BREAKING**: Remove `_add_bar_labels` (dead code, never called)
- Implement auto-contrast text color in `_make_bar_text` (white on dark bars, black on light bars, matching `_add_heatmap_text` pattern)
- Add responsive font sizing: scale `fontSize` based on bar count (fewer bars → larger text, more bars → smaller text)
- Add overlap strategy parameter: "hide" (default, skip labels on narrow bars), "rotate" (90° rotation for dense charts), "stagger" (alternating dy offsets)
- Unify label format across all bar charts: count via SI prefix (`~s` or `_si()` helper), percentage via `_si() + " (pct%)"` when `show_pct=True`
- Ensure all 21 bar chart functions use the optimized helper

### 5. BAM Alignment Stats Expansion
- Extract 4 new BAM metrics in `bam_stats.py`: duplication rate, coverage distribution bins (1×, 10×, 20×, 50×, 100×), properly paired fraction, insert size standard deviation
- Add `plot_bam_metrics_sample_wise`: faceted bar chart of BAM metrics per sample with set_number grouping
- Add `plot_bam_coverage_distribution`: coverage bin distribution as grouped bar chart

### 6. Dead Code Removal
- Remove `_add_bar_labels()` from `visualizer.py` (defined lines 255-264, never referenced)
- Remove any remaining references to deprecated `_RESCUE_COLOR_REGISTRY` (replaced by `_COLOR_REGISTRY["modality_evidence_caller"]` in Proposal 1)

## Capabilities

### New Capabilities

- `multi-allelic-visualization`: Charts for multi-allelic site classification, allele balance, VAF sum histogram, and category conflict summary. Consumes columns from Proposal 1 (`multiallelic_class`, `vaf_sum`, `allele_balance_ratio`, `category_conflict`).
- `per-caller-vaf-dp-scatters`: DNA vs RNA per-caller scatter plots with both VAF and DP metrics, 6-panel layout across 3 callers.
- `bam-alignment-visualization`: Sample-wise faceted BAM metrics and coverage distribution charts from expanded BAM statistics.
- `bar-text-optimization`: Unified text mark rendering with auto-contrast color, responsive sizing, overlap strategy, and consistent SI-prefix formatting.

### Modified Capabilities

- `variant-visualization`: Rescue charts migrated to `modality_evidence_caller`, `plot_cross_modality` redesigned, `_add_bar_labels` removed, `_make_bar_text` enhanced. Per-caller VAF scatters extended to VAF+DP. Multi-allelic and BAM charts added.
- `bam-statistics`: 4 new metrics extracted (duplication rate, coverage bins, properly paired fraction, insert size stddev).

## Impact

- `bin/vcf_stats/seq2neo/visualizer.py` — ~250 lines: `_make_bar_text` rewrite, `_add_bar_labels` removal, `plot_modality_evidence` redesign, rescue chart `evidence_col` wiring, `plot_dna_vs_rna_per_caller` extension, 4 new multi-allelic chart functions, 2 new BAM chart functions
- `bin/vcf_stats/seq2neo/bam_stats.py` — ~60 lines: 4 new metric extraction functions
- `bin/vcf_stats/seq2neo/cli.py` — ~20 lines: chart registry updates for new charts, disease-wise per-caller registration
- `bin/vcf_stats/tests/test_seq2neo_stats.py` — ~150 lines: tests for new charts, text optimization, BAM metrics
