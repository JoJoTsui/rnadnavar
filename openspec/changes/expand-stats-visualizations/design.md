## Context

The `refactor-variant-evidence-foundation` change (Proposal 1) added 15 new data columns to per-sample parquet files: `modality_evidence_caller`, `modality_evidence_dp`, multi-allelic metrics, and biological flags. The visualization code in `visualizer.py` currently references the deprecated `RESCUED` column in rescue charts, `plot_dna_vs_rna_per_caller` only generates VAF scatters (no DP), and text marks on bar charts lack auto-contrast and overlap avoidance. The `_add_bar_labels` function is dead code. BAM stats are limited to 6 basic metrics.

This change is purely additive to the visualization layer — it consumes columns that already exist (from Proposal 1) and adds new charts. It does not modify the data model.

## Goals / Non-Goals

**Goals:**
- Wire rescue charts to use `modality_evidence_caller` (already computed in Proposal 1's parquet output)
- Extend `plot_dna_vs_rna_per_caller` to generate both VAF and DP scatters
- Add 4 new multi-allelic charts and 2 new BAM charts
- Optimize text marks: auto-contrast, responsive sizing, overlap strategy, unified formatting
- Remove dead `_add_bar_labels` code
- Expand BAM stats with 4 new metrics

**Non-Goals:**
- Changing the data model or per-sample parquet schema — columns already exist from Proposal 1
- Modifying the unified filter or CLI flags — already done in Proposal 1
- Adding charts for biological flags (flag_*) — left for a future change
- Tumor heterogeneity subclonal analysis — requires external tools beyond visualization

## Decisions

### Decision 1: `_make_bar_text` rewrite with luminance-based contrast

**Chosen**: Compute perceived luminance of the bar fill color using the sRGB luminance formula `L = 0.299*R + 0.587*G + 0.521*B`. If L < 0.5, use white text; otherwise use black text. When the bar fill is not known (no color encoding), use `#333` (dark gray).

**Alternatives considered**:
- Always black text → Rejected: invisible on dark bar segments (stacked bars)
- Always white text → Rejected: invisible on light bar segments (yellow, light blue)
- User-specified color per chart → Rejected: adds burden to chart authors, inconsistent

**Rationale**: Luminance-based contrast is used by `_add_heatmap_text` (lines 235-252) in the same file, proving the approach works. Consistency with existing heatmap pattern.

### Decision 2: Responsive font sizing via category count

**Chosen**: Auto-scale fontSize based on the number of unique x-axis categories, computed from the DataFrame before chart construction.

**Alternatives considered**:
- Fixed font size with manual override → Rejected: what works for 8 bars doesn't work for 60
- Altair's built-in responsive sizing → Rejected: Altair doesn't provide auto-font-scaling for mark_text
- Compute from bar pixel width → Rejected: requires knowing output dimensions at render time, not composable with Vega-Lite

**Rationale**: Category count is a reliable proxy for bar width. The mapping (≤10→10px, 11-25→9px, 26-50→8px, >50→7px) is conservative and works for typical 800-1200px chart widths.

### Decision 3: Overlap strategies

**Chosen**: Three strategies controlled by `overlap` parameter: `"hide"` (default, skip on dense charts), `"rotate"` (90° rotation), `"stagger"` (alternating dy). Default `"hide"` for backward compatibility and because hidden labels are better than overlapping labels.

**Rationale**: No single strategy works for all charts. Dense sample-wise charts (>60 samples) benefit from "hide". Tier-wise charts (~14 tiers) benefit from "stagger". Category charts with long labels benefit from "rotate".

### Decision 4: VAF+DP layout as 2-row × 3-column grid

**Chosen**: `alt.vconcat(alt.hconcat(vaf_row), alt.hconcat(dp_row))` — VAF row on top, DP row below, 3 caller columns each. Shared color scale via `.resolve_scale(color="shared")`.

**Alternatives considered**:
- Side-by-side per caller (VAF|DP for each caller in 3 columns) → Rejected: harder to compare VAF across callers
- 3-row × 2-column → Rejected: VAF vs DP comparison per caller is secondary to cross-caller comparison
- Single scatter with VAF on x and DP on y → Rejected: loses the DNA-vs-RNA dimension

**Rationale**: This layout groups by metric (VAF row, DP row) allowing easy comparison across callers for the same metric. The DNA-vs-RNA relationship is shown within each scatter.

## Risks / Trade-offs

- **[Auto-contrast accuracy]** — Luminance-based contrast works for most colors but fails for colorblind users who can't distinguish red/green. → **Mitigation**: The 4 modality_evidence colors were chosen to be colorblind-friendly (green, blue, orange, red).
- **[Text removal on dense charts]** — "hide" strategy silently drops labels from charts with many bars. Users may think labels are broken. → **Mitigation**: Print a warning once when labels are hidden: "Text labels hidden on chart X due to N categories (threshold: 30). Use overlap='rotate' to force labels."
- **[BAM stats expansion requires BAM reprocessing]** — New metrics need BAM files to be re-scanned. → **Mitigation**: Fall back to null values for existing `bam_stats.tsv` files. The `--resume` flag reloads bam_stats.tsv; if it lacks new columns, they are null-filled.
