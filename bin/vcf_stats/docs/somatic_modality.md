# Somatic Modality Sub-Classification

## Overview

Somatic variants are sub-classified by their modality support pattern,
derived from the `caller_tier` column (CxDy tiers). This provides
insight into whether a somatic call is supported by both DNA and RNA
callers, or only one modality.

## Modality Categories

| Modality       | Caller Tiers | Description                          |
|----------------|-------------|--------------------------------------|
| MultiModality  | C1          | Both DNA and RNA callers agree       |
| DNA_only       | C2, C5      | Supported by DNA callers only        |
| RNA_only       | C3, C6      | Supported by RNA callers only        |
| Weak           | C4, C7      | Low caller support (single caller)   |
| Unknown        | Other       | Unrecognized tier value              |

## Interpretation

- **MultiModality** variants have the highest confidence — independent
  confirmation from both sequencing modalities.
- **DNA_only** variants may represent true somatic mutations not expressed
  in RNA, or low-expression variants below RNA detection threshold.
- **RNA_only** variants may include RNA editing events misclassified as
  somatic, or variants in regions with poor DNA coverage.
- **Weak** variants have minimal caller support and warrant manual review.

## Output Files

- `stats/tier/somatic_modality.tsv` — Per-modality summary with variant
  counts and mean VAF/DP.
- `plots/46_somatic_modality_pie.html` — Pie chart of modality distribution.
- `plots/47_somatic_modality_bars.html` — Bar chart with variant counts
  and mean DNA VAF by modality.

## Implementation

Function: `compute_somatic_modality()` in `statistics.py`

The function filters to `FILTER == "Somatic"` variants, maps `caller_tier`
to modality categories, then aggregates counts and mean VAF/DP per modality.
