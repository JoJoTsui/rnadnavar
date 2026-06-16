# FP Cross-Tabulation

## Overview

The FP (false positive) cross-tabulation computes a FILTER x
N_SUPPORT_CALLERS matrix for all non-Somatic variants. This reveals
how many callers support each type of non-somatic classification,
helping to identify systematic caller agreement patterns for
false-positive categories.

## Interpretation

The cross-tab answers:

- **How many callers support Germline calls?** High caller support
  for Germline variants is expected (germline variants are present
  in all cells and should be detected by most callers).
- **Are Artifact calls single-caller?** Artifacts often come from
  one caller's specific error mode. A high count at N_SUPPORT=1
  confirms this.
- **Do NoConsensus variants have low support?** By definition,
  NoConsensus variants failed the consensus threshold, so they
  should cluster at low N_SUPPORT_CALLERS values.

## Output Files

- `stats/threshold/fp_cross_tab.tsv` — Long-form table with columns:
  FILTER, N_SUPPORT_CALLERS, count.
- `plots/45_fp_cross_tab_heatmap.html` — Heatmap visualization with
  log-scaled color encoding.

## Implementation

Function: `compute_fp_cross_tab()` in `statistics.py`

Filters to `FILTER != "Somatic"`, groups by `[FILTER, N_SUPPORT_CALLERS]`,
and counts variants per group.
