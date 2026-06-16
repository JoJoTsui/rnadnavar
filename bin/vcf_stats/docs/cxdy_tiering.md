# CxDy Variant Tiering System

## Overview

The CxDy tiering system classifies rescue VCF variants into evidence-based tiers
for downstream filtering and analysis. Tiers are computed by the Rust
`stats_core` tiering engine (`tier.rs`) via PyO3.

## Tier Format

Each variant is assigned a **final_tier** in `CxDy` format:
- **C** = Caller evidence level (C1 strongest, C7 weakest)
- **D** = Database evidence level (D1 = present in COSMIC/gnomAD, D0 = absent)

Example: `C1D1` = strong caller evidence + database support

## Caller Tiers (C1–C7)

| Tier | Description |
|------|-------------|
| C1 | ≥3 callers support the variant across both modalities |
| C2 | ≥2 callers support with strong concordance |
| C3 | ≥2 callers support |
| C4 | Single caller support with high confidence |
| C5 | Single caller support |
| C6 | Low-confidence single caller |
| C7 | No consensus (callers disagree) |

The exact count thresholds depend on the number of available callers per sample.

## Database Tiers (D0–D1)

| Tier | Description |
|------|-------------|
| D1 | Variant has COSMIC_ID or GNOMAD_AF annotation — external database evidence supports somatic relevance |
| D0 | No database annotation — variant is novel or not in curated databases |

## Implications for Downstream Use

- **C1D1**: Highest confidence — suitable for training sets
- **C1D0, C2D1**: High confidence — suitable for discovery
- **C3–C5**: Moderate confidence — requires manual review
- **C6–C7**: Low confidence — likely artifacts or germline
- **D=0**: No database evidence by design — these variants should NOT be considered
  "failed" annotations. They are simply novel or uncurated.
