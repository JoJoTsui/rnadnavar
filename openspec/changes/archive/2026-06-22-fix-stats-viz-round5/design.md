## Context

Round 5 addresses 23 requirements from pipeline review. The core insight: most issues are NOT individual bugs but **missing infrastructure**. The codebase has ~50 chart functions that each independently handle colors, faceting, text, clipping, and scales — leading to inconsistency. Building 6 centralized helpers first, then applying them systematically, fixes all issues at once.

## Goals / Non-Goals

**Goals:**
- Build helper infrastructure that makes ALL chart functions consistent by construction
- Unify colors for 16 categorical entities across ~40+ callsites
- Fix faceting for 17 functions (11 using wrong `alt.Column()` pattern + 6 needing grid sizing)
- Replace broken violin density plots with boxplots
- Add per-FILTER sub-plots and variant-category-wise dimension
- Add text to all 7 heatmaps and stacked bars
- Apply symlog to all count axes

**Non-Goals:**
- Refactoring the wise registry architecture itself (incremental addition only)
- Switching from altair to another charting library
- The streaming/OOM architecture (separate change)
- Automatic variant filtering or ML training

## Decisions

### D1: Helper-first approach — build infrastructure before fixing charts

**Choice:** Create 6 helper functions first, then migrate all chart functions to use them. This ensures consistency by construction — every chart that calls `_color_scale("caller")` gets the same colors, period.

**Alternative:** Fix each chart individually. Rejected because it's error-prone (40+ callsites), fragile (colors drift over time), and slower.

### D2: Color registry — explicit DOMAIN + COLORS lists for 10 entities

**Choice:** Define pairs of `{ENTITY}_DOMAIN` (list of values) and `{ENTITY}_COLORS` (list of hex colors) at module top. The `_color_scale(entity)` helper returns the right `alt.Scale`. For entities with variable domains (disease, tier, chromosome), use `scheme="category10"` as fallback.

**Entities with explicit palettes:**
1. Classification/FILTER: existing `CLASSIFICATION_DOMAIN` + `CLASSIFICATION_COLORS`
2. Caller: new `CALLER_DOMAIN` + `CALLER_COLORS` — DNA/RNA paired (blue/light-blue, orange/light-orange, green/light-green)
3. Rescue: new `RESCUE_DOMAIN` + `RESCUE_COLORS` — YES=green, NO=red
4. Variant type: new `VARIANT_TYPE_DOMAIN` + `VARIANT_TYPE_COLORS`
5. Somatic modality: new `SOMATIC_MODALITY_DOMAIN` + `SOMATIC_MODALITY_COLORS`
6. BAM type: new `BAM_TYPE_DOMAIN` + `BAM_TYPE_COLORS`
7. Modality DNA/RNA: new `MODALITY_DOMAIN` + `MODALITY_COLORS`
8. Agreement level: new `AGREEMENT_DOMAIN` + `AGREEMENT_COLORS`

**Entities with dynamic domains (use scheme fallback):**
9. Set number, 10. Disease, 11. Tier, 12. Chromosome, 13. N_SUPPORT_CALLERS, 14. REDI_EVIDENCE, 15. caller_tier, 16. BAM metric

### D3: Faceting — `_apply_faceting()` with per-group_col rules

**Choice:** A single helper that applies faceting based on `group_col` value:
- `set_number` → `.facet(columns=2).resolve_scale(x="independent", y="shared")`
- `disease_normalized` → `.facet(columns=4).resolve_scale(y="shared")`
- `sample_id` → `.facet(columns=2).resolve_scale(x="independent", y="shared")`
- `FILTER` → `.facet(columns=3).resolve_scale(y="shared")`
- Others → `.facet(columns=3)` default

**Critical:** `x="independent"` for per-sample charts so each facet shows ONLY its own samples, not all samples.

### D4: Violin replacement — boxplot only (no KDE)

**Choice:** Remove all `transform_density` + `mark_area` code. Use `mark_boxplot()` alone. KDE in Python (scipy) would add a heavy dependency for marginal visual benefit.

### D5: Per-FILTER sub-plots — facet within chart, not separate charts

**Choice:** Instead of calling chart functions 3× with pre-filtered data, add `FILTER` as a facet column within the chart. This:
- Shares axis scales for valid comparison
- Requires zero new function signatures
- Uses the same `_apply_faceting()` helper

For threshold sweeps: stop filtering to `classification.is_null()` and add `color=classification:N` encoding.

### D6: Symlog — `constant=1` for count axes

**Choice:** `alt.Scale(type="symlog", constant=1)` — linear for |x| < 1, logarithmic beyond. Handles zeros in stacked bars without `-inf`. Applied via `_count_scale()` helper.

### D7: FP cross-tab — include Somatic, rename

**Choice:** Remove `FILTER != "Somatic"` filter in `compute_fp_cross_tab()`. Rename chart title to "Classification × Caller Support Cross-Tabulation". The TP (Somatic) provides valuable contrast against FP patterns.

### D8: Resume fallback — recompute from parquet

**Choice:** When `--resume` and `sample_summary.tsv` doesn't exist, recompute by reading each parquet file and running `sample_summary()`. This is fast (per-sample aggregation) and prevents silent chart omission.

## Risks / Trade-offs

- **[Risk]** Migrating 17 functions to `_apply_faceting()` in one batch may introduce layout regressions → mitigate with visual spot-checks on test run
- **[Risk]** Removing violin overlay loses density shape information → boxplots still show median, IQR, whiskers, outliers — sufficient for most analysis
- **[Risk]** Including Somatic in FP cross-tab may dominate the color scale → use log color scale (already present)
- **[Trade-off]** `symlog` count scale makes reading exact values harder than linear → tooltips remain available for precise values
- **[Trade-off]** 6 new helpers increase visualizer.py complexity → but reduce overall line count through deduplication
