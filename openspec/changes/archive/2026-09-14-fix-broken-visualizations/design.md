## Context

The seq2neo visualizer (`bin/vcf_stats/seq2neo/visualizer.py`, 3074 lines) generates ~40 unique chart types across 10 per-wise directories, plus a `dashboard.html` aggregating all charts. The current dashboard is 100MB and broken: `fig.to_html(inline=True)` embeds the full vega-embed JS library (~1MB) into each of 106 charts, and `_extract_body_content()` uses a regex that matches `<body>` strings inside the inlined JS library (in template literals for the "View Source" feature) rather than the real HTML `<body>` tag, causing raw JS fragments to appear as visible text in every chart div.

Additionally, 8 individual plot behaviors are incorrect: `coverage_distribution` renders empty when all cov columns are null; `metrics_sample_wise` uses `column=` encoding (single row) instead of grid faceting and doesn't facet by set; `20_bam_metrics` uses fixed `width=350` producing thin bars; `10_ti_tv_ratio` text layer is missing chromosome sort; disease-wise and set-wise registries don't pass `facet_col=` so 14/17 charts don't facet; tier-axis plots lack explicit `FINAL_TIER_ORDER` sort; and `_sample_if_large` uses `df.sample()` without a seed making output non-deterministic.

The `FINAL_TIER_ORDER` constant exists at `tiering_stats.py:67-68` (`["C1D1","C1D0","C2D1","C2D0",...,"C7D1","C7D0"]`) but is never imported by `visualizer.py`.

## Goals / Non-Goals

**Goals:**
- Fix dashboard.html to produce valid, renderable HTML under 1MB
- Fix all 8 broken/incorrect individual plot behaviors
- Make visualizations deterministic (seeded sampling)
- Make tier-axis sorting explicit and future-proof (C10+ safe)
- Make per-wise plots correctly facet by their wise dimension

**Non-Goals:**
- Redesigning the chart set or adding new chart types (that belongs in a separate proposal)
- Fixing the BAM stats data quality issues (cov_*_pct null, MAPQ=255) — those belong in `fix-bam-stats-rust`
- Fixing multi-allelic visualization logic (allele_balance, category_conflict) — that belongs in `fix-multiallelic-logic`
- Changing the Altair vs Plotly split between seq2neo and parent visualizer

## Decisions

### Decision 1: Dashboard uses `inline=False` with CDN script tags

**Chosen**: Replace `fig.to_html(inline=True)` with `fig.to_html(inline=False)` which produces `<script src="...vega-embed CDN...">` tags instead of inlining the library. The dashboard template includes one vega-embed/vega-lite/vega `<script>` block in `<head>`, and each chart div contains only its `<div id="vis-N">` + `<script>vegaEmbed('#vis-N', spec, ...)</script>` with the JSON spec.

**Alternatives considered**:
- Extract only the JSON spec via `fig.to_dict()` and build the vegaEmbed calls manually → Rejected: more code, harder to maintain, loses altair's HTML generation
- Keep `inline=True` but fix `_extract_body_content` regex → Rejected: still 100MB, the root cause is the per-chart library inlining

**Rationale**: `inline=False` is the standard altair pattern for multi-chart pages. Size drops from ~100MB to <1MB. The `_extract_body_content` function still needs fixing (to handle the `inline=False` output format correctly), but the regex problem with `<body>` inside JS library strings disappears because the library is no longer inlined.

### Decision 2: Fix `_extract_body_content` to skip `<script>` blocks

**Chosen**: Strip `<script>...</script>` blocks from the HTML before applying the `<body>` regex, OR search for `</head>\s*<body>` as the body-start marker. This is robust against `<body>` strings appearing inside JS template literals.

**Alternatives considered**:
- Use an HTML parser (BeautifulSoup) → Rejected: adds a dependency for a simple extraction task
- Search for the last `</head>` and first `</body>` → Rejected: still fragile if `</head>` appears in JS

### Decision 3: Wise registries pass `facet_col` explicitly

**Chosen**: Add `facet_col="disease_normalized"` to all disease-wise registry entries and `facet_col="set_number"` to all set-wise registry entries in `cli.py`. Chart functions that already call `_apply_faceting()` will use the passed value; functions that don't will receive it via the `facet_col` parameter.

**Alternatives considered**:
- Make all chart functions internally call `_apply_faceting(group_col)` → Rejected: some charts use group_col as x-axis, not as facet; making them all facet would break charts where group_col is the primary axis
- Add a `facet_col` parameter to every chart function signature → This is needed regardless; the registry must pass it

### Decision 4: Tier sort uses imported `FINAL_TIER_ORDER`

**Chosen**: Import `FINAL_TIER_ORDER` from `tiering_stats.py` at the top of `visualizer.py` and pass `sort=FINAL_TIER_ORDER` to all `alt.X("final_tier:N")` and `alt.Column("caller_tier:N")` encodings.

**Alternatives considered**:
- Define a separate tier order in visualizer.py → Rejected: DRY violation, could drift
- Use a natural sort plugin → Rejected: over-engineered for a fixed 14-element list

### Decision 5: Sampling seed is a module-level constant

**Chosen**: Add `SAMPLE_SEED = 42` as a module-level constant in `visualizer.py` and pass `seed=SAMPLE_SEED` to every `df.sample()` call in `_sample_if_large`.

**Alternatives considered**:
- CLI-configurable seed → Rejected: unnecessary complexity; determinism is the goal, not user control
- Random seed from system time → Rejected: defeats the purpose

## Risks / Trade-offs

- **[CDN dependency]**: `inline=False` requires internet access to load vega-embed from CDN. In air-gapped environments, charts won't render. → **Mitigation**: The individual `.html` files in `plots/` subdirectories already use `inline=False` (CDN links), so this is consistent with existing behavior. For air-gapped use, users can download the vega-embed JS locally and update the CDN path.
- **[Wise faceting with single group]**: When all samples have `set_number=0` or `disease_normalized=""`, faceting produces a single panel — same visual result as no faceting. No regression, but no benefit until multiple sets/diseases exist.
- **[Orphaned functions]**: Wiring up the 5 orphaned functions (charts 38-42) adds 5 new charts to the dashboard. If they have bugs (e.g., `plot_caller_tier_heatmap` has a known `break` bug), they could produce broken charts. → **Mitigation**: Fix the `break` bug in `plot_caller_tier_heatmap` before wiring, or remove the orphans entirely. The proposal leaves this decision open — either wire+fix or delete.
