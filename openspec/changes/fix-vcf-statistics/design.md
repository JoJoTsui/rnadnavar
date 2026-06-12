## Context

The seq2neo variant statistics pipeline (`bin/vcf_stats/seq2neo/`) currently aggregates per-variant data from rescue VCFs and 6 caller VCFs across 65 complete samples. Two bugs exist in the visualization layer, and the statistics lack integration with the existing CxDy variant tiering system in `bin/vcf_stats/tiering.py`. There are 14 known issues to fix.

**Current architecture:**
- `rescue_parser.py` → parses rescue VCF for 60+ INFO fields
- `caller_parser.py` → extracts DP, AD, GT, VAF from 6 individual caller VCFs
- `statistics.py` → computes per-variant VAF, DNA/RNA means, sample/set/disease aggregates
- `visualizer.py` → 12 altair charts, dashboard, PNG/SVG export
- `rescue_validator.py` → cross-validates rescue INFO means vs caller ground truth
- `cli.py` → CLI entry point, sequential sample processing with ThreadPoolExecutor

**Tiering system** (`bin/vcf_stats/tiering.py` + `tiering_engine.py`):
- CxDy hybrid tiering: C1-C7 (caller support based on concordant DNA/RNA count) × D0-D1 (database evidence)
- Uses `FILTERS_NORMALIZED` / `FILTER_NORMALIZED_*` INFO fields from rescue VCF
- `TieringEngine` computes tiers via `count_concordant_callers()` and database checks
- Tier examples: C1D1 (≥2 DNA + ≥2 RNA with DB), C7D0 (no callers, no DB)

## Goals / Non-Goals

**Goals:**
1. Fix dashboard.html rendering (only one plot shown) by extracting body fragments instead of full HTML pages
2. Fix VC color bar legend not appearing
3. Integrate CxDy tiering into variant statistics — tier each variant, then produce per-tier statistics
4. Replace VAF boxplots with violin plots, split by tier
5. Add REF_DP and ALT_DP statistics alongside existing DP/VAF
6. Add per-sample variant distribution (sample ID on x-axis)
7. Add percentage annotations to cross-modality chart
8. Statistics at three levels: per-sample, per-tier, whole-dataset
9. Per-sample per-modality BAM statistics
10. Per-variant BAM stats appended to variant_details.parquet
11. GT concordance broken down by tier

**Non-Goals:**
- maturin + Rust rewrite (defer until profiling shows polars is the bottleneck)
- Modifying the tiering engine itself (use it as-is)
- Changing the VCF file structure or pipeline output
- Real-time streaming statistics (batch processing only)

## Decisions

### D1: Dashboard fix — extract body content from to_html()
**Choice**: Parse each chart's `to_html()` output to extract only the `<body>` inner content, then wrap in a single HTML page.
**Alternative**: Use `chart.to_dict()` and serialize manually with vega-embed — more control but more code.
**Rationale**: `to_html()` returns full HTML documents. Concatenating them creates invalid HTML with multiple `<html>`/`<body>` tags. Browsers only render the first one. Extracting body fragments is the minimal fix.

### D2: Tiering integration — add tiering_stats.py as bridge module
**Choice**: Create `tiering_stats.py` that wraps the existing `tiering.py` and `tiering_engine.py` to compute tiers for parsed variants and merge tiers into the statistics pipeline.
**Alternative**: Modify rescue_parser to compute tiers during VCF parsing. Rejected because tiering depends on FILTERS_NORMALIZED fields that are already parsed by rescue_parser.
**Rationale**: The existing tiering module expects cyvcf2 VCF records directly. We'll adapt it to also accept rows from our polars DataFrame (parsed FILTERS_NORMALIZED info) to avoid re-parsing VCFs.

### D3: VIolin plots — use altair transform_density
**Choice**: Use altair's built-in `mark_area(interpolate='step')` with computed density for violin-like visualization per tier.
**Alternative**: Pre-compute density in polars. Rejected because altair's transforms are more maintainable.
**Rationale**: altair doesn't have a native `mark_violin()`. The approach is to use layered area charts mirroring each other around a central axis, organized by tier.

### D4: REF_DP / ALT_DP naming
**Choice**: REF_DP is already partially available as `{caller}_AD_REF`. Add explicit column naming `{caller}_DP` for total depth and keep `{caller}_AD_REF` / `{caller}_AD_ALT`. Add `DNA_REF_mean`, `RNA_REF_mean` alongside existing `DNA_ALT_mean`, `RNA_ALT_mean`.
**Alternative**: New column `{caller}_REF_DP` separate from AD_REF. Rejected — AD_REF already represents REF depth.

### D5: BAM statistics — deferred to BAM file parsing
**Choice**: Per-sample per-modality BAM statistics require parsing BAM files. Use pysam (already available or pip-installable) for BAM read.
**Rationale**: BAM stats (read counts, mapping quality, coverage) are orthogonal to VCF stats and need BAM file access. This is a new module.

### D6: Per-variant BAM stats — use existing caller FORMAT fields
**Choice**: Some per-variant BAM-level data already exists in caller FORMAT fields (strand bias SB, fragment depth FAD, etc.). Extract these and add to variant_details.parquet. Additional true BAM-pileup stats per-variant would require samtools mpileup which is too slow for all samples.
**Rationale**: Reuse existing data where possible. Full per-variant BAM pileup across all 65 samples × ~50K variants each is computationally impractical.

## Risks / Trade-offs

- **Tiering computation cost**: Computing tiers for millions of variants redundant with what the rescue VCF already does (INFO/FILTERS_NORMALIZED). Mitigation: compute tiers once during parsing and cache in a column.
- **Violin plot scalability**: Density computation on millions of rows is expensive. Mitigation: sample to 50K rows per plot (same pattern already used for boxplots).
- **BAM parsing memory**: Reading BAM files for 65 samples could be memory-intensive. Mitigation: stream BAM reads, compute stats incrementally.
- **API compatibility**: tiering.py uses pandas while seq2neo uses polars. Mitigation: bridge layer converts where needed, or tiering module's TieringEngine objects can be invoked with raw dicts (no pandas dependency).
