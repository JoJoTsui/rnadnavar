## Context

The seq2neo statistics pipeline (`bin/vcf_stats/seq2neo/`) processes ~12 samples per `--set` run, writing per-sample variant parquet files to `variant_details/` and aggregating them via `pl.scan_parquet("*")`. The pipeline produces wise-based statistics (set, disease, sample, tier, caller, chromosome, variant-category) and 22+ chart types across 7 wise directories, plus a `dashboard.html`.

The existing `fix-broken-visualizations` change fixed rendering-level bugs: dashboard body extraction, axis sorting, tier ordering, faceting toggle, and sampling seed. However, the underlying data quality and filtering problems remain unaddressed.

### Current Architecture (data flow)

```
Manifest (parquet)
  ├── set_number (may be 0/NULL)
  ├── disease (may be "")
  └── disease_normalized (may be "")
       │
       ▼
cli.py: main()
  ├── [1] Filter manifest by --set / --sample-ids
  ├── [2] Per-sample processing → variant_details/{sid}_variants.parquet
  │         └── disease/set_number written as-is from manifest
  ├── [3] combined_df = pl.scan_parquet("*")  ← SCANS ALL FILES, NO SAMPLE FILTER
  ├── [4] BAM stats → bam_stats.tsv (fresh: filtered; resume: all)
  ├── [5] sample_summary → sample_summary.tsv (fresh: filtered; resume: all)
  ├── [6] Wise summaries + charts use combined_df (unfiltered)
  └── [7] Confidence tiers computed but no wise dim / charts / ML export
```

### Key Constraints

- Per-sample parquet files must persist across `--set` runs (no cleanup after each run)
- The `--resume` flag must work with previously written parquet files
- Rust `stats_core` backend is compiled as a PyO3 extension module
- Altair/Vega-Lite is the only charting library (Plotly removed)
- Polars (eager + lazy) is the only DataFrame library (pandas only used for altair bridge)
- The pipeline runs on machines with 32-64 cores and 256GB+ RAM

## Goals / Non-Goals

**Goals:**
- Fix `combined_df` to only include the current run's sample_ids
- Fix disease/set data quality: no more empty string `""` diseases, warnings when all `set_number=0`
- Fix `--resume` path to filter `sample_summary.tsv` and `bam_stats.tsv` reloads
- Extend Rust BAM stats to compute `cov_*_pct` for whole-genome mode
- Fix chart faceting coverage (4 unfaceted functions) and disease registry asymmetry
- Verify `caller_tier` column exists in `combined_df` schema; warn/recompute if missing
- Fix BAM chart layouts (metrics_sample_wise grid, bam_metrics_bars top_n)
- Add confidence analysis pipeline: wise dimension, charts, filters, ML export

**Non-Goals:**
- Redesigning the chart set or adding new chart types beyond confidence analysis
- Changing the Altair vs Plotly split (the parent `vcf_stats/` still uses Plotly)
- Fixing multi-allelic visualization logic — that belongs in `fix-multiallelic-logic`
- Fixing BAM pileup data quality (MAPQ=255) — that belongs in `fix-bam-stats-rust`
- Re-architecting the parquet directory structure (one dir per `--set` is a separate design decision)
- Adding a new disease taxonomy or reclassifying existing disease names

## Decisions

### Decision 1: combined_df sample_id filter — filter after scan, not separate directories

**Chosen**: Add `combined_df = combined_df.filter(pl.col("sample_id").is_in(current_sample_ids))` immediately after `pl.scan_parquet("*")` at cli.py:1103. This is a single-line fix that ensures all downstream stats, charts, and validation use only the current-run's sample_ids.

**Alternatives considered**:
- Per-set parquet directories (`variant_details/set1/`, `set2/`, etc.) → Rejected: breaks `--resume` behavior, requires migration of existing parquet files, complicates glob pattern
- Delete old parquets before each run → Rejected: defeats the purpose of `--resume`, requires maintaining a parquet index

**Rationale**: The filter-at-scan approach is zero-overhead (predicate pushdown in polars parquet reader skips entire row groups from excluded samples), maintains backward compatibility with existing parquet files, and is a single line of code.

### Decision 2: Disease normalization — fix in-process, not upstream manifest

**Chosen**: Add "Unknown" sentinel in `process_single_sample` (cli.py) when `disease_normalized` is empty or whitespace-only. This is the last possible normalization point before data flows into parquet files. Also fix `parse_projects_to_json.py` to use the last seen disease or "Unknown" instead of `""`.

**Alternatives considered**:
- Fix only in `parse_projects_to_json.py` (upstream) → Rejected: doesn't handle existing manifests with empty disease, doesn't guard against future manifest sources
- Fix in `build_sample_manifest.py` → Rejected: good belt, but need suspenders too; the parquet writer is the final gate
- Map to a controlled vocabulary → Rejected: premature; "Unknown" is sufficient as a sentinel; taxonomy standardization is a separate project

**Rationale**: Belt-and-suspenders approach — normalize upstream (parse_projects_to_json.py), validate midstream (build_sample_manifest.py warning), and sanitize downstream (cli.py sentinel).

### Decision 3: WG coverage bins — streaming histogram in Rust

**Chosen**: Extend `whole_genome_stats_impl` in `bam.rs` to accumulate a depth histogram (array of 5 counters: ≥1x, ≥10x, ≥20x, ≥50x, ≥100x) during the existing whole-genome scan. Instead of per-position tracking (which would require a genome-sized array), increment the histogram during each read's CIGAR walk: for each alignment match base, increment the base position counter; when a position's depth passes a threshold, increment that threshold's counter once (not per-read). Use the BAM header's reference sequence lengths as the denominator.

**Alternatives considered**:
- BAI-indexed WG scan (like BED path) → Rejected: BAI indexes on 16Kbp bins; for WG you'd need ~200K queries, each with BAM seek overhead — prohibitively slow
- Python-side pileup (pysam) → Rejected: would duplicate the existing Rust BAM reading; much slower (Python GIL + per-position iteration)
- Use `mosdepth` subprocess → Rejected: external dependency, no Python API, requires BED file generation for target regions

**Rationale**: The whole-genome scan already walks every alignment match for mean_coverage. Adding 5 counter increments per base position is ~5ns overhead per base — negligible compared to the existing CIGAR walk. The depth histogram approach uses O(1) memory regardless of genome size.

### Decision 4: Chart faceting — add `_apply_faceting` to 4 unfaceted functions

**Chosen**: Add `_apply_faceting(chart, group_col)` call to `plot_filter_distribution`, `plot_redi_evidence`, `plot_ref_alt_dp_scatter`, and `plot_dna_vs_rna_dp`. These functions already receive `group_col` from the wise registries but never facet by it. The `_apply_faceting` function already has rules for `set_number` (2 columns), `disease_normalized` (4 columns), and `FILTER` (3 columns) — the 4 unfaceted functions will automatically inherit correct faceting.

**Alternatives considered**:
- Add `facet_col` parameter to all 4 functions → Rejected: `_apply_faceting` is already the canonical faceting mechanism; adding a separate `facet_col` parameter would create two competing faceting code paths
- Remove these charts from set/disease registries → Rejected: these charts are useful for per-set/disease comparison; removing them degrades the dashboard

**Rationale**: Consistent faceting across all chart functions. The `_apply_faceting` mechanism is already well-tested on `plot_caller_overlap` and `plot_gt_concordance`. Adding it to the remaining 4 functions aligns their behavior.

### Decision 5: Confidence analysis — new wise dimension, not a separate pipeline

**Chosen**: Add `confidence` as a new wise dimension in the existing `wise_configs` framework (`cli.py:1317`), create a `_WISE_CHART_REGISTRY["confidence"]` entry with 8-10 chart types, add `--min-confidence-tier` to `build_unified_filter()`, and add `--export-high-confidence` as a post-processing step that writes `high_confidence_variants.parquet`.

**Alternatives considered**:
- Separate CLI entry point (`run_confidence_analysis.py`) → Rejected: adds maintenance burden; the wise framework already handles grouping, metrics, and chart generation generically
- Add confidence as a FILTER value → Rejected: `confidence_tier` is an orthogonal dimension to `FILTER` (Somatic/Germline/Reference); cross-tabulating them is more informative than merging

**Rationale**: The wise framework is designed for this — `compute_wise_summary(df, ["confidence_tier"])` reuses the same 27 shared metrics. No new statistics computation code needed.

### Decision 6: --resume filtering — reload then filter, don't rescan

**Chosen**: After reloading `sample_summary.tsv` and `bam_stats.tsv` in `--resume` mode, filter their rows to only the current `sample_ids` from the manifest. For `sample_summary.tsv`, also recompute from parquet files if the TSV doesn't contain the necessary columns.

**Alternatives considered**:
- Always recompute from parquet on `--resume` → Rejected: defeats the purpose of `--resume` (speed); recomputing sample_summary for 100 samples takes 2-3 minutes
- Ignore: `--resume` users should accept all-sample data → Rejected: user explicitly asked for `--set 1`, getting all sets is a bug

**Rationale**: Filtering a loaded DataFrame is O(n) and near-instant. The filtered TSV is still written to the output directory (it represents the current run's scope), but the original unfiltered file is preserved on disk.

## Risks / Trade-offs

- **[combined_df filter may hide missing samples]**: If a user expects `--set 1` to include 12 samples but the manifest only has 11, the filter silently drops the 12th sample's parquet (even if it exists). → **Mitigation**: Log the number of samples included/excluded by the filter. If 0 samples match, error out with a clear message.
- **[WG coverage bins memory for >1B reads]**: A sample with 1B+ mapped reads could increment the depth histogram 1B+ times. This is fine (5× u64 counters = 40 bytes, trivially small). The per-base depth tracking during CIGAR walks uses `u64` accumulators that don't overflow until 1.8e19 bases — well beyond any real genome.
- **["Unknown" disease may group unrelated samples]**: If 3 different diseases are all empty-string, they will now all be "Unknown" — correctly grouped together rather than silently merged. This is an improvement over the current behavior.
- **[Confidence export may produce very large files]**: HIGH-confidence variants may still be millions of rows for large cohorts. → **Mitigation**: Export as parquet (columnar, compressed), not CSV. Log the count before writing.
- **[caller_tier verification requires schema check]**: Calling `collect_schema()` on a LazyFrame is O(1) (no data scan), but if `caller_tier` is missing, computing tiers for the full `combined_df` may take minutes (Rust fallback path is O(n) with n=millions). → **Mitigation**: Only recompute if missing; log a clear warning. In the common case (fresh run), `caller_tier` is present from per-sample processing.

## Open Questions

1. **Should `variant_details/` be cleaned between runs?** Currently parquets accumulate. The filter fix handles this, but disk usage grows unbounded. This is deferred — a separate `--clean-variant-details` flag could be added later.
2. **Should `disease_normalized` use a controlled vocabulary?** Currently it's free-text (lowercased). "Colorectal Cancer" vs "colorectal" vs "CRPC" are all distinct. A mapping table with canonical names would improve grouping but requires domain expertise to build. Deferred.
3. **Should `coverage_distribution` use a single chart or faceted by bam_type?** Currently one chart with bam_type on x-axis. Could also facet by bam_type for denser comparison. Deferred — user can request this after the fix.
