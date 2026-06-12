## Context

See `fix-cross-sample-aggregation-oom` and `fix-visualizer-lazy-frame-and-oom` for prior work. This change builds on those fixes and addresses remaining data quality and visualization bugs discovered during deep code review.

## Goals / Non-Goals

**Goals:** Fix 10 data bugs, reorganize stats/viz by wise dimensions, add threshold-focused analysis, migrate CSV→TSV, integrate BAM pileup.

**Non-Goals:** Changing rescue VCF semantics, modifying caller VCF parsing logic, adding new callers.

## Part 1: Data Bug Fixes

### 1.1 Insert Size — Supplementary Read Filter

**Bug**: `bam.rs:55` and `bam_stats.py:96` check `is_properly_segmented()` but not `is_supplementary()` or `is_secondary()`. Chimeric/supplementary alignments pass the proper-pair flag but have extreme TLEN values.

**Fix (Rust)**:
```rust
if flags.is_properly_segmented()
    && !flags.is_supplementary()
    && !flags.is_secondary()
{ ... }
```

**Fix (pysam)**:
```python
if read.is_proper_pair and not read.is_supplementary and not read.is_secondary
    and read.template_length and read.template_length > 0:
```

RNA insert sizes with introns (5K-15K) are correct after filtering — no cap needed.

### 1.2 Coverage — WES Denominator

**Bug**: `bam.rs:29` sums ALL reference sequence lengths (~3 GB). For WES data this under-reports coverage by ~100×.

**Fix**: Accept `--bed` flag. Read BED file, sum region lengths, use as denominator. Default to whole-genome if no BED. Requires new Rust function or Python-side coverage recalculation.

### 1.3 Sample Summary — Classification Fix

**Bug**: `sample_summary()` uses `VC` (all null) and `FILTER=="PASS"` (never matches).

**Fix**: Use `FILTER` column with correct categories:
```python
for cat in ["Somatic","Germline","Reference","Artifact","RNAedit","NoConsensus"]:
    result[f"n_{cat.lower()}"] = df.filter(pl.col("FILTER") == cat).height
```
Remove `pass_variants`/`pass_pct` (PASSES_CONSENSUS always YES — meaningless).

### 1.4 VC → FILTER Migration

**Bug**: 5 functions in visualizer.py check `"VC" in df.columns` — passes but all values null.

**Fix**: Replace "VC" with "FILTER" everywhere. Expand domain:
```python
CLASSIFICATION_DOMAIN = ["Somatic","Germline","Reference","Artifact","RNAedit","NoConsensus"]
CLASSIFICATION_COLORS = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b"]
```

### 1.5 BAM Pileup Integration

**Bug**: `pileup_variants()` exists but is never called. BAM validation `has_bam_data` always false.

**Fix**: Add to `process_single_sample`:
```python
if not no_pileup:
    pileup_df = pileup_variants(bam_path, positions)
    df = df.with_columns(pileup_df)
```
Controlled by `--no-pileup` (default: pileup enabled).

### 1.6 Strelka VAF Documentation

Strelka FORMAT/DP = tier1 filtered depth (not total depth). Strelka VAF = TAR[0]/DP_tier1 is correct but different metric from Mutect2/DeepSomatic VAF = AD/DP_total. Add `vaf_denominator` column to output.

## Part 2: Statistics Redesign

### 2.1 Wise Kernel

```python
_WISE_METRICS = {
    "n_variants": pl.len(),
    "n_somatic": (pl.col("FILTER")=="Somatic").sum(),
    "n_germline": (pl.col("FILTER")=="Germline").sum(),
    "n_reference": (pl.col("FILTER")=="Reference").sum(),
    "n_artifact": (pl.col("FILTER")=="Artifact").sum(),
    "n_rnaedit": (pl.col("FILTER")=="RNAedit").sum(),
    "n_noconsensus": (pl.col("FILTER")=="NoConsensus").sum(),
    "n_snv": (pl.col("variant_type")=="SNV").sum(),
    "n_ins": (pl.col("variant_type")=="INS").sum(),
    "n_del": (pl.col("variant_type")=="DEL").sum(),
    "n_mnv": (pl.col("variant_type")=="MNV").sum(),
    "mean_dna_vaf": pl.col("DNA_VAF_mean").mean(),
    "mean_rna_vaf": pl.col("RNA_VAF_mean").mean(),
    "mean_dna_dp": pl.col("DNA_DP_mean").mean(),
    "mean_rna_dp": pl.col("RNA_DP_mean").mean(),
    "n_cross_modality": (pl.col("CROSS_MODALITY")=="YES").sum(),
    "n_rescued": (pl.col("RESCUED")=="YES").sum(),
    "n_cosmic": pl.col("COSMIC_ID").is_not_null().sum(),
    "n_gnomad": pl.col("GNOMAD_AF").is_not_null().sum(),
}

def compute_wise_summary(df, group_cols):
    existing = {k: v for k, v in _WISE_METRICS.items()
                if all(c in df.columns for c in _extract_cols(v))}
    return df.group_by(group_cols).agg(list(existing.values()))
```

### 2.2 Threshold Analysis

VAF threshold sweep per caller:
```
threshold ∈ {0.05, 0.10, ..., 0.50}
  × caller ∈ {DNA_mutect2, DNA_deepsomatic, DNA_strelka, RNA_mutect2, RNA_deepsomatic, RNA_strelka}
  × classification ∈ {Somatic, Germline, Reference, Artifact, RNAedit, NoConsensus}
```

Filter effectiveness matrix:
```
FILTER × classification → count + pct
For each flag field (min_alt_reads, gnomad, blacklist, etc.):
  count where flag=True, broken down by classification
```

### 2.3 TSV Output

All `.write_csv()` → `write_tsv()` using `separator="\t"`. No embedded commas in VCF pipe/colon-delimited fields.

## Part 3: Visualization Redesign

### 3.1 Chart Factories

```python
def _plot_bar_wise(df, output_dir, group_col, wise_name, value_col, chart_id, title):
    """Generic bar chart for any wise. Adds percentage text marks."""

def _plot_box_violin_wise(df, output_dir, group_col, wise_name, value_cols, chart_id):
    """Box + violin overlay for VAF/DP distributions."""

def _plot_scatter_wise(df, output_dir, group_col, wise_name, x_col, y_col, color_col, chart_id):
    """Scatter plot for any wise with configurable coloring."""
```

### 3.2 Per-Chart Fixes

| Chart | Fix |
|-------|-----|
| VC distribution | FILTER column, 6 colors, all wises |
| VAF distribution | Box+violin overlay, color by caller+modality |
| DNA vs RNA VAF | Guard empty, FILTER coloring |
| DNA vs RNA DP | Coloring configurable per wise |
| GT concordance | Split by wise, facet |
| Cosmic/gnomad | Split by wise, % text annotations |
| Per-sample | All samples, horizontal scroll |
| Variant type | % marks, all wises |
| Cross-modality | All wises |
| VAF/DP per tier | Box+violin, all wises |
| ref_alt_dp | All wises |
| caller_overlap_per_tier | % marks, log scale |
| sample_tier_dist | Unique colors per tier |
| filter_dist | % marks, all wises |
| per_tier VAF | DNA+RNA merged grid |
| dna_vs_rna_per_caller | Guard nulls, all wises |
| chromosome | Natural sort (chr1..chrX) |
| tier_quality | Already fixed (_sample_if_large) |

### 3.3 New Charts

- VAF threshold sweep: multi-line retention% vs threshold, per caller
- Caller concordance vs VAF: box plot by # of supporting callers
- Filter effectiveness heatmap: FILTER × Classification
- Database enrichment by tier: COSMIC/gnomAD % per CxDy tier

## Part 4: Output Structure

```
stats/
├── set/set_summary.tsv, set_tier_summary.tsv
├── disease/disease_summary.tsv
├── sample/sample_summary.tsv
├── tier/tier_summary.tsv
├── caller/caller_summary.tsv
├── chromosome/chromosome_summary.tsv
├── threshold/vaf_sweep.tsv, filter_effectiveness.tsv
├── bam_stats.tsv
└── bam_validation.tsv
plots/
├── set/   (23 charts)
├── disease/ (VC, variant_type, filter, cosmic_gnomad, ...)
├── sample/  (per_sample, per_sample_tier)
├── tier/    (vaf, dp, gt_concordance, caller_overlap)
├── caller/  (vaf_distribution, dp_distribution)
└── chromosome/ (density)
```

## Risks / Trade-offs

- **BAM pileup is slow**: Indexed BAM queries per variant × 3 BAM types. Mitigation: `--no-pileup` flag.
- **Many charts**: 6 wises × 15 chart types = 90+ charts. Mitigation: generate on demand via `--wise` flag.
- **TSV type fidelity**: polars `write_csv(separator="\t")` preserves numeric types.
