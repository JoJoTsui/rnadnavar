## Why

The full-pipeline run (64 samples, June 2026) revealed 15 specific issues across statistics computation, BAM validation, and visualization rendering. Three are outright bugs (BAM validation column swap, empty coverage violin, missing BAM integrity check), several are missing features discovered during data review (depth threshold sweep, DP distribution plots, data leakage exclusions), and the rest are visualization quality issues (chromosome ordering gaps, VAF plot styling, grid-by-set layout). These fixes make the pipeline output suitable for downstream ML tasks (zero-shot disease experiments, low-VAF/DP validation studies).

## What Changes

### Bug Fixes (P0)
- **Fix `has_bam_data` column naming**: `bam_validation.py` checks `BAM_DP_{bt}` but actual columns are `BAM_{bt}_DP` — swapped prefix causes `has_bam_data=False` for all samples
- **Fix `bam_coverage_violin` empty**: relax the `len(dp_cols) < 2` gate to handle partial BAM type availability; add null-check after `transform_density`
- **Add BAM integrity check before processing**: validate BGZF EOF marker for all BAMs upfront, reject truncated files before hours of processing (raises clear error)

### Missing Features (P1)
- **Depth threshold sweep**: replicate the VAF threshold sweep pattern for total DP, REF DP, and ALT DP — compute retention % at each threshold for each caller
- **DP distribution plots**: per-caller total DP, BAM pileup total/REF/ALT DP distributions (violin + box), per-tier DP boxplots — mirror the VAF distribution chart set
- **Show all samples in BAM read counts**: replace top-20 filter with faceted grid by set_number, or 2D heatmap when sample count exceeds chart height

### Visualization Quality (P1-P2)
- **Chromosome ordering**: apply `_sort_chromosomes()` to all remaining charts that group by CHROM (plot_vc_distribution, plot_dna_vs_rna_per_caller, statistics aggregations)
- **VAF distribution styling**: consistent color palette across all VAF charts, clamp VAF axis to [0, 1] for visualization (Strelka VAF > 1 annotated separately), vertical reference lines at 0.005/0.01 thresholds
- **Strelka VAF > 1 handling**: add `vaf_clamped` column for visualization; document that raw Strelka VAF uses tier-1 depth denominator (correct by design)
- **Grid sample-wise charts by set_number**: facet `plot_bam_metrics_bars` and `plot_per_sample_tier_distribution` by set for cross-set comparison with shared axis scales
- **CxD0 COSMIC annotation**: skip D=0 tiers in `plot_database_enrichment_by_tier` or add annotation explaining D=0 means "no database evidence" by design

### Data Leakage Exclusions (P3)
- **Exclusion filters for downstream ML**: add `--exclude-disease`, `--min-vaf`, `--min-dp` flags that write filtered parquet subsets without modifying the full dataset

### Documentation (P3)
- Document rescue VCF validation methodology (what "rescued" means, how mismatch is calculated)
- Document CxDy tiering system (C1-C7 caller tiers, D0-D1 database tiers, how they combine)
- Document how BAM DP columns map to caller VAF/DP fields in validation

## Capabilities

### New Capabilities
- `bam-integrity-check`: Pre-processing BGZF EOF validation that rejects truncated BAM files before processing begins
- `depth-threshold-sweep`: DP-based threshold sweep (total/REF/ALT) computing variant retention at each threshold, mirroring the existing VAF sweep
- `data-leakage-exclusions`: CLI flags (`--exclude-disease`, `--min-vaf`, `--min-dp`) to produce filtered parquet subsets for downstream zero-shot/low-VAF/low-DP experiments

### Modified Capabilities
- `variant-visualization`: Chromosome ordering applied consistently; VAF/DP distribution plot styling unified; coverage violin plot fixed; charts faceted by set_number; database enrichment chart excludes D=0 tiers
- `bam-statistics`: BAM integrity check before stats/pileup; per-sample read counts chart shows all samples (grid by set)
- `bam-pileup-rust`: BAM validation column naming fix (BAM_DP_{bt} → BAM_{bt}_DP)

## Impact

- `bin/vcf_stats/seq2neo/bam_validation.py` — column name fix (1 line)
- `bin/vcf_stats/seq2neo/bam_stats.py` — BAM integrity check (~30 lines)
- `bin/vcf_stats/seq2neo/visualizer.py` — coverage violin fix, chromosome ordering, VAF styling, DP plots, grid-by-set, database enrichment fix (~150 lines)
- `bin/vcf_stats/seq2neo/statistics.py` — depth threshold sweep, DP distribution aggregations (~80 lines)
- `bin/vcf_stats/seq2neo/cli.py` — data leakage flags, BAM integrity call, DP sweep wiring (~40 lines)
- `docs/` — rescue validation, CxDy tiering, BAM-column mapping docs (~3 files)
