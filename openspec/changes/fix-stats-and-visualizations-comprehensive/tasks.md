## Phase 1: Data Bug Fixes

### 1.1 Insert size
- [ ] 1.1.1 Add `!is_supplementary() && !is_secondary()` filter in bam.rs `whole_genome_stats()`
- [ ] 1.1.2 Add `not read.is_supplementary and not read.is_secondary` in bam_stats.py `_compute_bam_stats_pysam()`

### 1.2 WES coverage
- [ ] 1.2.1 Add `--bed` CLI flag accepting path to BED file
- [ ] 1.2.2 Implement BED region length summation for coverage denominator (Python side, post-Rust)
- [ ] 1.2.3 Default to whole-genome if no BED provided

### 1.3 Sample summary classification
- [ ] 1.3.1 Replace `FILTER == "PASS"` with FILTER-based classification in `sample_summary()`
- [ ] 1.3.2 Add all 6 categories: Somatic, Germline, Reference, Artifact, RNAedit, NoConsensus
- [ ] 1.3.3 Remove `pass_variants`/`pass_pct` (PASSES_CONSENSUS meaningless)
- [ ] 1.3.4 Update `set_summary()` to aggregate all 6 categories

### 1.4 VC → FILTER migration
- [ ] 1.4.1 Replace `VC_DOMAIN`/`VC_COLORS` with `CLASSIFICATION_DOMAIN`/`CLASSIFICATION_COLORS` (6 categories)
- [ ] 1.4.2 Replace all `"VC" in df.columns` checks with `"FILTER"` in visualizer.py (5 functions)
- [ ] 1.4.3 Update `sample_summary()` in statistics.py to use FILTER
- [ ] 1.4.4 Update `dataset_summary()` in statistics.py to use FILTER

### 1.5 BAM pileup integration
- [ ] 1.5.1 Add pileup call to `process_single_sample()` after tiering
- [ ] 1.5.2 Add `--no-pileup` flag (pileup enabled by default)
- [ ] 1.5.3 Restore `plot_bam_coverage_violin` chart with pileup data

### 1.6 Strelka VAF documentation
- [ ] 1.6.1 Add `vaf_denominator` column to per-caller VAF output: "tier1_depth" for Strelka, "total_depth" for others
- [ ] 1.6.2 Add comment in statistics.py documenting the Strelka VAF difference

## Phase 2: CSV → TSV Migration

- [ ] 2.1 Add `write_tsv(df, path)` helper in statistics.py
- [ ] 2.2 Replace all `.write_csv()` calls with `write_tsv()` in cli.py
- [ ] 2.3 Replace `.write_csv()` in bam_stats.py
- [ ] 2.4 Update all test expectations

## Phase 3: Statistics Redesign

### 3.1 Wise kernel
- [ ] 3.1.1 Define `_WISE_METRICS` dict with all 30+ shared metrics
- [ ] 3.1.2 Implement `compute_wise_summary(df, group_cols)` generic kernel
- [ ] 3.1.3 Generate set-wise, disease-wise, sample-wise, tier-wise, caller-wise, chromosome-wise summaries
- [ ] 3.1.4 Output to `stats/{wise}/` directory structure

### 3.2 Threshold analysis
- [ ] 3.2.1 Implement VAF threshold sweep per caller per classification
- [ ] 3.2.2 Implement filter effectiveness matrix (FILTER × classification)
- [ ] 3.2.3 Output to `stats/threshold/`

### 3.3 Remove PASS-based stats
- [ ] 3.3.1 Remove `pass_variants`/`pass_pct` from sample_summary, set_summary, dataset_summary
- [ ] 3.3.2 Remove PASS-based charts or replace with classification-based

## Phase 4: Visualization Redesign

### 4.1 Chart factories
- [ ] 4.1.1 Implement `_plot_bar_wise()` — bar chart factory (vc, filter, variant_type, cross_modality)
- [ ] 4.1.2 Implement `_plot_box_violin_wise()` — box+violin overlay (vaf, dp)
- [ ] 4.1.3 Implement `_plot_scatter_wise()` — scatter factory (dna_vs_rna, ref_alt_dp)
- [ ] 4.1.4 Implement `_plot_pie_wise()` — pie chart factory (cosmic_gnomad)

### 4.2 Per-chart fixes
- [ ] 4.2.1 VC distribution → FILTER, 6 colors, all wises
- [ ] 4.2.2 VAF distribution → box+violin, color by caller+modality
- [ ] 4.2.3 DNA vs RNA VAF → guard empty, FILTER coloring, all wises
- [ ] 4.2.4 DNA vs RNA DP → configurable coloring per wise
- [ ] 4.2.5 GT concordance → split by wise
- [ ] 4.2.6 Cosmic/gnomad → split by wise, % annotations
- [ ] 4.2.7 Per-sample → all samples, horizontal scroll
- [ ] 4.2.8 Variant type → % marks, all wises
- [ ] 4.2.9 Cross-modality → all wises
- [ ] 4.2.10 VAF/DP per tier → box+violin, all wises
- [ ] 4.2.11 ref_alt_dp → all wises
- [ ] 4.2.12 caller_overlap_per_tier → % marks, log scale
- [ ] 4.2.13 variant_types_per_tier → % marks
- [ ] 4.2.14 sample_tier_dist → unique colors per tier, all wises
- [ ] 4.2.15 filter_dist → % marks, all wises
- [ ] 4.2.16 per_tier VAF → DNA+RNA merged grid
- [ ] 4.2.17 dna_vs_rna_per_caller → guard nulls, all wises
- [ ] 4.2.18 chromosome → natural sort (chr1..chrX, chrY, chrM)
- [ ] 4.2.19 tier_quality → already fixed (verify)

### 4.3 New threshold charts
- [ ] 4.3.1 VAF threshold sweep: multi-line retention% vs threshold per caller
- [ ] 4.3.2 Caller concordance vs VAF: box plot by # supporting callers
- [ ] 4.3.3 Filter effectiveness heatmap: FILTER × Classification
- [ ] 4.3.4 Database enrichment by tier: COSMIC/gnomAD % per CxDy tier

## Phase 5: Code Organization

- [ ] 5.1 Restructure output directories: `stats/{wise}/` and `plots/{wise}/`
- [ ] 5.2 Remove dead code: duplicate functions, unused aliases, old VC references
- [ ] 5.3 Add `--wise` flag to generate specific wises on demand
- [ ] 5.4 Update cli.py to loop over wises for both stats and charts

## Phase 6: Tests

- [ ] 6.1 Update tests for FILTER-based classification (was VC-based)
- [ ] 6.2 Add tests for compute_wise_summary with each wise
- [ ] 6.3 Add tests for VAF threshold sweep
- [ ] 6.4 Add tests for filter effectiveness matrix
- [ ] 6.5 Add tests for TSV output
- [ ] 6.6 Add tests for BAM pileup integration
- [ ] 6.7 Run full test suite, verify 0 regressions

## Phase 7: Verification

- [ ] 7.1 Run 12-sample pipeline with all fixes, verify correct sample_summary counts
- [ ] 7.2 Verify all 6 wises generate without error
- [ ] 7.3 Verify all charts render (no white pages, correct colors)
- [ ] 7.4 Verify BAM pileup columns present in parquet
- [ ] 7.5 Verify Strelka VAF documentation in output
