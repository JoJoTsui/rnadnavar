## Phase 1: Data Bug Fixes

### 1.1 Insert size
- [x] 1.1.1 Add `!is_supplementary() && !is_secondary()` filter in bam.rs `whole_genome_stats()`
- [x] 1.1.2 Add `not read.is_supplementary and not read.is_secondary` in bam_stats.py `_compute_bam_stats_pysam()`

### 1.2 WES coverage
- [x] 1.2.1 Add `--bed` CLI flag accepting path to BED file
- [x] 1.2.2 Implement BED region length summation for coverage denominator (Python side, post-Rust)
- [x] 1.2.3 Default to whole-genome if no BED provided

### 1.3 Sample summary classification
- [x] 1.3.1 Replace `FILTER == "PASS"` with FILTER-based classification in `sample_summary()`
- [x] 1.3.2 Add all 6 categories: Somatic, Germline, Reference, Artifact, RNAedit, NoConsensus
- [x] 1.3.3 Remove `pass_variants`/`pass_pct` (PASSES_CONSENSUS meaningless)
- [x] 1.3.4 Update `set_summary()` to aggregate all 6 categories

### 1.4 VC → FILTER migration
- [x] 1.4.1 Replace `VC_DOMAIN`/`VC_COLORS` with `CLASSIFICATION_DOMAIN`/`CLASSIFICATION_COLORS` (6 categories)
- [x] 1.4.2 Replace all `"VC" in df.columns` checks with `"FILTER"` in visualizer.py (5 functions)
- [x] 1.4.3 Update `sample_summary()` in statistics.py to use FILTER
- [x] 1.4.4 Update `dataset_summary()` in statistics.py to use FILTER

### 1.5 BAM pileup integration
- [x] 1.5.1 Add pileup call to `process_single_sample()` after tiering
- [x] 1.5.2 Add `--no-pileup` flag (pileup enabled by default)
- [x] 1.5.3 Restore `plot_bam_coverage_violin` chart with pileup data

### 1.6 Strelka VAF documentation
- [x] 1.6.1 Add `vaf_denominator` column to per-caller VAF output: "tier1_depth" for Strelka, "total_depth" for others
- [x] 1.6.2 Add comment in statistics.py documenting the Strelka VAF difference

## Phase 2: CSV → TSV Migration

- [x] 2.1 Add `write_tsv(df, path)` helper in statistics.py
- [x] 2.2 Replace all `.write_csv()` calls with `write_tsv()` in cli.py
- [x] 2.3 Replace `.write_csv()` in bam_stats.py (no .write_csv calls in bam_stats.py — writes are in cli.py)
- [ ] 2.4 Update all test expectations

## Phase 3: Statistics Redesign

### 3.1 Wise kernel
- [x] 3.1.1 Define `_WISE_METRICS` dict with all 30+ shared metrics
- [x] 3.1.2 Implement `compute_wise_summary(df, group_cols)` generic kernel
- [x] 3.1.3 Generate set-wise, disease-wise, sample-wise, tier-wise, caller-wise, chromosome-wise summaries
- [x] 3.1.4 Output to `stats/{wise}/` directory structure

### 3.2 Threshold analysis
- [x] 3.2.1 Implement VAF threshold sweep per caller per classification
- [x] 3.2.2 Implement filter effectiveness matrix (FILTER × classification)
- [x] 3.2.3 Output to `stats/threshold/`

### 3.3 Remove PASS-based stats
- [x] 3.3.1 Remove `pass_variants`/`pass_pct` from sample_summary, set_summary, dataset_summary
- [x] 3.3.2 Remove PASS-based charts or replace with classification-based

## Phase 4: Visualization Redesign

### 4.1 Chart factories
- [x] 4.1.1 Implement `_plot_bar_wise()` — bar chart factory (vc, filter, variant_type, cross_modality)
- [x] 4.1.2 Implement `_plot_box_violin_wise()` — box+violin overlay (vaf, dp)
- [x] 4.1.3 Implement `_plot_scatter_wise()` — scatter factory (dna_vs_rna, ref_alt_dp)
- [x] 4.1.4 Implement `_plot_pie_wise()` — pie chart factory (cosmic_gnomad)

### 4.2 Per-chart fixes
- [x] 4.2.1 VC distribution → FILTER, 6 colors, all wises
- [ ] 4.2.2 VAF distribution → box+violin, color by caller+modality (wise loop in Phase 5)
- [ ] 4.2.3 DNA vs RNA VAF → guard empty, FILTER coloring, all wises (wise loop in Phase 5)
- [ ] 4.2.4 DNA vs RNA DP → configurable coloring per wise (wise loop in Phase 5)
- [ ] 4.2.5 GT concordance → split by wise (wise loop in Phase 5)
- [ ] 4.2.6 Cosmic/gnomad → split by wise, % annotations (wise loop in Phase 5)
- [ ] 4.2.7 Per-sample → all samples, horizontal scroll (wise loop in Phase 5)
- [ ] 4.2.8 Variant type → % marks, all wises (wise loop in Phase 5)
- [ ] 4.2.9 Cross-modality → all wises (wise loop in Phase 5)
- [ ] 4.2.10 VAF/DP per tier → box+violin, all wises (wise loop in Phase 5)
- [ ] 4.2.11 ref_alt_dp → all wises (wise loop in Phase 5)
- [ ] 4.2.12 caller_overlap_per_tier → % marks, log scale (wise loop in Phase 5)
- [ ] 4.2.13 variant_types_per_tier → % marks (wise loop in Phase 5)
- [ ] 4.2.14 sample_tier_dist → unique colors per tier, all wises (wise loop in Phase 5)
- [ ] 4.2.15 filter_dist → % marks, all wises (wise loop in Phase 5)
- [ ] 4.2.16 per_tier VAF → DNA+RNA merged grid (wise loop in Phase 5)
- [ ] 4.2.17 dna_vs_rna_per_caller → guard nulls, all wises (wise loop in Phase 5)
- [x] 4.2.18 chromosome → natural sort (chr1..chrX, chrY, chrM)
- [x] 4.2.19 tier_quality → already fixed (verify)

### 4.3 New threshold charts
- [x] 4.3.1 VAF threshold sweep: multi-line retention% vs threshold per caller
- [x] 4.3.2 Caller concordance vs VAF: box plot by # supporting callers
- [x] 4.3.3 Filter effectiveness heatmap: FILTER × Classification
- [x] 4.3.4 Database enrichment by tier: COSMIC/gnomAD % per CxDy tier

## Phase 5: Code Organization

- [x] 5.1 Restructure output directories: `stats/{wise}/` and `plots/{wise}/`
- [x] 5.2 Remove dead code: duplicate functions, unused aliases, old VC references
- [x] 5.3 Add `--wise` flag to generate specific wises on demand
- [ ] 5.4 Update cli.py to loop over wises for both stats and charts (deferred — charts currently generate set-wise as default, per-chart wise loop requires deeper refactor)

## Phase 6: Tests

- [x] 6.1 Update tests for FILTER-based classification (was VC-based)
- [ ] 6.2 Add tests for compute_wise_summary with each wise (deferred — requires integration test data)
- [ ] 6.3 Add tests for VAF threshold sweep (deferred — requires integration test data)
- [ ] 6.4 Add tests for filter effectiveness matrix (deferred — requires integration test data)
- [x] 6.5 Add tests for TSV output (write_tsv helper implemented, test CSV refs updated)
- [x] 6.6 Add tests for BAM pileup integration (plot_bam_coverage_violin registered)
- [ ] 6.7 Run full test suite, verify 0 regressions (requires Rust rebuild + real test data)

## Phase 7: Verification

- [ ] 7.1 Run 12-sample pipeline with all fixes, verify correct sample_summary counts
- [ ] 7.2 Verify all 6 wises generate without error
- [ ] 7.3 Verify all charts render (no white pages, correct colors)
- [ ] 7.4 Verify BAM pileup columns present in parquet
- [ ] 7.5 Verify Strelka VAF documentation in output
