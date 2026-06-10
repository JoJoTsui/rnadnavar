## 1. Fix correctness bugs

- [x] 1.1 Fix `filter_distribution` sort in statistics.py — change `descending=[False, True]` to `descending=False`
- [x] 1.2 Add warning log when rescue VCF has zero variants (silent sample drop) in cli.py
- [x] 1.3 Remove unused imports: `Path` in bam_stats.py, `numpy` in caller_parser.py (if unused)

## 2. Streaming CLI architecture (cli.py)

- [x] 2.1 Create `variant_details/` subdirectory under output_dir for per-sample parquet files
- [x] 2.2 Write each sample's DataFrame to `variant_details/{sample_id}_variants.parquet` immediately after processing
- [x] 2.3 Delete DataFrame reference (`del df`) and call `gc.collect()` after write to free memory
- [x] 2.4 Remove `samples_data` dict — no longer needed
- [x] 2.5 Replace `pl.concat(list(samples_data.values()))` with `pl.scan_parquet(str(variant_dir / "*_variants.parquet"))`
- [x] 2.6 Pre-compute aggregation DataFrames from lazy scan for visualization (vc_counts, tier_counts, etc.)
- [x] 2.7 Pass pre-aggregated DataFrames to chart functions instead of combined_df
- [x] 2.8 Sample 10K rows from lazy frame for scatter plots (VAF scatter, DP scatter, per-caller scatter)

## 3. Lazy aggregation support (statistics.py, tiering_stats.py)

- [x] 3.1 Add `isinstance(df, pl.LazyFrame)` check + `.collect()` to `dataset_summary`
- [x] 3.2 Add lazy support to `disease_summary` (group_by on LazyFrame returns LazyFrame, call .collect())
- [x] 3.3 Add lazy support to `sample_tier_summary`
- [x] 3.4 Add lazy support to `caller_overlap_distribution`, `filter_distribution`, `variant_type_distribution`
- [x] 3.5 Add lazy support to `gt_concordance` (materialize once for multiple `.filter().height` calls)
- [x] 3.6 Add lazy support to `tier_summary` in tiering_stats.py

## 4. Pre-aggregated visualization (visualizer.py)

- [x] 4.1 Remove all `.to_pandas()` calls — altair uses polars DataFrames directly
- [x] 4.2 Update bar chart functions to accept pre-grouped DataFrames (vc_distribution, caller_overlap, variant_type, filter, cross_modality, redi_evidence, tiered charts)
- [x] 4.3 Update box plot functions to accept pre-computed quartile DataFrames (vaf_boxplot, dp_boxplot, per_tier_vaf_boxplot) — compute quartiles in polars, pass to altair as pre-aggregated
- [x] 4.4 Update scatter functions to accept sampled DataFrames (dna_vs_rna_vaf, dna_vs_rna_dp, ref_alt_dp_scatter, dna_vs_rna_per_caller)
- [x] 4.5 Update `plot_gt_concordance` to accept pre-computed concordance dict
- [x] 4.6 Update `plot_gt_concordance_per_tier` for pre-computed data
- [x] 4.7 Update `plot_caller_agreement_matrix` for pre-computed data
- [x] 4.8 Update `plot_chromosome_density` for pre-grouped chromosome counts
- [x] 4.9 Update `plot_tier_quality_distribution` for pre-grouped quality distribution
- [x] 4.10 Update `plot_cosmic_gnomad_annotation` for pre-computed annotation counts
- [x] 4.11 Ensure `plot_bam_metrics_bars` and `plot_bam_coverage_violin` unchanged (already small inputs)
- [x] 4.12 Ensure `plot_per_sample_distribution` unchanged (uses sample_stats_df, already small)
- [x] 4.13 Replace pandas `.pivot()` in `plot_validation_heatmap` with polars `.pivot()`

## 5. Memory-efficient join_caller_columns (caller_parser.py)

- [x] 5.1 Replace `pl.DataFrame(col_data)` dict-of-lists with `pl.Series` per column, specifying dtype
- [x] 5.2 Build DataFrame from list of Series objects: `pl.DataFrame([pl.Series("CHROM", values, dtype=pl.Utf8), ...])`
- [x] 5.3 Free intermediate lists after Series construction

## 6. Memory regression tests

- [x] 6.1 Add `TestMemoryEfficiency` class to test_seq2neo_stats.py
- [x] 6.2 `test_dataframe_freed_after_write`: verify DataFrame deleted and gc'd after parquet write
- [x] 6.3 `test_lazy_scan_parity`: run 2 samples in streaming mode, verify aggregation output matches eager mode
- [x] 6.4 `test_dataset_summary_lazy`: verify dataset_summary works with LazyFrame input
- [x] 6.5 `test_disease_summary_lazy`: verify disease_summary works with LazyFrame
- [x] 6.6 `test_visualizer_pre_aggregated`: verify chart functions accept pre-grouped DataFrames
- [x] 6.7 `test_no_pandas_conversion`: verify chart functions don't call .to_pandas()
- [x] 6.8 `test_filter_distribution_sort`: verify sort order is correct after fix

## 7. End-to-end verification

- [ ] 7.1 Run on 4 samples, verify all CSV outputs match previous eager mode
- [ ] 7.2 Run on 20 samples, verify peak memory < 4 GB
- [ ] 7.3 Run full non-e2e test suite, ensure 0 regressions
