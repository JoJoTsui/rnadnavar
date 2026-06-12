## 1. Add shared helpers

- [ ] 1.1 Add `_count_rows(df)` — efficient row count on LazyFrame via `.select(pl.len()).collect().item()`
- [ ] 1.2 Add `_sample_if_large(df, max_rows=5000)` — count, sample if large, collect, return eager DataFrame
- [ ] 1.3 Remove `_maybe_collect` — all charts use explicit `_sample_if_large` or `.collect()` instead

## 2. Fix sampled charts (height-before-sampling → _sample_if_large)

- [ ] 2.1 `plot_vaf_distribution` — replace `melted.height` + manual sample + `_maybe_collect` with `_sample_if_large(melted, 50000)`
- [ ] 2.2 `plot_vaf_boxplot_per_tier` — same fix
- [ ] 2.3 `plot_dp_boxplot_per_tier` — same fix
- [ ] 2.4 `plot_dna_vs_rna_vaf` — `_sample_if_large(pdf, 5000)`
- [ ] 2.5 `plot_dna_vs_rna_dp` — `_sample_if_large(pdf, 5000)`
- [ ] 2.6 `plot_dna_vs_rna_per_caller` — `_sample_if_large(pdf, 5000)` for each pair
- [ ] 2.7 `plot_ref_alt_dp_scatter` — `_sample_if_large(pdf, 5000)` for each subchart
- [ ] 2.8 `plot_per_tier_vaf_boxplot` — `_sample_if_large(pdf, 50000)`
- [ ] 2.9 `plot_bam_coverage_violin` — `_sample_if_large(melted, 50000)` (or remove — see task 5.1)

## 3. Fix count-based charts (height → _count_rows)

- [ ] 3.1 `plot_cosmic_gnomad_annotation` — replace `df.height` and `df.filter(...).height` with `_count_rows()`
- [ ] 3.2 `plot_caller_agreement_matrix` — replace `df.height` and `df.filter(...).height` with `_count_rows()`

## 4. Fix iter_rows charts (iter_rows → collect + iter_rows)

- [ ] 4.1 `plot_gt_concordance` — add `.collect()` before `.iter_rows()`
- [ ] 4.2 `plot_gt_concordance_per_tier` — add `.collect()` before subscript access and `.iter_rows()`

## 5. Remove broken / dead code

- [ ] 5.1 Remove `plot_bam_coverage_violin` from cli.py's chart list (`BAM_DP_*` columns don't exist in parquet files)
- [ ] 5.2 Remove its import from cli.py
- [ ] 5.3 Remove or comment out `plot_bam_coverage_violin` function in visualizer.py

## 6. Reorganize visualizer.py structure

- [ ] 6.1 Reorder functions into sections: Helpers → Aggregate charts → Sampled charts → Count-based charts → Small-data charts → Dashboard
- [ ] 6.2 Add section header comments for each group
- [ ] 6.3 Remove `plot_per_sample_violin` alias (backward compat, unused)

## 7. cli.py updates

- [ ] 7.1 Remove `plot_bam_coverage_violin(combined_df, ...)` from the figs.append list
- [ ] 7.2 Remove its import

## 8. Tests and verification

- [ ] 8.1 Run existing test suite, verify 0 regressions
- [ ] 8.2 Verify all 22 charts generate without crashing on lazy input
- [ ] 8.3 Verify no chart exceeds 5 GB peak RSS during generation
