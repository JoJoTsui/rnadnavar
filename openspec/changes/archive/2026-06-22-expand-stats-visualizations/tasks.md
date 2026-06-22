## 1. Bar Text Mark Optimization

- [x] 1.1 Rewrite `_make_bar_text()` in `visualizer.py` — add auto-contrast text color via luminance calculation, responsive fontSize based on x-axis category count, `overlap` parameter ("hide"/"rotate"/"stagger"), unified SI-prefix formatting via internal `_si()` helper for ALL text modes (not just pct mode)
- [x] 1.2 Remove `_add_bar_labels()` function from `visualizer.py` — delete lines 255-264, remove from any imports or exports
- [x] 1.3 Update all 21 bar chart call sites to use optimized `_make_bar_text` with appropriate overlap strategy per chart type (hide for sample-wise >50 categories, default for set/disease/tier wise)
- [x] 1.4 Tests: `test_make_bar_text_auto_contrast_dark`, `test_make_bar_text_auto_contrast_light`, `test_responsive_font_size`, `test_overlap_hide_dense`, `test_overlap_rotate`, `test_si_formatting`, `test_add_bar_labels_removed`

## 2. Multi-Allelic Visualization

- [x] 2.1 Implement `plot_multiallelic_classification()` — grouped bar chart of multiallelic_class counts per set_number, with count labels
- [x] 2.2 Implement `plot_allele_balance_scatter()` — scatter of major vs minor VAF at multi-allelic sites, colored by multiallelic_class, with y=x reference line
- [x] 2.3 Implement `plot_vaf_sum_histogram()` — histogram of vaf_sum with reference line at 1.0, distinct color for >1.0 bins
- [x] 2.4 Implement `plot_category_conflict_summary()` — horizontal bar chart of top 10 conflict types, ordered by count
- [x] 2.5 Register all 4 multi-allelic charts in `cli.py` global chart list (after rescue analytics section), output to `plots/multi_allelic/`
- [x] 2.6 Tests: `test_multiallelic_classification_chart`, `test_allele_balance_scatter`, `test_vaf_sum_histogram`, `test_category_conflict_summary`

## 3. Per-Caller VAF+DP Scatters

- [x] 3.1 Extend `plot_dna_vs_rna_per_caller()` to accept `metrics` parameter (tuple, default `("VAF", "DP")`), loop over metrics for each caller pair
- [x] 3.2 Implement 2-row × 3-column layout: VAF row top, DP row bottom, shared color via `.resolve_scale(color="shared")`
- [x] 3.3 DP scatters: clip to [0, 2000] with axis limits, consistent with `plot_dna_vs_rna_dp`
- [x] 3.4 Register in disease-wise chart list in `cli.py` (set and caller wises already registered in Proposal 1)
- [x] 3.5 Tests: `test_per_caller_vaf_dp_dual`, `test_per_caller_layout`, `test_per_caller_dp_clip`

## 4. BAM Alignment Stats Expansion

- [x] 4.1 Add `_compute_duplication_rate()` to `bam_stats.py` — count reads with 0x400 flag / total reads in pysam fallback
- [x] 4.2 Add `_compute_properly_paired_pct()` — count properly paired mapped reads / total mapped reads
- [x] 4.3 Add `_compute_coverage_bins()` — compute fraction of BED bases at 1×, 10×, 20×, 50×, 100× thresholds (requires BED)
- [x] 4.4 Add `_compute_insert_size_stddev()` — stddev of insert sizes for properly paired reads with positive insert size
- [x] 4.5 Update `compute_bam_stats()` and `_compute_bam_stats_pysam()` to call new metric functions, add columns to output DataFrame
- [x] 4.6 Tests: `test_duplication_rate`, `test_properly_paired_pct`, `test_coverage_bins`, `test_insert_size_stddev` — covered by rust-bam-expanded-metrics (P3) tests

## 5. BAM Visualization

- [x] 5.1 Implement `plot_bam_metrics_sample_wise()` — faceted bar chart of BAM metrics per sample, colored by bam_type, ordered by set_number
- [x] 5.2 Implement `plot_bam_coverage_distribution()` — grouped bar chart of mean coverage bin percentages per bam_type + set_number
- [x] 5.3 Register both BAM charts in `cli.py` global chart list, output to `plots/bam/`
- [x] 5.4 Tests: `test_bam_metrics_sample_wise`, `test_bam_coverage_distribution`

## 6. Integration & Verification

- [x] 6.1 Run existing 140 core tests — verify no regressions from text mark changes
- [x] 6.2 Run pipeline on a single sample — verify new multi-allelic charts generate successfully
- [x] 6.3 Run pipeline on a single sample — verify per-caller VAF+DP chart generates with 6 panels
- [x] 6.4 Run pipeline — verify BAM metrics charts generate with expanded stats
- [x] 6.5 Verify text labels are readable on all bar charts (no overlap on dense charts, visible contrast on stacked bars)
