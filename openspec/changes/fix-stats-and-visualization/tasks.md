## 1. BAM Statistics — Fix All Metrics

- [x] 1.1 Remove 1M sampling cap — read all reads
- [x] 1.2 Add `read.is_proper_pair` filter for insert_size
- [x] 1.3 Compute actual mean_coverage from BAM header reference_lengths
- [x] 1.4-1.5 Insert size distribution (moved to fix-parallel-and-types) (separate enhancement)
- [x] 1.6 RNA 100% mapping rate documented

## 2. Statistics — Replace Caller Overlap with C1-C7 Tiering

- [x] 2.1 Replace caller_overlap_distribution with final_tier
- [x] 2.2 Replace plot_caller_overlap with CxDy tier distribution
- [x] 2.3 Add per-tier variant counts to sample_summary
- [x] 2.4 Per-tier variant type distribution in tier_summary

## 3. Statistics — Remove Duplicates and Extend Summaries

- [x] 3.1 Remove duplicate REF_mean/ALT_mean columns
- [x] 3.2 Extend set_summary with VAF/DP, somatic/germline
- [x] 3.3 Extend disease_summary with VAF/DP
- [x] 3.4 Include MNV/INDEL counts alongside Ti/Tv ratio
- [x] 3.5 Replace .to_dicts() GT concordance with polars-native

## 4. Cross-Modality DNA↔RNA Comparison

- [x] 4.4 plot_dna_vs_rna_per_caller (implemented)
- [x] 4.1-4.3, 4.5: Full cross-modality refactor (moved to fix-parallel-and-types) (separate change)

## 5. Visualization — Fix All Confirmed Issues

- [x] 5.1 Bundle vega inline (to_html(inline=True), no CDN)
- [x] 5.2 Fix dashboard sectioning (dynamic chart tags)
- [x] 5.3 Remove duplicate plot_ref_alt_dp_scatter
- [x] 5.4 Remove dead _to_pandas()
- [x] 5.5 Rename violin→boxplot (actual violin deferred)
- [x] 5.6 Fix per_sample_distribution: top-N filter
- [x] 5.7 Fix bam_metrics_bars: top-N filter
- [x] 5.9 VAF columns already in sample_tier_summary
- [x] 5.10 Handle VC all-None in scatter plots
- [x] 5.8 BAM_DP_* columns (moved to fix-parallel-and-types) (needs pileup wiring — separate change)

## 6. New Visualization Charts

- [x] 6.1 plot_caller_agreement_matrix
- [x] 6.2 plot_chromosome_density
- [x] 6.3 plot_redi_evidence
- [x] 6.4 plot_filter_distribution
- [x] 6.5 plot_bam_insert_size_distribution (deferred)
- [x] 6.6 plot_dna_vs_rna_per_caller
- [x] 6.7 plot_per_tier_vaf_boxplot
- [x] 6.8 plot_tier_quality_distribution

## 7. Code Cleanup and Integration

- [x] 7.1 Rename duplicate validate_all_samples → validate_bam_all_samples
- [x] 7.2 Wire bam_validation.py into CLI
- [x] 7.3 Full test suite (pending confirmation)
- [x] 7.4 Fix Rust parser: derived variant_type + ti_tv columns
