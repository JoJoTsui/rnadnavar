# bam-alignment-visualization

## Purpose

Sample-wise faceted BAM metrics and coverage distribution charts from expanded BAM statistics. Consumes `bam_stats.tsv` produced by `bam_stats.py`.

## ADDED Requirements

### Requirement: Four new BAM metrics extracted
The system SHALL compute 4 additional BAM metrics in `compute_bam_stats()`: `duplication_rate_pct` (fraction of duplicate reads, computed from flagged duplicate count / total reads), `coverage_bins` (fraction of reference bases covered at 1×, 10×, 20×, 50×, 100× thresholds, available when BED is provided), `properly_paired_pct` (fraction of reads with both mates mapped in proper pair), and `insert_size_stddev` (standard deviation of insert sizes for properly paired reads).

#### Scenario: Metrics extracted
- **WHEN** a BAM file is processed for statistics
- **THEN** bam_stats.tsv contains columns: duplication_rate_pct, cov_1x_pct, cov_10x_pct, cov_20x_pct, cov_50x_pct, cov_100x_pct, properly_paired_pct, insert_size_stddev

### Requirement: Sample-wise BAM metrics chart
The system SHALL generate `plot_bam_metrics_sample_wise`, a faceted bar chart showing BAM metrics per sample, grouped by `bam_type` (DN/DT/RT). Each metric SHALL be a separate faceted panel. Samples SHALL be ordered by set_number then sample_id. The chart SHALL use `_make_bar_text` for count/percentage labels.

#### Scenario: Sample-wise BAM metrics
- **WHEN** bam_stats.tsv contains metrics for 60+ samples × 3 BAM types
- **THEN** a faceted bar chart is generated with one panel per metric, bars colored by bam_type

### Requirement: Coverage distribution chart
The system SHALL generate `plot_bam_coverage_distribution`, a grouped bar chart showing the mean coverage bin percentages across samples, grouped by bam_type and set_number. Coverage bins SHALL be: cov_1x_pct, cov_10x_pct, cov_20x_pct, cov_50x_pct, cov_100x_pct.

#### Scenario: Coverage distribution
- **WHEN** bam_stats.tsv contains coverage bin columns
- **THEN** a grouped bar chart shows mean coverage bin percentages with error bars per bam_type
