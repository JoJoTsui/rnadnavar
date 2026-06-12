## ADDED Requirements

### Requirement: Three BAM types per sample
The system SHALL discover and process three BAM files per sample: DNA normal (DN, suffix DN), DNA tumor (DT, suffix DT), and RNA tumor (RT, suffix RT). BAM files SHALL be located in `preprocessing/mapped/{prefix}{suffix}/` subdirectories.

#### Scenario: All three BAM types found
- **WHEN** a sample has DN, DT, and RT BAM files
- **THEN** all three are processed for whole-genome stats and variant-wise pileup

#### Scenario: Missing BAM type
- **WHEN** a BAM type directory does not exist
- **THEN** that BAM type's metrics are null-filled and a warning is logged

### Requirement: Whole-genome BAM statistics
The system SHALL compute per-BAM whole-genome statistics: total reads, mapped reads, mapping rate, mean coverage, mean insert size, and mean MAPQ. Statistics SHALL be computed by sampling up to 1M reads for performance.

#### Scenario: Whole-genome BAM stats
- **WHEN** a BAM file is processed for whole-genome stats
- **THEN** a row in `bam_stats.csv` is generated with read counts, mapping rate, coverage, insert size, and MAPQ

### Requirement: Variant-wise BAM pileup
The system SHALL perform per-position BAM pileup at each variant position (CHROM, POS, REF, ALT) using Rust noodles-bam. For each position, the system SHALL compute: total DP at position, REF-supporting read depth, ALT-supporting read depth, strand bias per allele (F1R2_ref, F2R1_ref, F1R2_alt, F2R1_alt), mean base quality at position, and mean mapping quality at position.

#### Scenario: Single position pileup
- **WHEN** a BAM file is queried at chr1:633987 with REF=A, ALT=G
- **THEN** DP, REF_depth, ALT_depth, strand counts per allele, mean BQ, and mean MQ are returned

#### Scenario: Position with zero coverage
- **WHEN** a position has no reads in the BAM
- **THEN** all pileup fields are null-filled for that position

### Requirement: Two pileup modes
The system SHALL support two pileup modes: `all` (process all variants, default) and `filtered` (exclude variants with FILTER=NoConsensus). The filtered mode SHALL reduce processing volume by approximately 76% while retaining all actionable variants.

#### Scenario: All-variants mode
- **WHEN** `--pileup-mode all` is specified
- **THEN** pileup is performed for every variant in the rescue VCF

#### Scenario: Filtered mode
- **WHEN** `--pileup-mode filtered` is specified
- **THEN** variants with FILTER=NoConsensus are excluded from pileup processing
