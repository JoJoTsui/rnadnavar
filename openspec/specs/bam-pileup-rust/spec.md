# bam-pileup-rust

## Purpose

Rust-based per-position BAM pileup at variant positions using noodles-bam indexed queries. Releases the GIL for parallel execution across multiple BAM files.

## Requirements

### Requirement: Per-position BAM pileup with Rust noodles-bam
The system SHALL perform per-position BAM pileup at variant positions using Rust noodles-bam indexed queries (`reader.query()` with CSI index). The pileup SHALL release the GIL via `py.detach()`.

#### Scenario: Pileup at a variant position
- **WHEN** pileup is performed at a (CHROM, POS, REF, ALT) position with a valid BAI/CSI index
- **THEN** the system returns DP (total depth), REF_DP (reference allele depth), ALT_DP (alternate allele depth), F1R2_ref, F2R1_ref, F1R2_alt, F2R1_alt (strand counts), mean_BQ (mean base quality), and mean_MQ (mean mapping quality) for that position

#### Scenario: Position with no reads
- **WHEN** pileup is performed at a position where no reads align
- **THEN** all metrics are null-filled for that position

#### Scenario: Missing BAM index
- **WHEN** the BAM index file (.bai or .csi) does not exist
- **THEN** pileup falls back to pysam-based pileup with a warning

#### Scenario: BAM file not found
- **WHEN** the BAM file does not exist at the expected path
- **THEN** all positions return null-filled metrics and a warning is logged

### Requirement: Rust BAM pileup releases the GIL
The system SHALL release the Python GIL during BAM pileup so multiple samples' pileups can run concurrently.

#### Scenario: Parallel pileup
- **WHEN** two threads call `pileup_variants()` on different BAM files simultaneously
- **THEN** total wall time is approximately max(thread1, thread2), not sum(thread1, thread2)

### Requirement: Pileup output matches pysam reference
Rust pileup output SHALL match pysam pileup output for the same positions within a tolerance of 1 count for depth metrics and 0.1 for quality means.

#### Scenario: Output parity with pysam
- **WHEN** 100 random variant positions are pileup'd by both Rust and pysam
- **THEN** DP, REF_DP, and ALT_DP match exactly, and mean_BQ/mean_MQ differ by no more than 0.1
