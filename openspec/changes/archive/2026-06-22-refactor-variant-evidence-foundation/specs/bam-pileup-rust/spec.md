# bam-pileup-rust (delta)

## MODIFIED Requirements

### Requirement: Per-position BAM pileup with Rust noodles-bam
The system SHALL perform per-position BAM pileup at variant positions using Rust noodles-bam indexed queries. The position list SHALL be determined by the unified filter pipeline — when filters exclude variants, their positions SHALL be excluded from the pileup position list. There is no separate `--pileup-mode` flag; pileup filtering is a consequence of the unified filter.

#### Scenario: Pileup respects unified filter
- **WHEN** unified filter excludes NoConsensus variants
- **THEN** NoConsensus positions are excluded from the pileup position list
- **AND** BAM pileup columns (BAM_DT_DP, etc.) are null for those positions

#### Scenario: Pileup with no filtering
- **WHEN** no filters are specified (all variants pass through)
- **THEN** all variant positions are included in pileup (equivalent to legacy `--pileup-mode all`)

## ADDED Requirements

### Requirement: Pileup position list from combined filter
The pileup position list construction in `process_single_sample()` SHALL use the same filter expression that applies to the combined_df. The position list SHALL be derived from the filtered DataFrame, ensuring pileup data is only computed for variants that survive the unified filter. This eliminates the need for a separate `--pileup-mode` flag.

#### Scenario: Position list constructed from filtered data
- **WHEN** `process_single_sample()` constructs pileup positions
- **THEN** positions are extracted from the DataFrame AFTER caller joining and filtering
- **AND** the resulting pileup columns are consistent with the filtered variant set

## REMOVED Requirements

None. All existing pileup quality and parity requirements remain unchanged.
