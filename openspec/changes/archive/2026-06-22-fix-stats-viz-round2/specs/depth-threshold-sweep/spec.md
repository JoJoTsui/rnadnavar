## ADDED Requirements

### Requirement: Depth threshold sweep computation
The system SHALL compute a depth threshold sweep that calculates variant retention percentage at each depth threshold for each caller. The sweep SHALL cover total DP thresholds `[1, 2, 5, 10, 20, 50, 100, 200]`, REF DP thresholds `[0, 1, 2, 5, 10, 20, 50]`, and ALT DP thresholds `[0, 1, 2, 5, 10, 20, 50]`. Results SHALL be returned as a long-format DataFrame with columns: caller, metric (DP/REF_DP/ALT_DP), threshold, retention_pct.

#### Scenario: Total DP threshold sweep
- **WHEN** variant data has per-caller DP columns (DNA_mutect2_DP, DNA_strelka_DP, etc.)
- **THEN** a sweep DataFrame is generated with retention % at each DP threshold for each caller
- **AND** retention at DP=1 is near 100% and retention at DP=200 reflects high-depth variants only

#### Scenario: ALT DP threshold sweep for rescue validation
- **WHEN** variant data has ALT_DP columns from BAM pileup or caller FORMAT
- **THEN** a sweep DataFrame is generated showing ALT allele support retention at each threshold
