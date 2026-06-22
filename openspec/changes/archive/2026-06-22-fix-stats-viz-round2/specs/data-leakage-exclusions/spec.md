## ADDED Requirements

### Requirement: CLI flags for data leakage exclusions
The system SHALL support `--exclude-disease`, `--min-vaf`, and `--min-dp` CLI flags. When specified, filtered parquet files SHALL be written to a separate `variant_details_filtered/` directory containing only variants that pass the exclusion criteria. The full dataset in `variant_details/` SHALL remain unchanged. CSV and chart outputs SHALL be generated from the full dataset only. Filtered files are intended for downstream ML tasks (zero-shot disease experiments, low-VAF validation, low-DP studies).

#### Scenario: Exclude specific diseases for zero-shot experiments
- **WHEN** `--exclude-disease BLCA,BRCA` is specified
- **THEN** variants from BLCA and BRCA disease samples SHALL be excluded from `variant_details_filtered/`
- **AND** all other diseases' variants SHALL be present in the filtered parquet

#### Scenario: Filter by minimum VAF
- **WHEN** `--min-vaf 0.05` is specified
- **THEN** variants with DNA_VAF_mean < 0.05 AND RNA_VAF_mean < 0.05 SHALL be excluded from the filtered parquet
- **AND** variants passing the threshold in either modality SHALL be retained

#### Scenario: Filter by minimum DP
- **WHEN** `--min-dp 10` is specified
- **THEN** variants with DNA_DP_mean < 10 AND RNA_DP_mean < 10 SHALL be excluded from the filtered parquet

#### Scenario: No exclusion flags
- **WHEN** none of `--exclude-disease`, `--min-vaf`, or `--min-dp` are specified
- **THEN** no `variant_details_filtered/` directory SHALL be created
