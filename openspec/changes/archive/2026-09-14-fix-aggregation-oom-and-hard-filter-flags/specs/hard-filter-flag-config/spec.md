# hard-filter-flag-config

## Purpose

Hard filter conditions are defined as a Python config module stored alongside the code, providing a single source of truth for variant quality flags that are observed (not dropped) during the statistics stage.

## Requirements

### Requirement: Hard filter conditions defined in Python config module
The system SHALL store hard filter conditions in `bin/vcf_stats/seq2neo/hard_filter_config.py` as a Python module following the project's existing config pattern (matching `bin/common/tier_config.py`). Each condition SHALL include: name, description, severity level, flag column name, and a polars expression builder function.

#### Scenario: Config module is importable
- **WHEN** `from vcf_stats.seq2neo.hard_filter_config import HARD_FILTER_CONDITIONS` is executed
- **THEN** `HARD_FILTER_CONDITIONS` SHALL be a list of dicts, each containing keys `name`, `description`, `severity`, `flag_column`, `required_columns`, and `build_expr`

#### Scenario: Conditions have valid expression builders
- **WHEN** `build_expr(columns)` is called for any condition with available columns
- **THEN** it SHALL return a valid `pl.Expr` that evaluates to a boolean value
- **AND** it SHALL return `None` when required columns are missing from the input column set

#### Scenario: All eight hard filter conditions are present
- **WHEN** the config module is loaded
- **THEN** `HARD_FILTER_CONDITIONS` SHALL contain exactly 8 conditions: `no_caller_support`, `vaf_overflow`, `noise_allele`, `no_coverage`, `no_alt_evidence`, `germline_low_vaf`, `somatic_loh`, `reference_with_signal`

### Requirement: Flag column naming convention
Each hard filter condition SHALL generate a boolean flag column named `flag_hard_<name>` where `<name>` matches the condition's `name` field. Two summary columns SHALL also be generated: `hard_filter_flags` (comma-joined string of active flag names) and `n_hard_flags` (integer count of active flags).

#### Scenario: Flag columns named correctly
- **WHEN** the `no_caller_support` condition triggers for a variant
- **THEN** `flag_hard_no_support` SHALL be True
- **AND** the `hard_filter_flags` column SHALL contain "no_caller_support" among its values
- **AND** `n_hard_flags` SHALL be at least 1

#### Scenario: No flags triggered
- **WHEN** a variant matches none of the hard filter conditions
- **THEN** all `flag_hard_*` columns SHALL be False
- **AND** `hard_filter_flags` SHALL be an empty string
- **AND** `n_hard_flags` SHALL be 0

### Requirement: Severity levels
Each condition SHALL be assigned one of three severity levels: `high`, `medium`, or `low`, stored in the `severity` field.

#### Scenario: Severity assignment
- **WHEN** conditions `no_caller_support`, `vaf_overflow`, `noise_allele` are defined
- **THEN** their severity SHALL be `high`
- **WHEN** conditions `no_coverage`, `no_alt_evidence`, `germline_low_vaf`, `somatic_loh` are defined
- **THEN** their severity SHALL be `medium`
- **WHEN** condition `reference_with_signal` is defined
- **THEN** its severity SHALL be `low`

### Requirement: Conditions use only existing columns
Each condition's `required_columns` SHALL list the exact columns it references, and `build_expr` SHALL check column availability before building the expression. If any required column is missing from the parquet schema, the condition SHALL be silently skipped.

#### Scenario: Missing column gracefully skipped
- **WHEN** `build_hard_filter_flag_exprs(columns)` is called with a column set that lacks `flag_vaf_overflow`
- **THEN** the `vaf_overflow` condition SHALL be omitted from the returned expressions
- **AND** no error SHALL be raised
- **AND** other conditions with available columns SHALL still be included
