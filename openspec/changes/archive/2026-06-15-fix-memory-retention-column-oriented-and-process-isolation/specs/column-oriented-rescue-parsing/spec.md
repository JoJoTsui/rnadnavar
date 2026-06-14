## ADDED Requirements

### Requirement: Rust rescue parser returns column-oriented data
The Rust VCF rescue parser SHALL return variant data as a dictionary mapping column names
to lists of values (`{column_name: [values]}`), instead of a list of per-record dictionaries.
Each column list SHALL have the same length, with `None` representing missing values.

#### Scenario: All INFO fields present
- **WHEN** a rescue VCF contains all INFO fields declared in its header
- **THEN** the returned dictionary SHALL contain one key per INFO field, one key per fixed
  column (CHROM, POS, REF, ALT, FILTER), and one key per derived column (variant_type, ti_tv)
- **AND** all value lists SHALL have identical length equal to the number of parsed variants

#### Scenario: Missing INFO field in a record
- **WHEN** a variant record lacks a particular INFO field declared in the header
- **THEN** the corresponding column SHALL contain `None` at that row position

#### Scenario: Empty VCF
- **WHEN** a rescue VCF contains no variant records
- **THEN** the parser SHALL return an empty dictionary

### Requirement: Column types are preserved from Rust to polars
The Rust parser SHALL return values with types matching the VCF INFO field declarations:
Integer fields as Python `int`, Float fields as Python `float`, Flag fields as Python `bool`,
and String/Character fields as Python `str`. Type casting in Python (`_cast_columns`) SHALL
NOT be required.

#### Scenario: Integer INFO field
- **WHEN** an INFO field is declared as `Type=Integer` in the VCF header
- **THEN** the corresponding Python list SHALL contain `int` values (not strings)

#### Scenario: Float INFO field
- **WHEN** an INFO field is declared as `Type=Float` in the VCF header
- **THEN** the corresponding Python list SHALL contain `float` values (not strings)

#### Scenario: Flag INFO field
- **WHEN** an INFO field is declared as `Type=Flag` in the VCF header and is present in a record
- **THEN** the corresponding Python list SHALL contain `True` at that position
- **AND** when absent, the value SHALL be `False`

### Requirement: Derived columns computed during parsing
The Rust parser SHALL compute `variant_type` (SNV/INS/DEL/MNV) and `ti_tv`
(transition/transversion) columns during the VCF parse pass, without requiring a separate
Python iteration step.

#### Scenario: SNV classification
- **WHEN** REF and ALT are both single nucleotides
- **THEN** variant_type SHALL be "SNV"
- **AND** ti_tv SHALL be `True` for transitions (A↔G, C↔T) and `False` for transversions

#### Scenario: Insertion classification
- **WHEN** ALT is longer than REF
- **THEN** variant_type SHALL be "INS"
- **AND** ti_tv SHALL be `None`

#### Scenario: Deletion classification
- **WHEN** REF is longer than ALT
- **THEN** variant_type SHALL be "DEL"
- **AND** ti_tv SHALL be `None`

### Requirement: Column-oriented output matches row-oriented output
The column-oriented parser SHALL produce the same polars DataFrame (same values, same column
names, same dtypes) as the existing row-oriented parser when given the same VCF input, after
accounting for type-correct output.

#### Scenario: Parity with row-oriented parser
- **WHEN** the same rescue VCF is parsed by both `parse_rescue` (row-oriented) and
  `parse_rescue_columns` (column-oriented)
- **THEN** the resulting DataFrames SHALL have identical column sets
- **AND** all values SHALL be identical for each column and row

### Requirement: Python fallback uses polars-native derived columns
When the Rust parser is unavailable (cyvcf2 fallback path), the Python code SHALL use polars
vectorized `when/then/otherwise` expressions for `variant_type` and `ti_tv` computation instead
of `.to_list()` + Python list comprehensions.

#### Scenario: polars-native variant_type in Python fallback
- **WHEN** the cyvcf2 parser is used (Rust unavailable)
- **THEN** variant_type SHALL be computed using `pl.when().then().otherwise()` chain
- **AND** `.to_list()` SHALL NOT be called on REF or ALT columns

#### Scenario: polars-native ti_tv in Python fallback
- **WHEN** the cyvcf2 parser is used (Rust unavailable)
- **THEN** ti_tv SHALL be computed using vectorized string operations
- **AND** `.to_list()` SHALL NOT be called on REF or ALT columns
