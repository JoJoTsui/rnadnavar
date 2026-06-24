## ADDED Requirements

### Requirement: Rust tier engine returns C7 for zero concordant callers
The Rust tier engine SHALL return tier C7 (caller_tier="C7") when the concordant DNA and RNA caller counts are both zero, WITHOUT falling back to raw `N_DNA_CALLERS_SUPPORT` / `N_RNA_CALLERS_SUPPORT`. This matches the Python `TieringEngine` behavior and correctly assigns RNAedit variants (which have 0 concordant callers by design) to the C7 tier.

#### Scenario: RNAedit variant with 0 concordant callers
- **WHEN** a variant has FILTER=RNAedit, 0 concordant DNA callers, 0 concordant RNA callers, but 3 raw RNA callers
- **THEN** the Rust tier engine returns caller_tier="C7" (NOT "C3")
- **AND** the variant is assigned to C7D1 (with database support) or C7D0 (without)

#### Scenario: Unclassified variant with 0 concordant callers
- **WHEN** a variant has 0 concordant DNA callers and 0 concordant RNA callers and no database support
- **THEN** the Rust tier engine returns caller_tier="C7", database_tier="D0"
- **AND** the variant is assigned to C7D0
