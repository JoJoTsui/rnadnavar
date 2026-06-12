## ADDED Requirements

### Requirement: Rescue VCF parsing via noodles-vcf
The system SHALL parse rescue VCF files using Rust noodles-vcf library, extracting CHROM, POS, REF, ALT, FILTER, and all INFO fields (60+ fields) into a polars DataFrame via Arrow interop. Parsing SHALL be at least 5× faster than the current cyvcf2-based implementation.

#### Scenario: Rescue VCF parsed to DataFrame
- **WHEN** a rescue VCF path is provided
- **THEN** a polars DataFrame is returned with CHROM, POS, REF, ALT, FILTER columns and all INFO fields as typed columns (Int64, Float64, Utf8, Boolean)

#### Scenario: Large rescue VCF handling
- **WHEN** a rescue VCF has 1M+ variants
- **THEN** parsing completes within the memory budget (streaming reads, no full-file buffering in Rust)

### Requirement: Caller VCF FORMAT extraction at target positions
The system SHALL extract FORMAT fields (GT, AD, DP, AF/VAF, TAR/TIR/TOR for Strelka, SB, FAD for Mutect2, AU/CU/GU/TU for Strelka) from 6 caller VCFs at positions specified by the rescue VCF. Extraction SHALL be caller-aware — Strelka callers lack GT/AD and use TAR/TOR as AD proxies.

#### Scenario: Position-targeted extraction
- **WHEN** a set of (chrom, pos) target positions is provided from the rescue VCF
- **THEN** each caller VCF is scanned once, returning FORMAT fields only for matching positions

#### Scenario: Strelka caller handling
- **WHEN** extracting from a Strelka VCF with no GT or AD in FORMAT
- **THEN** TAR[0] is returned as AD_ALT, TOR[0] as AD_REF, and GT column is null-filled

### Requirement: Python fallback when Rust unavailable
The system SHALL provide a Python fallback using cyvcf2 when the Rust _core module cannot be imported, with identical output schema. A warning SHALL be emitted indicating the fallback is in use.

#### Scenario: Rust unavailable, fallback activated
- **WHEN** `import _core` fails (no Rust build, wrong platform)
- **THEN** cyvcf2-based parsing is used with a printed warning, producing the same DataFrame schema
