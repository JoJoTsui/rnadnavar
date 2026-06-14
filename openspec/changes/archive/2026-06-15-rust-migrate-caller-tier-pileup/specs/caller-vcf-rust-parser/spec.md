## ADDED Requirements

### Requirement: Parse normalized caller VCFs with Rust
The system SHALL parse 6 caller VCFs (DNA/RNA × Mutect2/DeepSomatic/Strelka) using Rust noodles-vcf from bcftools-normalized files at `normalized/<caller>/*.dec.norm.vcf.gz` and `vcf_realignment/normalized/<caller>/*.dec.norm.vcf.gz`. Parsing SHALL release the GIL via `py.detach()`.

#### Scenario: DNA Mutect2 FORMAT extraction
- **WHEN** a DNA Mutect2 normalized VCF is parsed with target positions from the rescue VCF
- **THEN** FORMAT fields GT, AD (REF+ALT), AF, DP are extracted for the tumor sample (suffix DT) at matching (CHROM, POS, REF, ALT) positions

#### Scenario: DNA Strelka FORMAT extraction
- **WHEN** a DNA Strelka normalized VCF is parsed
- **THEN** FORMAT fields DP, TAR, TIR, TOR, AU, CU, GU, TU (tier1) are extracted for the TUMOR sample, and GT/AD fields are NOT extracted (Strelka lacks them)

#### Scenario: Sample suffix matching
- **WHEN** a caller VCF contains multiple samples
- **THEN** the tumor sample is identified by suffix matching: DT for DNA tumor, RT for RNA tumor, TUMOR for Strelka

### Requirement: Match caller variants on all 4 coordinate columns
The system SHALL match caller VCF records to rescue VCF variants on (CHROM, POS, REF, ALT), not just (CHROM, POS). Target positions SHALL be a set of 4-tuples.

#### Scenario: Multiallelic site with different REF/ALT
- **WHEN** a caller VCF has two records at the same position with different ALT alleles
- **THEN** each record is matched to the rescue variant with the corresponding REF/ALT, not conflated by position alone

#### Scenario: Missing caller variant
- **WHEN** a target position exists in the rescue VCF but not in a caller VCF
- **THEN** all FORMAT fields for that caller are null-filled at that position

### Requirement: Early termination with target position set
The system SHALL stop scanning a caller VCF when all target positions have been found.

#### Scenario: All targets found mid-scan
- **WHEN** scanning a caller VCF and the set of remaining target positions becomes empty
- **THEN** the scan terminates immediately, returning results for all found positions

### Requirement: Rust caller parser releases the GIL
The system SHALL release the Python GIL during caller VCF parsing so multiple samples' caller VCFs can be parsed concurrently.

#### Scenario: Parallel caller parsing
- **WHEN** two threads call `parse_caller_vcf()` on different VCF files simultaneously
- **THEN** total wall time is approximately max(thread1, thread2), not sum(thread1, thread2)
