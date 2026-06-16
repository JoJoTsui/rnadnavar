## ADDED Requirements

### Requirement: Pre-processing BAM integrity validation
The system SHALL validate BAM file integrity before any variant processing or statistics computation. Validation SHALL check the BGZF EOF marker (magic bytes `1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00`) in the last 28 bytes of each BAM file. If ALL BAM files for a sample are truncated (missing valid EOF), the sample SHALL be skipped with an error message listing the affected files. If some BAMs are valid and some truncated, the valid BAMs SHALL be processed and the truncated ones SHALL be logged as warnings.

#### Scenario: All BAMs valid
- **WHEN** all BAM files for a sample contain a valid BGZF EOF marker
- **THEN** the integrity check passes silently and processing continues

#### Scenario: All BAMs truncated
- **WHEN** ALL BAM files for a sample (DN, DT, RT) are missing the BGZF EOF marker
- **THEN** the sample SHALL be skipped with a clear error message listing each truncated file path
- **AND** no variant processing or pileup SHALL be attempted for that sample

#### Scenario: Some BAMs truncated
- **WHEN** one or two BAM files are truncated but at least one is valid
- **THEN** a warning SHALL be printed for each truncated BAM
- **AND** processing SHALL continue with the valid BAMs only
