## ADDED Requirements

### Requirement: Pre-processing BAM integrity check
The system SHALL check every BAM file's BGZF EOF marker before computing statistics or pileup. If ALL BAM files for a sample are truncated, the sample SHALL be skipped with an error. If some BAMs are valid and some truncated, the valid ones SHALL be processed and truncated ones logged as warnings.

#### Scenario: Sample with all truncated BAMs
- **WHEN** a sample's DN, DT, and RT BAM files all lack a valid BGZF EOF marker
- **THEN** the sample SHALL be skipped before any statistics computation
- **AND** an error message SHALL list each truncated file path

### Requirement: Per-sample BAM read counts show all samples
The BAM metrics bar chart SHALL show all samples, faceted by set_number when multiple sets are present. No sample SHALL be excluded by a top-N filter.

#### Scenario: 65 samples across 4 sets
- **WHEN** 65 samples from 4 sets have BAM statistics
- **THEN** the BAM metrics chart SHALL render all 65 samples in a faceted grid
- **AND** y-axis scales SHALL be consistent across facets
