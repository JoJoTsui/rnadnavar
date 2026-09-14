## ADDED Requirements

### Requirement: Dashboard uses shared CDN-loaded vega-embed runtime
The system SHALL generate `dashboard.html` with the vega-embed, vega-lite, and vega JS libraries loaded once via `<script src>` CDN tags in the `<head>` section. Each chart SHALL be embedded as a `<div>` container with a `<script>vegaEmbed(...)` call containing only the vega-lite JSON spec. The system SHALL NOT inline the full vega-embed JS library per chart.

#### Scenario: Dashboard with 106 charts
- **WHEN** 106 charts are generated and passed to generate_dashboard()
- **THEN** dashboard.html contains exactly one vega-embed `<script src>` tag in `<head>`
- **AND** each chart div contains only a vegaEmbed call with JSON spec
- **AND** the total file size is under 2MB

#### Scenario: Dashboard body content extraction
- **WHEN** individual chart HTML is extracted via _extract_body_content()
- **THEN** the function SHALL skip `<script>` blocks before searching for the `<body>` tag
- **AND** no raw JavaScript library fragments SHALL appear as visible text in chart divs

### Requirement: Dashboard charts are deterministic
The system SHALL use a fixed seed for all sampling operations in chart generation. The seed SHALL be a module-level constant (default: 42). Two pipeline runs over identical input SHALL produce byte-identical chart specs (up to vega-embed version differences).

#### Scenario: Reproducible sampling
- **WHEN** _sample_if_large is called on a DataFrame with more than max_rows
- **THEN** df.sample(max_rows, seed=SAMPLE_SEED) is called
- **AND** two consecutive runs produce identical sampled data
