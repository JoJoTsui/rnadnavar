# bar-text-optimization

## Purpose

Unified text mark rendering for all bar charts with auto-contrast color, responsive fontSize, overlap avoidance strategy, and consistent SI-prefix formatting. Removes dead `_add_bar_labels` code.

## ADDED Requirements

### Requirement: Auto-contrast text color
The system SHALL compute text color based on bar fill luminance in `_make_bar_text`. When a `color_col` or fill encoding is specified, the system SHALL determine whether the bar fill is dark (luminance < 0.5) or light and set text color to white or black accordingly. When no fill is specified, text SHALL default to `"#333"` (dark gray, visible on white/default background).

#### Scenario: White text on dark bars
- **WHEN** a bar segment uses dark fill color (#1f77b4, blue)
- **THEN** text is rendered in white

#### Scenario: Black text on light bars
- **WHEN** a bar segment uses light fill color (#ff7f0e, orange)
- **THEN** text is rendered in black

### Requirement: Responsive font sizing
The system SHALL scale `fontSize` in `_make_bar_text` based on the number of x-axis categories. Default mapping: ≤10 categories → fontSize=10, 11-25 → fontSize=9, 26-50 → fontSize=8, >50 → fontSize=7. The `fontSize` parameter SHALL still be overridable by callers.

#### Scenario: Few bars get larger text
- **WHEN** a bar chart has 8 x-axis categories
- **THEN** fontSize defaults to 10

#### Scenario: Many bars get smaller text
- **WHEN** a bar chart has 60 x-axis categories
- **THEN** fontSize defaults to 7

### Requirement: Overlap avoidance strategy
The system SHALL accept an `overlap` parameter in `_make_bar_text` with values: `"hide"` (skip labels on bars narrower than estimated text width, default), `"rotate"` (rotate labels 90° when bars would overlap), `"stagger"` (alternate dy=-8 and dy=-16 for adjacent bars). When `overlap="hide"` and the bar chart is faceted or has many categories (>30), the system SHALL skip rendering text labels to prevent visual clutter.

#### Scenario: Hide on dense chart
- **WHEN** a bar chart has 60+ categories and overlap="hide"
- **THEN** no text labels are rendered

#### Scenario: Rotate on narrow bars
- **WHEN** overlap="rotate" and bars are estimated to be <40px wide
- **THEN** text labels are rendered with angle=90, align="left", dy=-4

### Requirement: Unified SI-prefix formatting
All text labels in `_make_bar_text` SHALL use the same formatting path: count values SHALL be formatted via an internal `_si()` helper (≥1M → "1.5M", ≥1K → "1.5K", else integer). When `show_pct=True`, the label SHALL be `"{count} ({pct}%)"` using the SI-formatted count. When `show_pct=False`, the label SHALL be the SI-formatted count only.

#### Scenario: SI formatting for large counts
- **WHEN** a bar segment has count=1500000 and show_pct=True with 75%
- **THEN** label is "1.5M (75%)"

#### Scenario: SI formatting for small counts
- **WHEN** a bar segment has count=42 and show_pct=False
- **THEN** label is "42"

### Requirement: Dead code removal
The system SHALL remove `_add_bar_labels()` function (lines 255-264 in visualizer.py) which has zero callers in the codebase. Any references to this function in imports or exports SHALL also be removed.

#### Scenario: Dead code removed
- **WHEN** visualizer.py is inspected
- **THEN** _add_bar_labels function does not exist
- **AND** all bar charts continue to use _make_bar_text
