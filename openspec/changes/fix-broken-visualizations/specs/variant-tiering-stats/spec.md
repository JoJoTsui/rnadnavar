## ADDED Requirements

### Requirement: FINAL_TIER_ORDER is exported for downstream consumers
The `tiering_stats` module SHALL export the `FINAL_TIER_ORDER` constant (already defined at module level) for import by `visualizer.py` and other downstream consumers. The constant SHALL contain the full ordered list of tier identifiers from C1D1 through C7D0.

#### Scenario: Visualizer imports tier order
- **WHEN** visualizer.py imports from tiering_stats
- **THEN** FINAL_TIER_ORDER is available as a list of 14 tier strings
- **AND** the order is C1D1, C1D0, C2D1, C2D0, C3D1, C3D0, C4D1, C4D0, C5D1, C5D0, C6D1, C6D0, C7D1, C7D0
