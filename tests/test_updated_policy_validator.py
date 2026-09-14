"""Unit tests for the read-only validation report contract."""
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "examples/seqc2/scripts/validate_updated_policy.py"
spec = importlib.util.spec_from_file_location("policy_validator", SCRIPT)
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


def test_metric_tuple_uses_record_counts_only():
    metrics = {"records": {"tp": 1, "fp": 2, "fn": 3, "precision": 0.33}}
    assert validator.metric_tuple(metrics) == (1, 2, 3)


def test_expected_matrix_has_all_four_target_cells():
    assert set(validator.EXPECTED) == {
        "wes_ll/ukb", "wes_ll/medexome", "wgs_il/ukb", "wgs_il/medexome"
    }


def test_check_can_mark_inconclusive_without_calling_it_pass():
    result = validator.check("unknown", False, "insufficient evidence", severity="inconclusive")
    assert result["status"] == "inconclusive"
