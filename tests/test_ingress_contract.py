"""Tests for explicit ingress assignment parsing."""
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "examples/seqc2/scripts/validate_ingress_contract.py"
spec = importlib.util.spec_from_file_location("ingress", SCRIPT)
ingress = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ingress)


def test_assignment_parses_groovy_and_yaml_like_values():
    text = "step = \'variant_calling\'\ntools = \'consensus,rescue\'\n"
    assert ingress.assignment(text, "step") == "variant_calling"
    assert ingress.assignment(text, "tools") == "consensus,rescue"


def test_missing_assignment_is_explicit():
    assert ingress.assignment("# step omitted\n", "step") is None
