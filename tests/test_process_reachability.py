"""Regression checks for the process reachability assay."""
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "examples/seqc2/scripts/validate_process_reachability.py"
spec = importlib.util.spec_from_file_location("reachability", SCRIPT)
reachability = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reachability)


def test_script_has_explicit_static_contract():
    assert "static_only" in SCRIPT.read_text()
