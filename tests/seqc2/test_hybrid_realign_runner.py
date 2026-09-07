"""Public-contract tests for the isolated SEQC2 hybrid realignment runner."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


pytest.importorskip("yaml", reason="SEQC2 runner requires PyYAML")


ROOT = Path(__file__).resolve().parents[2]


def load_runner():
    spec = importlib.util.spec_from_file_location(
        "seqc2_runner", ROOT / "examples/seqc2/scripts/run_pipeline.py"
    )
    runner = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(runner)
    return runner


def test_completion_validator_rejects_stale_artifacts(tmp_path):
    runner = load_runner()
    stale = tmp_path / "stale.vcf.gz"
    stale.write_bytes(b"not a bgzip VCF")
    validator = tmp_path / "reject_stale.py"
    validator.write_text("import sys; raise SystemExit(1)\n")
    cfg = {"completion_validator": ["python3", str(validator), "{outdir}"]}

    assert not runner.is_complete(tmp_path, ["stale.vcf.gz"], cfg)


def test_completion_validator_is_part_of_success_contract(tmp_path):
    runner = load_runner()
    artifact = tmp_path / "final.vcf.gz"
    artifact.write_bytes(b"artifact")
    validator = tmp_path / "accept.py"
    validator.write_text("import sys; raise SystemExit(0)\n")
    cfg = {"completion_validator": ["python3", str(validator), "{outdir}"]}

    assert runner.is_complete(tmp_path, ["final.vcf.gz"], cfg)


def test_hybrid_realign_config_requires_the_second_round_artifact():
    runner = load_runner()
    config = runner.load_config(
        ROOT / "examples/seqc2/hybrid/config_realign.yaml", {}
    )

    assert Path(config["rdv_conf"]).name == "seqc2.hybrid.realign.config"
    assert any("vcf_realignment/rescue" in artifact
               for artifact in config["completion_artifacts"])
    assert config["completion_validator"]
