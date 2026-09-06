"""Focused public checks for SEQC2 hybrid-realignment validation."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "hybrid_realign_validation",
    ROOT / "examples/seqc2/hybrid/scripts/validate_realign_run.py",
)
validator = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(validator)


def row(tmp_path: Path, status: int, sample: str, library: str) -> dict[str, str]:
    fastq_1 = tmp_path / "r1.fastq.gz"
    fastq_2 = tmp_path / "r2.fastq.gz"
    fastq_1.touch(exist_ok=True)
    fastq_2.touch(exist_ok=True)
    return {
        "patient": "P",
        "status": str(status),
        "sample": sample,
        "input_stage": "raw_reads",
        "library": library,
        "fastq_1": str(fastq_1),
        "fastq_2": str(fastq_2),
        "bam": "",
        "bai": "",
    }


def test_manifest_accepts_one_logical_pooled_rt_sample(tmp_path):
    rows = [row(tmp_path, 0, "DN", "dn"), row(tmp_path, 1, "DT", "dt")]
    rows.extend(row(tmp_path, 2, "RT", library) for library in ("A", "B", "C"))

    validator.validate_manifest(rows)


def test_manifest_rejects_cross_sample_pooled_rt(tmp_path):
    rows = [row(tmp_path, 0, "DN", "dn"), row(tmp_path, 1, "DT", "dt"), row(tmp_path, 2, "RT_A", "A"), row(tmp_path, 2, "RT_B", "B")]

    with pytest.raises(ValueError, match="one logical sample"):
        validator.validate_manifest(rows)


def test_manifest_rejects_duplicate_rna_library(tmp_path):
    rows = [row(tmp_path, 0, "DN", "dn"), row(tmp_path, 1, "DT", "dt"), row(tmp_path, 2, "RT", "A"), row(tmp_path, 2, "RT", "A")]

    with pytest.raises(ValueError, match="unique and non-empty"):
        validator.validate_manifest(rows)
