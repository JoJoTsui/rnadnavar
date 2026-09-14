"""Tests for descriptive pileup handling; no alignment is run by the tests."""
import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "examples/seqc2/scripts/audit_rescue_read_evidence.py"
spec = importlib.util.spec_from_file_location("read_evidence", SCRIPT)
read_evidence = importlib.util.module_from_spec(spec)
spec.loader.exec_module(read_evidence)


def test_failed_pileup_is_inconclusive():
    result = read_evidence.pileup("/bin/false", Path("ref.fa"), Path("reads.bam"), "chr1", 10)
    assert result["status"] == "inconclusive"


def test_empty_pileup_is_inconclusive(tmp_path):
    tool = tmp_path / "empty-samtools"
    tool.write_text("#!/bin/sh\nexit 0\n")
    tool.chmod(0o755)
    result = read_evidence.pileup(str(tool), Path("ref.fa"), Path("reads.bam"), "chr1", 10)
    assert result["status"] == "inconclusive"
