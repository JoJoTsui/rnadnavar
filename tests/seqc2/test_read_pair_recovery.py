"""Original RNA recovery CLI: exact identifiers and synchronized mate outputs."""
import json

import pytest
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def fastq(path, names):
    path.write_text("".join(f"@{name}\nACGT\n+\nIIII\n" for name in names))


def recover(tmp_path, names, mate1, mate2, raw_r1=None):
    ids = tmp_path / "read_ids.txt"
    ids.write_text("\n".join(names) + "\n")
    r1, r2 = tmp_path / "input_R1.fastq", tmp_path / "input_R2.fastq"
    fastq(r1, mate1)
    fastq(r2, mate2)
    if raw_r1 is not None:
        r1.write_text(raw_r1)
    return subprocess.run([
        sys.executable, str(ROOT / "bin/recover_rna_read_pairs.py"),
        "--read-ids", str(ids), "--fastq-1", str(r1), "--fastq-2", str(r2),
        "--library", "RT_library_A", "--out-prefix", str(tmp_path / "selected"),
        "--stats", str(tmp_path / "stats.json"),
    ], capture_output=True, text=True)


def headers(path):
    return path.read_text().splitlines()[::4]


def test_preserves_numeric_read_names_and_separates_singletons(tmp_path):
    result = recover(tmp_path, ["read11", "read22", "orphan", "absent"],
                     ["read11/1", "read22/1", "orphan/1"], ["read22/2", "read11/2"])
    assert result.returncode == 0, result.stderr
    assert headers(tmp_path / "selected_R1.fastq") == ["@read11/1", "@read22/1"]
    assert headers(tmp_path / "selected_R2.fastq") == ["@read11/2", "@read22/2"]
    assert headers(tmp_path / "selected_singleton_R1.fastq") == ["@orphan/1"]
    assert headers(tmp_path / "selected_singleton_R2.fastq") == []
    stats = json.loads((tmp_path / "stats.json").read_text())
    assert stats["paired"] == 2
    assert stats["singleton_r1"] == 1
    assert stats["missing"] == 1
    assert stats["library"] == "RT_library_A"


def test_duplicate_selected_names_fail_without_outputs(tmp_path):
    result = recover(tmp_path, ["read11"], ["read11/1", "read11/1"], ["read11/2"])
    assert result.returncode != 0
    assert "duplicate" in result.stderr.lower()
    assert "RT_library_A" in result.stderr
    assert not (tmp_path / "selected_R1.fastq").exists()
    assert not (tmp_path / "stats.json").exists()


@pytest.mark.parametrize("raw", [
    "@read11/1\nACGT\n+\n",
    "@read11/1\nACGT\n+\nIII\n",
    "read11/1\nACGT\n+\nIIII\n",
    "@read11/1\nACGT\ninvalid\nIIII\n",
])
def test_malformed_fastq_fails_before_publishing(tmp_path, raw):
    result = recover(tmp_path, ["read11"], [], ["read11/2"], raw_r1=raw)
    assert result.returncode != 0
    assert "FASTQ" in result.stderr
    assert not (tmp_path / "selected_R1.fastq").exists()
    assert not (tmp_path / "stats.json").exists()
