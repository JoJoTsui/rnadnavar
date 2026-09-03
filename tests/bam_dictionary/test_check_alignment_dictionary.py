import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "bin" / "check_alignment_dictionary.py"


def write_inputs(tmp_path, alignment_sq, reference_sq, idxstats=None):
    header = tmp_path / "header.sam"
    dictionary = tmp_path / "reference.dict"
    alignment = tmp_path / "input.bam"
    stats = tmp_path / "idxstats.tsv"
    header.write_text("@HD\tVN:1.6\n" + "".join(
        f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in alignment_sq
    ))
    dictionary.write_text("@HD\tVN:1.6\n" + "".join(
        f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in reference_sq
    ))
    alignment.write_bytes(b"alignment-placeholder")
    rows = idxstats or [
        (name, length, 10, 0) for name, length in alignment_sq
    ] + [("*", 0, 0, 2)]
    stats.write_text("".join("\t".join(map(str, row)) + "\n" for row in rows))
    return header, dictionary, alignment, stats


def run_checker(tmp_path, alignment_sq, reference_sq, policy="normalize", idxstats=None):
    header, dictionary, alignment, stats = write_inputs(
        tmp_path, alignment_sq, reference_sq, idxstats
    )
    audit = tmp_path / "audit.json"
    result = subprocess.run(
        [
            sys.executable,
            str(CHECKER),
            "--header", str(header),
            "--idxstats", str(stats),
            "--reference-dict", str(dictionary),
            "--alignment", str(alignment),
            "--sample-id", "sample",
            "--stage", "input",
            "--policy", policy,
            "--output", str(audit),
        ],
        capture_output=True,
        text=True,
    )
    return result, json.loads(audit.read_text())


def test_exact_dictionary_is_compatible(tmp_path):
    sequence = [("chr1", 100), ("chr2", 80)]
    result, audit = run_checker(tmp_path, sequence, sequence)
    assert result.returncode == 0
    assert result.stdout.strip() == "false"
    assert audit["decision"] == "compatible"
    assert audit["alignment"]["sha256"]
    assert audit["reference"]["sha256"]


def test_extra_and_missing_contigs_are_safely_normalizable_and_audited(tmp_path):
    alignment = [("chr1", 100), ("virus", 40)]
    reference = [("chr1", 100), ("chr2", 80)]
    idxstats = [("chr1", 100, 10, 0), ("virus", 40, 7, 1), ("*", 0, 0, 2)]
    result, audit = run_checker(
        tmp_path, alignment, reference, idxstats=idxstats
    )
    assert result.returncode == 0
    assert result.stdout.strip() == "true"
    assert audit["decision"] == "normalize"
    assert audit["differences"]["extra_contigs"] == [
        {"name": "virus", "length": 40, "mapped_reads": 7, "unmapped_reads": 1}
    ]
    assert audit["differences"]["missing_reference_contigs"] == [
        {"name": "chr2", "length": 80}
    ]
    assert audit["normalization"]["candidate_reads_affected"] == 8


def test_strict_policy_rejects_any_dictionary_difference(tmp_path):
    result, audit = run_checker(
        tmp_path,
        [("chr1", 100), ("virus", 40)],
        [("chr1", 100)],
        policy="strict",
    )
    assert result.returncode == 2
    assert audit["decision"] == "rejected"
    assert "policy is strict" in audit["reason"]


def test_shared_contig_length_conflict_is_never_normalized(tmp_path):
    result, audit = run_checker(
        tmp_path,
        [("chr1", 99)],
        [("chr1", 100)],
    )
    assert result.returncode == 2
    assert audit["decision"] == "rejected"
    assert audit["differences"]["shared_length_mismatches"] == [
        {"contig": "chr1", "alignment_length": 99, "reference_length": 100}
    ]


def test_order_only_difference_requires_normalization(tmp_path):
    result, audit = run_checker(
        tmp_path,
        [("chr2", 80), ("chr1", 100)],
        [("chr1", 100), ("chr2", 80)],
    )
    assert result.returncode == 0
    assert result.stdout.strip() == "true"
    assert audit["differences"]["shared_order_mismatch"] is True


def test_malformed_dictionary_is_rejected(tmp_path):
    result, audit = run_checker(tmp_path, [], [("chr1", 100)])
    assert result.returncode == 2
    assert audit["decision"] == "rejected"
    assert audit["differences"]["malformed_records"]
