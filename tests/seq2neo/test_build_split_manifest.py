"""Tests for examples/seq2neo/scripts/build_split_manifest.py.

Run with: .venv/bin/python -m pytest tests/seq2neo/ -v

Tests exercise the build on tiny synthetic gzipped VCFs + mini manifests in
tmp_path, following the tests/label_qc conventions (in-process main() calls,
pytest.raises(SystemExit) for build failures). Mini manifests always include
all RESERVED_SAMPLE_IDS rows so reserved validation passes in happy-path tests.
"""

import gzip
import importlib.util
import sys
from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_SCRIPT = _REPO_ROOT / "examples" / "seq2neo" / "scripts" / "build_split_manifest.py"
_spec = importlib.util.spec_from_file_location("build_split_manifest", _SCRIPT)
bsm = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = bsm  # required for multiprocessing pickling
_spec.loader.exec_module(bsm)

VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

MANIFEST_COLUMNS = [
    "sample_id", "project_id", "patient_id", "set_number", "disease_normalized",
    "label_qc_verdict", "training_label_vcf", "known_limitations",
]

RESERVED_IDS = sorted(bsm.RESERVED_SAMPLE_IDS)


def rec(chrom, pos, ref, alt, filt, info=None):
    info_s = ";".join(f"{k}={v}" for k, v in (info or {}).items()) or "."
    return "\t".join([chrom, str(pos), ".", ref, alt, ".", filt, info_s])


def write_vcf(path, records):
    with gzip.open(path, "wt") as fh:
        fh.write(VCF_HEADER)
        for r in records:
            fh.write(r + "\n")
    return path


def make_row(sid, verdict, vcf, disease="colon cancer"):
    return {
        "sample_id": sid,
        "project_id": "PRJNA999999",
        "patient_id": sid.rsplit("_", 1)[-1],
        "set_number": "1",
        "disease_normalized": disease,
        "label_qc_verdict": verdict,
        "training_label_vcf": str(vcf),
        "known_limitations": bsm.KNOWN_LIMITATIONS,
    }


def reserved_rows(tmp_path):
    """All 5 reserved samples, PASS verdicts, expected diseases, tiny VCFs."""
    rows = []
    for sid in RESERVED_IDS:
        vcf = write_vcf(tmp_path / f"{sid}.vcf.gz", [
            rec("chr2", 100, "A", "T", "Somatic", {"RESCUED": "NO"}),
        ])
        rows.append(make_row(sid, "PASS", vcf,
                             disease=bsm._EXPECTED_RESERVED_DISEASES[sid]))
    return rows


def write_manifest(path, rows):
    text = ("\t".join(MANIFEST_COLUMNS) + "\n"
            + "".join("\t".join(r[c] for c in MANIFEST_COLUMNS) + "\n"
                      for r in rows))
    Path(path).write_text(text)
    return path


def run_build(manifest, outdir, *extra):
    return bsm.main(["--manifest", str(manifest), "--outdir", str(outdir),
                     "--workers", "1", *extra])


def build_cohort(tmp_path, extra_rows):
    """Happy-path cohort: 5 reserved + caller-supplied extra rows."""
    rows = reserved_rows(tmp_path) + extra_rows
    manifest = write_manifest(tmp_path / "manifest.tsv", rows)
    outdir = tmp_path / "out"
    assert run_build(manifest, outdir) == 0
    return manifest, outdir


def read_parquet(outdir):
    return pl.read_parquet(Path(outdir) / "selected_variants.parquet")


# ── whitelist filtering ──────────────────────────────────────────────────────

def test_whitelist_filtering(tmp_path):
    vcf = write_vcf(tmp_path / "S1.vcf.gz", [
        rec("chr2", 10, "A", "T", "Somatic"),
        rec("chr2", 20, "A", "T", "Germline"),
        rec("chr2", 30, "A", "T", "Reference"),
        rec("chr2", 40, "A", "T", "Artifact"),
        rec("chr2", 50, "A", "T", "NoConsensus"),
        rec("chr2", 60, "A", "G", "RNAedit"),
    ])
    _, outdir = build_cohort(tmp_path, [make_row("S1", "PASS", vcf)])
    df = read_parquet(outdir).filter(pl.col("sample_id") == "S1")
    assert sorted(df["FILTER"].to_list()) == ["Germline", "Reference", "Somatic"]

    # sample-level TSV: n_variants_selected counts only whitelist rows
    tsv = (Path(outdir) / "sample_split.tsv").read_text().splitlines()
    row = [l for l in tsv if l.startswith("S1\t")][0].split("\t")
    assert row[6] == "train_pool"
    assert row[8] == "3"


# ── reserved routing ─────────────────────────────────────────────────────────

def test_reserved_samples_route_to_reserved(tmp_path):
    _, outdir = build_cohort(tmp_path, [])
    df = read_parquet(outdir)
    reserved = df.filter(pl.col("sample_id").is_in(RESERVED_IDS))
    assert reserved.height == 5
    assert set(reserved["split"].to_list()) == {"reserved"}
    assert reserved["is_zero_shot"].all()
    others = df.filter(~pl.col("sample_id").is_in(RESERVED_IDS))
    assert (~others["is_zero_shot"]).all() if others.height else True

    tsv = (Path(outdir) / "sample_split.tsv").read_text().splitlines()
    assert len(tsv) == 1 + 5  # header + 5 reserved rows, sorted by sample_id
    ids = [l.split("\t")[0] for l in tsv[1:]]
    assert ids == sorted(ids)
    for l in tsv[1:]:
        assert l.split("\t")[6] == "reserved"


# ── reserved validation failures ─────────────────────────────────────────────

def test_reserved_validation_missing_id(tmp_path, capsys):
    rows = reserved_rows(tmp_path)[1:]  # drop one reserved sample
    manifest = write_manifest(tmp_path / "manifest.tsv", rows)
    with pytest.raises(SystemExit) as exc:
        run_build(manifest, tmp_path / "out")
    assert exc.value.code != 0
    err = capsys.readouterr().err
    assert RESERVED_IDS[0] in err
    assert "reserved" in err


def test_reserved_validation_warn_verdict(tmp_path, capsys):
    rows = reserved_rows(tmp_path)
    rows[0]["label_qc_verdict"] = "WARN"
    manifest = write_manifest(tmp_path / "manifest.tsv", rows)
    with pytest.raises(SystemExit) as exc:
        run_build(manifest, tmp_path / "out")
    assert exc.value.code != 0
    assert "label_qc_verdict" in capsys.readouterr().err


def test_reserved_validation_wrong_disease(tmp_path, capsys):
    rows = reserved_rows(tmp_path)
    rows[0]["disease_normalized"] = "colon cancer"
    manifest = write_manifest(tmp_path / "manifest.tsv", rows)
    with pytest.raises(SystemExit) as exc:
        run_build(manifest, tmp_path / "out")
    assert exc.value.code != 0
    assert "disease_normalized" in capsys.readouterr().err


# ── other build failures ─────────────────────────────────────────────────────

def test_missing_vcf_fails(tmp_path, capsys):
    rows = reserved_rows(tmp_path) + [
        make_row("S1", "PASS", tmp_path / "does_not_exist.vcf.gz")]
    manifest = write_manifest(tmp_path / "manifest.tsv", rows)
    with pytest.raises(SystemExit) as exc:
        run_build(manifest, tmp_path / "out")
    assert exc.value.code != 0
    assert "missing/unreadable" in capsys.readouterr().err


def test_duplicate_sample_id_fails(tmp_path, capsys):
    rows = reserved_rows(tmp_path)
    dup = write_vcf(tmp_path / "dup.vcf.gz", [rec("chr2", 1, "A", "T", "Somatic")])
    rows += [make_row("S1", "PASS", dup), make_row("S1", "PASS", dup)]
    manifest = write_manifest(tmp_path / "manifest.tsv", rows)
    with pytest.raises(SystemExit) as exc:
        run_build(manifest, tmp_path / "out")
    assert exc.value.code != 0
    assert "duplicate" in capsys.readouterr().err


# ── split routing ────────────────────────────────────────────────────────────

def test_warn_sample_uses_chromosome_split(tmp_path):
    vcf = write_vcf(tmp_path / "W1.vcf.gz", [
        rec("chr1", 10, "A", "T", "Somatic"),    # would be test for PASS
        rec("chr21", 20, "A", "T", "Germline"),  # would be val for PASS
        rec("chr22", 30, "A", "T", "Reference"),
    ])
    _, outdir = build_cohort(tmp_path, [make_row("W1", "WARN", vcf)])
    df = read_parquet(outdir).filter(pl.col("sample_id") == "W1")
    assert df.height == 3
    by_pos = dict(zip(df["POS"].to_list(), df["split"].to_list()))
    assert by_pos == {10: "test", 20: "val", 30: "val"}
    assert set(df["label_verdict"].to_list()) == {"WARN"}


def test_pass_sample_chromosome_split(tmp_path):
    vcf = write_vcf(tmp_path / "P1.vcf.gz", [
        rec("chr1", 10, "A", "T", "Somatic"),
        rec("chr21", 20, "A", "T", "Somatic"),
        rec("chr22", 30, "A", "T", "Somatic"),
        rec("chr2", 40, "A", "T", "Somatic"),
        rec("chr20", 50, "A", "T", "Somatic"),
    ])
    _, outdir = build_cohort(tmp_path, [make_row("P1", "PASS", vcf)])
    df = read_parquet(outdir).filter(pl.col("sample_id") == "P1")
    by_pos = dict(zip(df["POS"].to_list(), df["split"].to_list()))
    assert by_pos == {10: "test", 20: "val", 30: "val", 40: "train", 50: "train"}


def test_chromosome_normalization(tmp_path):
    vcf = write_vcf(tmp_path / "N1.vcf.gz", [
        rec("1", 10, "A", "T", "Somatic"),            # → chr1 → test
        rec("Chr1", 20, "A", "T", "Somatic"),         # → chr1 → test
        rec("chromosome1", 30, "A", "T", "Somatic"),  # → chr1 → test
        rec("chrX", 40, "A", "T", "Somatic"),         # → train
        rec("chrY", 50, "A", "T", "Somatic"),         # → train
        rec("chrM", 60, "A", "T", "Somatic"),         # → train
        rec("chrUn_gl000220", 70, "A", "T", "Somatic"),  # unrecognized → train
    ])
    _, outdir = build_cohort(tmp_path, [make_row("N1", "PASS", vcf)])
    df = read_parquet(outdir).filter(pl.col("sample_id") == "N1")
    by_pos = dict(zip(df["POS"].to_list(), df["split"].to_list()))
    assert by_pos[10] == by_pos[20] == by_pos[30] == "test"
    assert all(by_pos[p] == "train" for p in (40, 50, 60, 70))


# ── sub-pool tag edge cases ──────────────────────────────────────────────────

def test_tag_edge_cases(tmp_path):
    vcf = write_vcf(tmp_path / "T1.vcf.gz", [
        # missing VAF_RNA_MEAN → rna_nonzero False → is_low_vaf_a False
        rec("chr2", 10, "A", "T", "Somatic",
            {"VAF_DNA_MEAN": 0.10, "DP_RNA_MEAN": 30, "RESCUED": "NO"}),
        # low VAF-A with RNA evidence → True
        rec("chr2", 20, "A", "T", "Somatic",
            {"VAF_DNA_MEAN": 0.10, "VAF_RNA_MEAN": 0.2, "DP_RNA_MEAN": 30,
             "RESCUED": "NO"}),
        # low VAF-B band → is_low_vaf_b True, is_low_vaf_a False
        rec("chr2", 30, "A", "T", "Somatic",
            {"VAF_DNA_MEAN": 0.20, "VAF_RNA_MEAN": 0.2, "DP_RNA_MEAN": 30,
             "RESCUED": "NO"}),
        # missing DP_DNA_MEAN → is_low_dp False (null semantics)
        rec("chr2", 40, "A", "T", "Somatic",
            {"VAF_DNA_MEAN": 0.5, "VAF_RNA_MEAN": 0.2, "DP_RNA_MEAN": 30,
             "RESCUED": "NO"}),
        # DP_DNA_MEAN < 20 with RNA evidence → is_low_dp True
        rec("chr2", 50, "A", "T", "Somatic",
            {"VAF_DNA_MEAN": 0.5, "DP_DNA_MEAN": 10, "VAF_RNA_MEAN": 0.2,
             "DP_RNA_MEAN": 30, "RESCUED": "YES"}),
        # missing RESCUED → both rescue tags False; also an indel
        rec("chr2", 60, "AT", "A", "Germline", {"VAF_DNA_MEAN": 0.5}),
    ])
    _, outdir = build_cohort(tmp_path, [make_row("T1", "PASS", vcf)])
    df = read_parquet(outdir).filter(pl.col("sample_id") == "T1")
    rows = {r["POS"]: r for r in df.iter_rows(named=True)}

    assert rows[10]["is_low_vaf_a"] is False
    assert rows[10]["VAF_RNA_MEAN"] is None
    assert rows[20]["is_low_vaf_a"] is True
    assert rows[20]["is_low_vaf_b"] is False
    assert rows[30]["is_low_vaf_a"] is False
    assert rows[30]["is_low_vaf_b"] is True
    assert rows[40]["is_low_dp"] is False
    assert rows[40]["DP_DNA_MEAN"] is None
    assert rows[50]["is_low_dp"] is True
    assert rows[50]["is_rescued"] is True
    assert rows[50]["is_non_rescued"] is False
    assert rows[10]["is_rescued"] is False
    assert rows[10]["is_non_rescued"] is True
    assert rows[60]["is_rescued"] is False
    assert rows[60]["is_non_rescued"] is False
    assert rows[60]["RESCUED"] is None
    assert rows[60]["is_indel"] is True
    assert rows[60]["variant_type"] == "DEL"
    assert rows[10]["variant_type"] == "SNV"
    assert rows[10]["is_indel"] is False
    assert (~df["is_zero_shot"]).all()


def test_variant_type_insertion_and_other(tmp_path):
    vcf = write_vcf(tmp_path / "V1.vcf.gz", [
        rec("chr2", 10, "A", "AT", "Somatic"),     # INS
        rec("chr2", 20, "AT", "GC", "Somatic"),    # same len > 1 → OTHER
    ])
    _, outdir = build_cohort(tmp_path, [make_row("V1", "PASS", vcf)])
    df = read_parquet(outdir).filter(pl.col("sample_id") == "V1")
    by_pos = {r["POS"]: r for r in df.iter_rows(named=True)}
    assert by_pos[10]["variant_type"] == "INS"
    assert by_pos[10]["is_indel"] is True
    assert by_pos[20]["variant_type"] == "OTHER"
    assert by_pos[20]["is_indel"] is False


# ── determinism ──────────────────────────────────────────────────────────────

def _two_extra_samples(tmp_path):
    vcf_a = write_vcf(tmp_path / "A.vcf.gz", [
        rec("chr2", 10, "A", "T", "Somatic", {"RESCUED": "YES"}),
        rec("chr1", 20, "A", "T", "Germline", {"RESCUED": "NO"}),
    ])
    vcf_b = write_vcf(tmp_path / "B.vcf.gz", [
        rec("chr21", 10, "A", "T", "Reference"),
        rec("chr1", 30, "A", "T", "Somatic"),   # WARN → train
    ])
    return [make_row("AAA", "PASS", vcf_a), make_row("BBB", "WARN", vcf_b)]


def test_determinism(tmp_path):
    rows = reserved_rows(tmp_path) + _two_extra_samples(tmp_path)
    manifest = write_manifest(tmp_path / "manifest.tsv", rows)
    out1, out2 = tmp_path / "out1", tmp_path / "out2"
    assert run_build(manifest, out1) == 0
    assert run_build(manifest, out2) == 0
    assert (out1 / "sample_split.tsv").read_bytes() == \
           (out2 / "sample_split.tsv").read_bytes()
    assert_frame_equal(read_parquet(out1), read_parquet(out2))


def test_workers_parallel_matches_sequential(tmp_path):
    rows = reserved_rows(tmp_path) + _two_extra_samples(tmp_path)
    manifest = write_manifest(tmp_path / "manifest.tsv", rows)
    out_seq, out_par = tmp_path / "seq", tmp_path / "par"
    assert run_build(manifest, out_seq) == 0
    assert bsm.main(["--manifest", str(manifest), "--outdir", str(out_par),
                     "--workers", "4"]) == 0
    assert (out_seq / "sample_split.tsv").read_bytes() == \
           (out_par / "sample_split.tsv").read_bytes()
    assert_frame_equal(read_parquet(out_seq), read_parquet(out_par))


# ── --check mode ─────────────────────────────────────────────────────────────

def test_check_up_to_date(tmp_path):
    manifest, outdir = build_cohort(tmp_path, _two_extra_samples(tmp_path))
    assert run_build(manifest, outdir, "--check") == 0


def test_check_stale_on_modified_output(tmp_path, capsys):
    manifest, outdir = build_cohort(tmp_path, _two_extra_samples(tmp_path))
    tsv = outdir / "sample_split.tsv"
    tsv.write_text(tsv.read_text().replace("train_pool", "train_poolX", 1))
    with pytest.raises(SystemExit) as exc:
        run_build(manifest, outdir, "--check")
    assert exc.value.code == 1
    assert "STALE" in capsys.readouterr().err


def test_check_stale_on_missing_outputs(tmp_path, capsys):
    rows = reserved_rows(tmp_path) + _two_extra_samples(tmp_path)
    manifest = write_manifest(tmp_path / "manifest.tsv", rows)
    with pytest.raises(SystemExit) as exc:
        run_build(manifest, tmp_path / "empty_out", "--check")
    assert exc.value.code == 1
    assert "STALE" in capsys.readouterr().err
