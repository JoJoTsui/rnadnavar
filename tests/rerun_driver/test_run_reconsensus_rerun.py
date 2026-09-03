"""Tests for examples/seq2neo/scripts/run_reconsensus_rerun.py — rerun driver.

Run with: .venv/bin/python -m pytest tests/rerun_driver/ -v

Tests build a tiny fake cohort tree (2 samples x 6 caller VCFs) under tmp_path
and exercise the driver CLI seam: dry-run output, samplesheet generation,
checksum guard record/verify, and the no-caller-tools / new-outdir guards.
"""

import gzip
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DRIVER = _REPO_ROOT / "examples" / "seq2neo" / "scripts" / "run_reconsensus_rerun.py"

spec = importlib.util.spec_from_file_location("run_reconsensus_rerun", DRIVER)
rr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rr)

CALLER_SUFFIX = rr.CALLER_VCF_SUFFIX  # {caller: suffix}
MODALITIES = rr.MODALITY_CONFIG       # {modality: {status, pair, subdir}}

MANIFEST_FIELDS = [
    "sample_id", "project_id", "patient_id", "set_number", "disease",
    "status", "base_output_dir", "dir_name", "vcf_prefix", "is_complete",
]


def _tiny_vcf_gz(path: Path, body: str = "chr1\t2\t.\tC\tT\t50\tPASS\t.\n"):
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write(header + body)


def make_sample(source_root: Path, sample_id: str, skip: str = None) -> dict:
    """Create the 6 caller VCFs for one sample; return its manifest row."""
    prefix = sample_id
    for modality, mcfg in MODALITIES.items():
        pair = mcfg["pair"].format(p=prefix)
        for caller, suffix in CALLER_SUFFIX.items():
            if skip == f"{modality}_{caller}":
                continue
            vcf = (
                source_root / sample_id / mcfg["subdir"] / caller / pair
                / f"{pair}{suffix}"
            )
            _tiny_vcf_gz(vcf, body=f"chr1\t2\t.\tC\t{'T' if modality == 'dna' else 'A'}\t50\tPASS\t.\n")
    return {
        "sample_id": sample_id,
        "project_id": sample_id.split("_")[0],
        "patient_id": sample_id.split("_")[-1],
        "set_number": "2",
        "disease": "Colon cancer",
        "status": "standard",
        "base_output_dir": str(source_root),
        "dir_name": sample_id,
        "vcf_prefix": prefix,
        "is_complete": "True",
    }


def write_manifest(path: Path, rows: list):
    path.parent.mkdir(parents=True, exist_ok=True)
    import csv

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def fake_cohort(tmp_path):
    """2-sample fake cohort tree + manifest; returns (manifest_path, rows)."""
    source_root = tmp_path / "output"
    rows = [
        make_sample(source_root, "PRJNA298330_4032"),
        make_sample(source_root, "PRJNA298376_3942"),
    ]
    manifest = tmp_path / "data" / "processed" / "sample_manifest.tsv"
    write_manifest(manifest, rows)
    return manifest, rows


def run_driver(*argv, cwd=None):
    return subprocess.run(
        [sys.executable, str(DRIVER), *argv],
        capture_output=True,
        text=True,
        cwd=cwd,
    )


# ---------------------------------------------------------------------------
# Dry-run
# ---------------------------------------------------------------------------


def test_dry_run_prints_plan_and_writes_nothing(fake_cohort, tmp_path):
    manifest, _ = fake_cohort
    rerun_root = tmp_path / "rerun"
    res = run_driver(
        "--manifest", str(manifest),
        "--seq2neo", str(rerun_root),
        "--dry-run",
    )
    assert res.returncode == 0, res.stderr
    out = res.stdout
    assert "Selected 2 sample(s)" in out
    # consensus-step entry, no caller names in tools
    assert "--step consensus" in out
    assert "--tools consensus,rescue,filtering,vep" in out
    # samplesheet preview: 6 VCF rows per sample, no fastq/bam/cram
    assert out.count("patient,sample,status,variantcaller,vcf") == 2
    assert "PRJNA298330_4032DT_vs_PRJNA298330_4032DN,1,mutect2" in out
    assert "PRJNA298330_4032RT_realign_vs_PRJNA298330_4032DN,2,strelka" in out
    for forbidden in ("fastq", ".bam", ".cram"):
        assert forbidden not in out
    # checksums printed
    assert "input md5" in out
    # dry-run writes nothing
    assert not (rerun_root / "runs").exists()
    assert not (rerun_root / "output_reconsensus").exists()


def test_dry_run_reports_missing_inputs(fake_cohort, tmp_path):
    manifest, rows = fake_cohort
    # Select the asserted modality explicitly; filesystem traversal order is
    # not a contract and may return the RNA file first.
    victim = Path(rr.locate_caller_vcfs(rows[0])["dna_deepsomatic"])
    victim.unlink()
    res = run_driver(
        "--manifest", str(manifest),
        "--seq2neo", str(tmp_path / "rerun"),
        "--sample", rows[0]["sample_id"],
        "--dry-run",
    )
    assert res.returncode == 0, res.stderr
    assert "missing caller VCFs" in res.stdout
    assert "dna_deepsomatic" in res.stdout


# ---------------------------------------------------------------------------
# Checksum guard (function-level)
# ---------------------------------------------------------------------------


def test_checksum_guard_roundtrip_and_tamper_detection(tmp_path):
    source_root = tmp_path / "output"
    row = make_sample(source_root, "TEST_1")
    inputs = rr.locate_caller_vcfs(row)
    assert sorted(inputs) == sorted(rr.expected_inputs())

    manifest_path = tmp_path / "checks" / "TEST_1.input_checksums.json"
    manifest = rr.record_checksums(inputs, manifest_path)
    # VCFs plus any indices are recorded
    assert len(manifest["files"]) >= 6

    ok, problems = rr.verify_checksums(manifest_path)
    assert ok and not problems

    # tamper with one input -> guard must trip
    victim = inputs["dna_mutect2"]
    with open(victim, "ab") as fh:
        fh.write(b"tampered")
    ok, problems = rr.verify_checksums(manifest_path)
    assert not ok
    assert any("checksum mismatch" in p and str(victim) in p for p in problems)

    # deleted input -> guard must trip
    victim2 = inputs["rna_strelka"]
    victim2.unlink()
    ok, problems = rr.verify_checksums(manifest_path)
    assert not ok
    assert any("missing input" in p for p in problems)


# ---------------------------------------------------------------------------
# Samplesheet generation
# ---------------------------------------------------------------------------


def test_samplesheet_layout(tmp_path):
    source_root = tmp_path / "output"
    row = make_sample(source_root, "PRJNA298330_4032")
    inputs = rr.locate_caller_vcfs(row)
    csv_path = tmp_path / "csv" / "x.csv"
    rr.write_sample_csv(row, inputs, csv_path)
    lines = csv_path.read_text().strip().splitlines()
    assert lines[0] == "patient,sample,status,variantcaller,vcf"
    assert len(lines) == 7  # 3 DNA + 3 RNA
    dna = [l for l in lines[1:] if ",1," in l]
    rna = [l for l in lines[1:] if ",2," in l]
    assert len(dna) == 3 and len(rna) == 3
    assert all("DT_vs_" in l for l in dna)
    assert all("RT_realign_vs_" in l for l in rna)
    # every referenced file exists
    for l in lines[1:]:
        assert Path(l.split(",")[-1]).is_file()


def test_samplesheet_numeric_vcf_prefix_uses_sample_id_as_patient(tmp_path):
    """Regression: a numeric-only vcf_prefix (e.g. "4278") is coerced to an
    integer by nf-schema and fails patient-is-string validation — the patient
    column must fall back to sample_id (PRJNA298376_4278-class samples)."""
    source_root = tmp_path / "output"
    row = make_sample(source_root, "PRJNA298376_4278")
    row["vcf_prefix"] = "4278"  # numeric-only prefix, as in the real manifest
    inputs = {f"{m}_{c}": "/x.vcf.gz" for m in MODALITIES for c in CALLER_SUFFIX}
    content = rr.render_sample_csv(row, inputs)
    lines = content.strip().splitlines()
    assert len(lines) == 7
    for l in lines[1:]:
        patient = l.split(",")[0]
        assert patient == "PRJNA298376_4278"
        assert not patient.isdigit()
        # sample column still uses the original vcf_prefix pair naming
        assert "4278DT_vs_4278DN" in l or "4278RT_realign_vs_4278DN" in l


# ---------------------------------------------------------------------------
# Config guards
# ---------------------------------------------------------------------------


def test_validate_config_rejects_caller_tools():
    cfg = dict(rr.DEFAULTS, main_nf="/x/main.nf", rdv_conf="/x/c.config",
               tools="consensus,rescue,mutect2")
    with pytest.raises(SystemExit):
        rr.validate_config(cfg)


def test_validate_config_rejects_non_consensus_step():
    cfg = dict(rr.DEFAULTS, main_nf="/x/main.nf", rdv_conf="/x/c.config",
               step="mapping")
    with pytest.raises(SystemExit):
        rr.validate_config(cfg)


def test_outdir_must_not_live_inside_source_root(tmp_path):
    source_root = tmp_path / "output"
    row = make_sample(source_root, "TEST_1")
    with pytest.raises(SystemExit):
        rr.check_outdir_separation([row], source_root / "rerun")
    # a sibling location is fine
    rr.check_outdir_separation([row], tmp_path / "output_reconsensus")


# ---------------------------------------------------------------------------
# Non-dry run without nextflow: missing inputs fail before execution
# ---------------------------------------------------------------------------


def test_missing_inputs_mark_failed_without_execution(fake_cohort, tmp_path):
    manifest, rows = fake_cohort
    victim = next(
        (Path(rows[0]["base_output_dir"]) / rows[0]["dir_name"]).rglob(
            "*.strelka.variants.vcf.gz"
        )
    )
    victim.unlink()
    rerun_root = tmp_path / "rerun"
    res = run_driver(
        "--manifest", str(manifest),
        "--seq2neo", str(rerun_root),
        "--sample", rows[0]["sample_id"],
        "--main-nf", "/nonexistent/main.nf",
        "--rdv-conf", "/nonexistent/c.config",
    )
    assert res.returncode == 0, res.stderr
    assert "missing caller VCFs" in res.stdout
    state = json.loads((rerun_root / "runs" / "rerun_state.json").read_text())
    assert state[rows[0]["sample_id"]]["status"] == "failed"
    assert "missing caller VCFs" in state[rows[0]["sample_id"]]["reason"]
    # nothing was executed: no outdir, no samplesheet, no checksum manifest
    assert not (rerun_root / "output_reconsensus").exists()
    assert not (rerun_root / "runs" / "rerun_csv").exists()
    assert not (rerun_root / "runs" / "rerun_checksums").exists()


# ---------------------------------------------------------------------------
# Verify-only mode
# ---------------------------------------------------------------------------


def test_verify_only_mode(fake_cohort, tmp_path):
    manifest, rows = fake_cohort
    rerun_root = tmp_path / "rerun"
    checksum_dir = rerun_root / "runs" / "rerun_checksums"
    for row in rows:
        inputs = rr.locate_caller_vcfs(row)
        rr.record_checksums(inputs, checksum_dir / f"{row['sample_id']}.input_checksums.json")

    res = run_driver(
        "--manifest", str(manifest),
        "--seq2neo", str(rerun_root),
        "--verify-only",
    )
    assert res.returncode == 0, res.stderr
    assert res.stdout.count("inputs untouched") == 2

    # tamper one input -> verify-only exits non-zero
    victim = rr.locate_caller_vcfs(rows[1])["dna_strelka"]
    with open(victim, "ab") as fh:
        fh.write(b"x")
    res = run_driver(
        "--manifest", str(manifest),
        "--seq2neo", str(rerun_root),
        "--verify-only",
    )
    assert res.returncode == 1
    assert "checksum mismatch" in res.stdout


# ---------------------------------------------------------------------------
# Parallel-cohort seams: comma-separated --sample, per-group --state-file
# ---------------------------------------------------------------------------


def test_sample_filter_accepts_comma_separated_list(fake_cohort, tmp_path):
    """The cohort launcher hands each parallel driver a disjoint comma-
    separated sample list via --sample."""
    manifest, rows = fake_cohort
    res = run_driver(
        "--manifest", str(manifest),
        "--seq2neo", str(tmp_path / "rerun"),
        "--sample", ",".join(r["sample_id"] for r in rows),
        "--dry-run",
    )
    assert res.returncode == 0, res.stderr
    assert "Selected 2 sample(s)" in res.stdout

    # a single id still works (backward compatible)
    res = run_driver(
        "--manifest", str(manifest),
        "--seq2neo", str(tmp_path / "rerun"),
        "--sample", rows[1]["sample_id"],
        "--dry-run",
    )
    assert res.returncode == 0, res.stderr
    assert "Selected 1 sample(s)" in res.stdout


def test_state_file_cli_override(fake_cohort, tmp_path):
    """--state-file redirects state writes (per-group state files keep
    concurrent drivers from racing on the shared state file)."""
    manifest, rows = fake_cohort
    victim = next(
        (Path(rows[0]["base_output_dir"]) / rows[0]["dir_name"]).rglob(
            "*.strelka.variants.vcf.gz"
        )
    )
    victim.unlink()
    rerun_root = tmp_path / "rerun"
    group_state = rerun_root / "runs" / "cohort_state" / "group1.json"
    res = run_driver(
        "--manifest", str(manifest),
        "--seq2neo", str(rerun_root),
        "--sample", rows[0]["sample_id"],
        "--state-file", str(group_state),
        "--main-nf", "/nonexistent/main.nf",
        "--rdv-conf", "/nonexistent/c.config",
    )
    assert res.returncode == 0, res.stderr
    # state written to the override path, not the default rerun_state.json
    assert group_state.is_file()
    state = json.loads(group_state.read_text())
    assert state[rows[0]["sample_id"]]["status"] == "failed"
    assert not (rerun_root / "runs" / "rerun_state.json").exists()
