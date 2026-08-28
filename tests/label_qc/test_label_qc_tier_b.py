"""Tests for bin/label_qc.py — Tier B (BAM-based verification via samtools).

Run with: .venv/bin/python -m pytest tests/label_qc/test_label_qc_tier_b.py -v

Tier B verifies Tier A labels against the BAMs without any pysam dependency:
  S4  normal contamination — alt-VAF in the DNA normal at Somatic truth sites,
      measured by a single batched `samtools mpileup -l sites` subprocess
  B1  strand/orientation bias — per-site strand counts of alt reads in the DNA
      tumor BAM (variant-level flag; optional WARN gate S8)

Fixtures are tiny synthetic SAMs compiled to BAM in-test with samtools, so the
suite is hermetic. Tests needing samtools are skipped when it is not on PATH.
"""

import gzip
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from test_label_qc import (CLI, CONFIG, clean_records, make_cohort,  # noqa: F401
                           read_flagged, read_summary, read_verdicts, rec,
                           run_cli, somatic, write_vcf)

SAMTOOLS = shutil.which("samtools")
requires_samtools = pytest.mark.skipif(SAMTOOLS is None,
                                       reason="samtools not available")

REF_LEN = 10000
SITE_STEP = 300
SITE_START = 1000
READ_LEN = 100
ALT_OFFSET = 50  # 0-based offset of the variant base inside each read


def site_positions(n):
    return [SITE_START + i * SITE_STEP for i in range(n)]


def make_reference(tmp_path):
    ref = Path(tmp_path) / "ref.fa"
    ref.write_text(">chr1\n" + "A" * REF_LEN + "\n")
    subprocess.run([SAMTOOLS, "faidx", str(ref)], check=True,
                   capture_output=True)
    return ref


def site_reads(pos, n_alt_fwd, n_alt_rev, n_ref_fwd, n_ref_rev, prefix):
    """Reads of READ_LEN with the variant base at ALT_OFFSET covering pos.

    Per the SAM spec, SEQ is stored reverse-complemented when flag 16 is set,
    i.e. in forward orientation relative to the reference — so the aligned
    sequence is written identically for both strands; mpileup encodes strand
    by letter case.
    """
    reads = []
    start = pos - ALT_OFFSET
    for i in range(n_alt_fwd + n_alt_rev + n_ref_fwd + n_ref_rev):
        if i < n_alt_fwd:
            flag, base = 0, "G"
        elif i < n_alt_fwd + n_alt_rev:
            flag, base = 16, "G"
        elif i < n_alt_fwd + n_alt_rev + n_ref_fwd:
            flag, base = 0, "A"
        else:
            flag, base = 16, "A"
        seq = "A" * ALT_OFFSET + base + "A" * (READ_LEN - ALT_OFFSET - 1)
        reads.append(("%s_%d_%d" % (prefix, pos, i), flag, start, seq))
    return reads


def make_bam(tmp_path, name, reads):
    """Compile a read list [(qname, flag, pos, seq)] into an indexed BAM."""
    tmp_path = Path(tmp_path)
    sam = tmp_path / (name + ".sam")
    lines = ["@HD\tVN:1.6\tSO:coordinate",
             "@SQ\tSN:chr1\tLN:%d" % REF_LEN]
    for qname, flag, pos, seq in reads:
        lines.append("\t".join([qname, str(flag), "chr1", str(pos), "60",
                                "%dM" % len(seq), "*", "0", "0", seq,
                                "I" * len(seq)]))
    sam.write_text("\n".join(lines) + "\n")
    unsorted = tmp_path / (name + ".unsorted.bam")
    bam = tmp_path / (name + ".bam")
    subprocess.run([SAMTOOLS, "view", "-b", "-o", str(unsorted), str(sam)],
                   check=True, capture_output=True)
    subprocess.run([SAMTOOLS, "sort", "-o", str(bam), str(unsorted)],
                   check=True, capture_output=True)
    subprocess.run([SAMTOOLS, "index", str(bam)], check=True, capture_output=True)
    unsorted.unlink()
    return bam


def truth_vcf(tmp_path, n_sites=25):
    return write_vcf(Path(tmp_path) / "truth.vcf",
                     [somatic(p) for p in site_positions(n_sites)])


def write_bam_manifest(path, rows):
    """rows: (sample_id, normal_bam, tumor_bam, truth_vcf) — ClairS header."""
    lines = ["#sample_id\tnormal_bam\ttumor_bam\ttruth_vcf"]
    lines += ["\t".join(str(c) for c in r) for r in rows]
    Path(path).write_text("\n".join(lines) + "\n")
    return path


def gate_status(outdir, sid, gate):
    summary = read_summary(outdir)
    for s in summary["samples"]:
        if s["sample"] == sid:
            return {g["gate"]: g["status"] for g in s["gates"]}.get(gate)
    return None


# ---------------------------------------------------------------- S4: normal contamination

@requires_samtools
class TestS4NormalContamination:
    def test_contaminated_normal_fails(self, tmp_path):
        """Normal with 80% alt reads at every truth site (4081 signature)."""
        make_reference(tmp_path)
        vcf = truth_vcf(tmp_path)
        reads = []
        for pos in site_positions(25):
            reads += site_reads(pos, 8, 8, 2, 2, "dn")  # VAF 0.8
        normal = make_bam(tmp_path, "normal", reads)
        manifest = write_bam_manifest(tmp_path / "m.tsv",
                                      [("contam", normal, "", vcf)])
        out = tmp_path / "out"
        run_cli("--samples", manifest, "--out", out, "--verify-bam")
        v = read_verdicts(out)
        assert v["contam"]["verdict"] == "FAIL"
        assert "S4" in v["contam"]["gates_failed"]
        assert float(v["contam"]["normal_contam_fraction"]) >= 0.9

    def test_clean_normal_passes(self, tmp_path):
        make_reference(tmp_path)
        vcf = truth_vcf(tmp_path)
        reads = []
        for pos in site_positions(25):
            reads += site_reads(pos, 0, 0, 10, 10, "dn")  # all ref
        normal = make_bam(tmp_path, "normal", reads)
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "clean", "--out", out,
                "--verify-bam", "--normal-bam", normal)
        v = read_verdicts(out)
        assert v["clean"]["verdict"] == "PASS"
        assert gate_status(out, "clean", "S4") == "PASS"
        assert float(v["clean"]["normal_contam_fraction"]) == 0.0

    def test_fail_fraction_threshold_override(self, tmp_path):
        """Raising fail_fraction above the measured fraction demotes to WARN."""
        make_reference(tmp_path)
        cfg = json.loads(CONFIG.read_text())
        cfg["tier_b"]["normal_contamination"]["fail_fraction"] = 1.5
        cfg_path = tmp_path / "cfg.json"
        cfg_path.write_text(json.dumps(cfg))
        vcf = truth_vcf(tmp_path)
        reads = []
        for pos in site_positions(25):
            reads += site_reads(pos, 8, 8, 2, 2, "dn")
        normal = make_bam(tmp_path, "normal", reads)
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out,
                "--verify-bam", "--normal-bam", normal, "--config", cfg_path)
        v = read_verdicts(out)
        assert "S4" in v["s1"]["gates_warned"]
        assert v["s1"]["verdict"] == "WARN"

    def test_verify_bam_without_bam_skips_s4(self, tmp_path):
        vcf = truth_vcf(tmp_path)
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "nobam", "--out", out,
                "--verify-bam")
        v = read_verdicts(out)
        assert v["nobam"]["verdict"] == "PASS"
        assert gate_status(out, "nobam", "S4") == "SKIP"

    def test_samtools_missing_skips_s4(self, tmp_path):
        vcf = truth_vcf(tmp_path)
        normal = tmp_path / "fake.bam"
        normal.write_bytes(b"not a bam")
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out,
                "--verify-bam", "--normal-bam", normal,
                "--samtools", "/nonexistent/samtools")
        v = read_verdicts(out)
        assert v["s1"]["verdict"] == "PASS"
        assert gate_status(out, "s1", "S4") == "SKIP"


# ---------------------------------------------------------------- B1: strand bias

@requires_samtools
class TestStrandBias:
    def _tumor_bam(self, tmp_path, biased_pos):
        reads = []
        for pos in site_positions(25):
            if pos == biased_pos:
                reads += site_reads(pos, 10, 0, 5, 5, "dt")  # alt all forward
            else:
                reads += site_reads(pos, 5, 5, 5, 5, "dt")   # balanced
        return make_bam(tmp_path, "tumor", reads)

    def test_strand_biased_site_flagged(self, tmp_path):
        make_reference(tmp_path)
        vcf = truth_vcf(tmp_path)
        tumor = self._tumor_bam(tmp_path, SITE_START)
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "sb", "--out", out,
                "--verify-bam", "--tumor-bam", tumor)
        rows = read_flagged(out)
        b1 = [r for r in rows if "B1" in r["rules"].split("+")]
        assert len(b1) == 1
        assert b1[0]["pos"] == str(SITE_START)
        assert b1[0]["confidence"] == "low"
        assert b1[0]["action"] == "REPORT_ONLY"
        v = read_verdicts(out)
        assert v["sb"]["n_strand_bias_sites"] == "1"
        # strand bias alone must not change the verdict (S8 disabled by default)
        assert v["sb"]["verdict"] == "PASS"

    def test_balanced_strands_not_flagged(self, tmp_path):
        make_reference(tmp_path)
        vcf = truth_vcf(tmp_path)
        reads = []
        for pos in site_positions(25):
            reads += site_reads(pos, 5, 5, 5, 5, "dt")
        tumor = make_bam(tmp_path, "tumor", reads)
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "bal", "--out", out,
                "--verify-bam", "--tumor-bam", tumor)
        assert read_flagged(out) == []
        v = read_verdicts(out)
        assert v["bal"]["n_strand_bias_sites"] == "0"

    def test_strand_share_threshold_override(self, tmp_path):
        """Raising max_strand_share to 1.0 + epsilon disables the flag."""
        make_reference(tmp_path)
        cfg = json.loads(CONFIG.read_text())
        cfg["tier_b"]["strand_bias"]["max_strand_share"] = 1.01
        cfg_path = tmp_path / "cfg.json"
        cfg_path.write_text(json.dumps(cfg))
        vcf = truth_vcf(tmp_path)
        tumor = self._tumor_bam(tmp_path, SITE_START)
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "sb", "--out", out,
                "--verify-bam", "--tumor-bam", tumor, "--config", cfg_path)
        assert read_flagged(out) == []


# ---------------------------------------------------------------- Tier A isolation

class TestTierAIsolation:
    def test_no_bam_output_byte_identical(self, tmp_path):
        """Without BAMs, --verify-bam output equals plain Tier A output."""
        manifest = make_cohort(tmp_path,
                               [("a", clean_records(30)), ("b", clean_records(30))])
        out_a = tmp_path / "out_a"
        out_b = tmp_path / "out_b"
        run_cli("--samples", manifest, "--out", out_a)
        run_cli("--samples", manifest, "--out", out_b, "--verify-bam")
        assert (out_a / "samples_qc.tsv").read_bytes() == \
               (out_b / "samples_qc.tsv").read_bytes()
        with gzip.open(out_a / "flagged_sites.tsv.gz", "rt") as fh:
            fa = fh.read()
        with gzip.open(out_b / "flagged_sites.tsv.gz", "rt") as fh:
            fb = fh.read()
        assert fa == fb

    def test_default_run_has_no_tier_b_columns(self, tmp_path):
        manifest = make_cohort(tmp_path, [("a", clean_records(30))])
        out = tmp_path / "out"
        run_cli("--samples", manifest, "--out", out)
        header = (out / "samples_qc.tsv").read_text().split("\n")[0]
        assert "normal_contam_fraction" not in header
        assert "n_strand_bias_sites" not in header
        summary = read_summary(out)
        assert "tier_b" not in summary["samples"][0]

    def test_bundled_config_has_tier_b_section(self):
        cfg = json.loads(CONFIG.read_text())
        tb = cfg["tier_b"]
        assert tb["enabled"] is False
        assert "samtools" in tb
        nc = tb["normal_contamination"]
        assert nc["alt_vaf_min"] == 0.05
        assert nc["fail_fraction"] > nc["warn_fraction"] > 0
        sb = tb["strand_bias"]
        assert 0.5 < sb["max_strand_share"] <= 1.0
