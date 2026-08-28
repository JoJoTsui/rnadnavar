"""Tests for bin/label_qc.py — Tier A (VCF-only) label QC CLI.

Run with: .venv/bin/python -m pytest tests/label_qc/ -v

Tests exercise the CLI seam on tiny synthetic VCF fixtures, one per rule,
following the TruthQC calibration contract (see
dev_docs/audit/2026-08-28_adversarial_review.md).
"""

import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CLI = _REPO_ROOT / "bin" / "label_qc.py"
CONFIG = _REPO_ROOT / "bin" / "label_qc_config.json"

BASE_HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=Somatic,Description="Somatic variant">
##FILTER=<ID=Germline,Description="Germline variant">
##FILTER=<ID=Reference,Description="Reference/wildtype">
##FILTER=<ID=Artifact,Description="Artifact/technical error">
##FILTER=<ID=NoConsensus,Description="No consensus">
##FILTER=<ID=RNAedit,Description="RNA editing variant">
##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="Alternate allele frequency">
##INFO=<ID=COSMIC_ID,Number=1,Type=String,Description="COSMIC identifier">
##INFO=<ID=COSMIC_CNT,Number=1,Type=Integer,Description="COSMIC count">
##INFO=<ID=FILTERS_CATEGORY,Number=.,Type=String,Description="Filter categories per caller">
##INFO=<ID=UNIFIED_FILTER,Number=1,Type=String,Description="Unified filter status">
##INFO=<ID=N_SUPPORT_CALLERS,Number=1,Type=Integer,Description="Supporting callers">
##INFO=<ID=CALLERS,Number=.,Type=String,Description="All callers">
##INFO=<ID=VAF_MEAN,Number=1,Type=Float,Description="Mean VAF">
##INFO=<ID=DP_MEAN,Number=1,Type=Float,Description="Mean depth">
##INFO=<ID=DP_DNA_MEAN,Number=1,Type=Float,Description="Mean DNA depth">
##INFO=<ID=REDI_ACCESSION,Number=1,Type=String,Description="REDIportal accession">
##INFO=<ID=REDI_EVIDENCE,Number=1,Type=String,Description="RNA editing evidence level">
##INFO=<ID=REDI_CANONICAL,Number=1,Type=String,Description="Canonical A>G or T>C">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""


def rec(chrom, pos, ref, alt, filt, info=None):
    info_s = ";".join(f"{k}={v}" for k, v in (info or {}).items()) or "."
    return "\t".join([chrom, str(pos), ".", ref, alt, ".", filt, info_s])


def write_vcf(path, records):
    path = Path(path)
    path.write_text(BASE_HEADER + "\n".join(records) + "\n")
    return path


def write_manifest(path, rows):
    """rows: list of (sample_id, vcf_path). TruthQC-style header."""
    lines = ["#sample_id\ttruth_vcf"]
    lines += [f"{sid}\t{vcf}" for sid, vcf in rows]
    Path(path).write_text("\n".join(lines) + "\n")
    return path


def run_cli(*argv, expect_rc=0):
    cmd = [sys.executable, str(CLI)] + [str(a) for a in argv]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == expect_rc, (
        f"rc={proc.returncode}\nstdout:{proc.stdout}\nstderr:{proc.stderr}"
    )
    return proc


def read_flagged(outdir):
    path = Path(outdir) / "flagged_sites.tsv.gz"
    assert path.is_file(), "flagged_sites.tsv.gz missing"
    with gzip.open(path, "rt") as fh:
        lines = [ln.rstrip("\n").split("\t") for ln in fh]
    header, rows = lines[0], lines[1:]
    return [dict(zip(header, r)) for r in rows if r and r[0]]


def read_verdicts(outdir):
    path = Path(outdir) / "samples_qc.tsv"
    assert path.is_file(), "samples_qc.tsv missing"
    lines = path.read_text().strip().split("\n")
    header = lines[0].split("\t")
    return {r["sample"]: r for r in
            (dict(zip(header, ln.split("\t"))) for ln in lines[1:])}


def read_summary(outdir):
    return json.loads((Path(outdir) / "summary.json").read_text())


def somatic(pos, info=None, chrom="chr1", ref="A", alt="G"):
    base = {"FILTERS_CATEGORY": "DNA_strelka:Somatic|DNA_mutect2:Somatic",
            "UNIFIED_FILTER": "Somatic", "N_SUPPORT_CALLERS": "2",
            "CALLERS": "DNA_strelka|DNA_mutect2",
            "VAF_MEAN": "0.25", "DP_MEAN": "40", "DP_DNA_MEAN": "38"}
    base.update(info or {})
    return rec(chrom, pos, ref, alt, "Somatic", base)


# ---------------------------------------------------------------- R1

class TestR1GnomadAF:
    def test_common_af_flagged(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"GNOMAD_AF": "0.005"}),  # >= af_common, < af_strong
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        rows = read_flagged(out)
        assert len(rows) == 1
        assert "R1" in rows[0]["rules"].split("+")
        assert rows[0]["confidence"] == "low"
        assert rows[0]["action"] == "REPORT_ONLY"

    def test_strong_af_is_mid(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"GNOMAD_AF": "0.05"}),  # >= af_strong
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        rows = read_flagged(out)
        assert rows[0]["confidence"] == "mid"
        assert rows[0]["action"].startswith("RELABEL")

    def test_rare_af_not_flagged(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"GNOMAD_AF": "0.0001"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        assert read_flagged(out) == []

    def test_missing_af_not_flagged(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [somatic(1000)])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        assert read_flagged(out) == []


# ---------------------------------------------------------------- R2

class TestR2NoDnaSomatic:
    def test_rna_only_category_flagged(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"FILTERS_CATEGORY": "RNA_strelka:Somatic|RNA_mutect2:Somatic"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        rows = read_flagged(out)
        assert "R2" in rows[0]["rules"].split("+")

    def test_dna_somatic_present_not_flagged(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"FILTERS_CATEGORY": "DNA_strelka:Somatic|RNA_mutect2:Artifact"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        assert read_flagged(out) == []

    def test_dna_nonsomatic_entry_flagged(self, tmp_path):
        # DNA caller saw the site but judged it Reference — still zero DNA Somatic
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"FILTERS_CATEGORY": "DNA_deepsomatic:Reference|RNA_strelka:Somatic"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        rows = read_flagged(out)
        assert "R2" in rows[0]["rules"].split("+")


# ---------------------------------------------------------------- R3

class TestR3CosmicCommon:
    def test_cosmic_hit_with_common_af_is_high(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"GNOMAD_AF": "0.2", "COSMIC_ID": "COSV58987646",
                           "COSMIC_CNT": "5"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        rows = read_flagged(out)
        rules = rows[0]["rules"].split("+")
        assert "R1" in rules and "R3" in rules
        assert rows[0]["confidence"] == "high"
        assert rows[0]["action"] == "DROP"

    def test_cosmic_hit_rare_af_not_flagged(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"GNOMAD_AF": "0.00001", "COSMIC_ID": "COSV58987646"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        assert read_flagged(out) == []


# ---------------------------------------------------------------- R4

class TestR4ColocatedConflict:
    def test_colocated_germline_record_is_high(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000),
            rec("chr1", 1000, "A", "G", "Germline", {}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        rows = read_flagged(out)
        assert len(rows) == 1
        assert "R4" in rows[0]["rules"].split("+")
        assert rows[0]["confidence"] == "high"

    def test_colocated_reference_record_is_high(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000),
            rec("chr1", 1000, "A", "G", "Reference", {}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        assert "R4" in read_flagged(out)[0]["rules"].split("+")

    def test_different_allele_no_conflict(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000),
            rec("chr1", 1000, "A", "T", "Germline", {}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        assert read_flagged(out) == []


# ---------------------------------------------------------------- R6

class TestR6SelfContradiction:
    def test_unified_reference_is_high(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"UNIFIED_FILTER": "Reference"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        rows = read_flagged(out)
        assert "R6" in rows[0]["rules"].split("+")
        assert rows[0]["confidence"] == "high"
        assert rows[0]["action"] == "DROP"

    def test_unified_germline_is_high(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"UNIFIED_FILTER": "Germline"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        assert "R6" in read_flagged(out)[0]["rules"].split("+")

    def test_unified_somatic_not_flagged(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [somatic(1000)])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        assert read_flagged(out) == []


# ---------------------------------------------------------------- R7 (RNA editing)

class TestR7RnaEditing:
    def test_rediportal_high_evidence_flagged(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"REDI_EVIDENCE": "HIGH", "REDI_ACCESSION": "REDI123",
                           "REDI_CANONICAL": "YES"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        rows = read_flagged(out)
        assert "R7" in rows[0]["rules"].split("+")

    def test_no_rediportal_not_flagged(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [somatic(1000)])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        assert read_flagged(out) == []


# ---------------------------------------------------------------- R8 (clustered)

class TestR8Clustered:
    def test_three_somatic_in_window_flagged(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000), somatic(1010), somatic(1020),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        rows = read_flagged(out)
        assert len(rows) == 3
        assert all("R8" in r["rules"].split("+") for r in rows)
        assert all(r["confidence"] == "low" for r in rows)

    def test_spread_out_not_flagged(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000), somatic(5000), somatic(9000),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        assert read_flagged(out) == []


# ---------------------------------------------------------------- R9 (caller agreement)

class TestR9CallerAgreement:
    def test_single_caller_support_flagged(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"N_SUPPORT_CALLERS": "1",
                           "FILTERS_CATEGORY": "DNA_strelka:Somatic",
                           "CALLERS": "DNA_strelka"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        rows = read_flagged(out)
        assert "R9" in rows[0]["rules"].split("+")
        assert rows[0]["confidence"] == "low"


# ---------------------------------------------------------------- disposition tiers

class TestDisposition:
    def test_strong_pair_is_high_drop(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"GNOMAD_AF": "0.005",
                           "FILTERS_CATEGORY": "RNA_strelka:Somatic"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        rows = read_flagged(out)
        assert rows[0]["confidence"] == "high"
        assert rows[0]["action"] == "DROP"

    def test_only_target_filter_records_evaluated(self, tmp_path):
        # Germline records with scary INFO must not be flagged
        vcf = write_vcf(tmp_path / "s.vcf", [
            rec("chr1", 1000, "A", "G", "Germline",
                {"GNOMAD_AF": "0.5", "UNIFIED_FILTER": "Reference"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        assert read_flagged(out) == []


# ---------------------------------------------------------------- rule toggling / config

class TestConfig:
    def test_rule_can_be_disabled(self, tmp_path):
        cfg = json.loads(CONFIG.read_text())
        cfg["variant_rules"]["R6"]["enabled"] = False
        cfg_path = tmp_path / "cfg.json"
        cfg_path.write_text(json.dumps(cfg))
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"UNIFIED_FILTER": "Reference"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out,
                "--config", cfg_path)
        assert read_flagged(out) == []

    def test_threshold_override(self, tmp_path):
        cfg = json.loads(CONFIG.read_text())
        cfg["variant_rules"]["R1"]["af_common"] = 0.5
        cfg_path = tmp_path / "cfg.json"
        cfg_path.write_text(json.dumps(cfg))
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"GNOMAD_AF": "0.1"}),  # below overridden threshold
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out,
                "--config", cfg_path)
        assert read_flagged(out) == []


# ---------------------------------------------------------------- sample gates

def clean_records(n, start=100000):
    return [somatic(start + i * 500) for i in range(n)]


def make_cohort(tmp_path, specs):
    """specs: list of (sid, records). Returns manifest path."""
    rows = []
    for sid, records in specs:
        vcf = write_vcf(tmp_path / f"{sid}.vcf", records)
        rows.append((sid, vcf))
    return write_manifest(tmp_path / "manifest.tsv", rows)


class TestSampleGates:
    def test_s2_contradiction_rate_fail(self, tmp_path):
        bad = [somatic(1000 + i * 100, {"UNIFIED_FILTER": "Reference"})
               for i in range(30)]
        good = clean_records(30)
        manifest = make_cohort(tmp_path, [("bad", bad), ("good", good)])
        out = tmp_path / "out"
        run_cli("--samples", manifest, "--out", out)
        v = read_verdicts(out)
        assert v["bad"]["verdict"] == "FAIL"
        assert "S2" in v["bad"]["gates_failed"]
        assert v["good"]["verdict"] == "PASS"

    def test_s3_rna_only_fraction_fail(self, tmp_path):
        # RNA-only Somatic calls at common population frequency = germline leakage
        rna_only = [somatic(1000 + i * 100,
                            {"FILTERS_CATEGORY": "RNA_strelka:Somatic|RNA_mutect2:Somatic",
                             "CALLERS": "RNA_strelka|RNA_mutect2",
                             "GNOMAD_AF": "0.2"})
                    for i in range(30)]
        good = clean_records(30)
        manifest = make_cohort(tmp_path, [("rnaonly", rna_only), ("good", good)])
        out = tmp_path / "out"
        run_cli("--samples", manifest, "--out", out)
        v = read_verdicts(out)
        assert v["rnaonly"]["verdict"] == "FAIL"
        assert "S3" in v["rnaonly"]["gates_failed"]
        assert v["good"]["verdict"] == "PASS"

    def test_s3_rna_only_rare_af_not_failed(self, tmp_path):
        # RNA-only support without common AF is the norm for rescued records
        rna_only = [somatic(1000 + i * 100,
                            {"FILTERS_CATEGORY": "RNA_strelka:Somatic|RNA_mutect2:Somatic",
                             "CALLERS_SUPPORT": "RNA_strelka|RNA_mutect2"})
                    for i in range(30)]
        good = clean_records(30)
        manifest = make_cohort(tmp_path, [("rnaonly", rna_only), ("good", good)])
        out = tmp_path / "out"
        run_cli("--samples", manifest, "--out", out)
        v = read_verdicts(out)
        assert "S3" not in v["rnaonly"]["gates_failed"].split(";")

    def test_s5_no_dna_participation_fail(self, tmp_path):
        rna_only = [somatic(1000 + i * 100,
                            {"FILTERS_CATEGORY": "RNA_strelka:Somatic|RNA_mutect2:Somatic",
                             "CALLERS": "RNA_strelka|RNA_mutect2"})
                    for i in range(10)]
        manifest = make_cohort(tmp_path, [("nodna", rna_only)])
        out = tmp_path / "out"
        run_cli("--samples", manifest, "--out", out)
        v = read_verdicts(out)
        assert v["nodna"]["verdict"] == "FAIL"
        assert "S5" in v["nodna"]["gates_failed"]

    def test_s1_count_outlier_fail(self, tmp_path):
        specs = [(f"n{i}", clean_records(100)) for i in range(7)]
        specs.append(("huge", clean_records(20000)))
        manifest = make_cohort(tmp_path, specs)
        out = tmp_path / "out"
        run_cli("--samples", manifest, "--out", out)
        v = read_verdicts(out)
        assert v["huge"]["verdict"] == "FAIL"
        assert "S1" in v["huge"]["gates_failed"]
        assert v["n0"]["verdict"] == "PASS"

    def test_s6_indel_fraction_warn(self, tmp_path):
        indels = [somatic(1000 + i * 100, chrom="chr1", ref="AT", alt="A")
                  for i in range(30)]
        manifest = make_cohort(tmp_path, [("indels", indels)])
        out = tmp_path / "out"
        run_cli("--samples", manifest, "--out", out)
        v = read_verdicts(out)
        assert v["indels"]["verdict"] in ("WARN", "FAIL")
        assert "S6" in v["indels"]["gates_warned"] + v["indels"]["gates_failed"]

    def test_s6_het_vaf_fraction_warn(self, tmp_path):
        het = [somatic(1000 + i * 100, {"VAF_MEAN": "0.5"}) for i in range(30)]
        manifest = make_cohort(tmp_path, [("het", het)])
        out = tmp_path / "out"
        run_cli("--samples", manifest, "--out", out)
        v = read_verdicts(out)
        assert "S6" in v["het"]["gates_warned"] + v["het"]["gates_failed"]

    def test_s7_low_coverage_warn(self, tmp_path):
        low = [somatic(1000 + i * 100, {"DP_MEAN": "3", "DP_DNA_MEAN": "3"})
               for i in range(30)]
        manifest = make_cohort(tmp_path, [("lowdp", low)])
        out = tmp_path / "out"
        run_cli("--samples", manifest, "--out", out)
        v = read_verdicts(out)
        assert "S7" in v["lowdp"]["gates_warned"] + v["lowdp"]["gates_failed"]


# ---------------------------------------------------------------- output contract

class TestOutputContract:
    def test_four_part_contract(self, tmp_path):
        specs = [("a", clean_records(20)), ("b", clean_records(20))]
        manifest = make_cohort(tmp_path, specs)
        out = tmp_path / "out"
        run_cli("--samples", manifest, "--out", out)
        assert (out / "report.md").is_file()
        assert (out / "summary.json").is_file()
        assert (out / "samples_qc.tsv").is_file()
        assert (out / "flagged_sites.tsv.gz").is_file()
        # no cleaned VCFs without --apply
        assert not (out / "cleaned_vcf").exists()
        summary = read_summary(out)
        assert summary["n_samples"] == 2
        assert {s["sample"] for s in summary["samples"]} == {"a", "b"}
        report = (out / "report.md").read_text()
        assert "PASS" in report

    def test_missing_input_is_error(self, tmp_path):
        proc = subprocess.run(
            [sys.executable, str(CLI), "--truth-vcf",
             str(tmp_path / "nope.vcf"), "--out", str(tmp_path / "o")],
            capture_output=True, text=True)
        assert proc.returncode != 0


# ---------------------------------------------------------------- apply mode

class TestApply:
    def test_apply_drops_high_relabels_mid_inputs_untouched(self, tmp_path):
        records = [
            somatic(1000, {"UNIFIED_FILTER": "Reference"}),          # HIGH -> DROP
            somatic(2000, {"GNOMAD_AF": "0.05"}),                    # MID  -> RELABEL
            somatic(3000, {"GNOMAD_AF": "0.002"}),                   # LOW  -> keep as-is
            somatic(4000),                                           # clean
        ]
        vcf = write_vcf(tmp_path / "s.vcf", records)
        digest_before = hashlib.sha256(vcf.read_bytes()).hexdigest()
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out, "--apply")

        # input untouched
        assert hashlib.sha256(vcf.read_bytes()).hexdigest() == digest_before

        cleaned = out / "cleaned_vcf" / "s1.vcf.gz"
        assert cleaned.is_file()
        with gzip.open(cleaned, "rt") as fh:
            body = [ln for ln in fh if not ln.startswith("#")]
        positions = [int(ln.split("\t")[1]) for ln in body]
        assert 1000 not in positions            # dropped
        assert 2000 in positions and 3000 in positions and 4000 in positions
        by_pos = {int(ln.split("\t")[1]): ln for ln in body}
        assert by_pos[2000].split("\t")[6] == "Germline"   # relabeled
        assert "LABEL_QC" in by_pos[2000].split("\t")[7]   # provenance tag
        assert by_pos[3000].split("\t")[6] == "Somatic"    # report-only: kept verbatim
        assert (out / "samples_cleaned.tsv").is_file()

    def test_no_apply_means_no_cleaned(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", [
            somatic(1000, {"UNIFIED_FILTER": "Reference"}),
        ])
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "s1", "--out", out)
        assert not (out / "cleaned_vcf").exists()


# ---------------------------------------------------------------- single-sample mode

class TestSingleSampleMode:
    def test_single_sample_pass(self, tmp_path):
        vcf = write_vcf(tmp_path / "s.vcf", clean_records(50))
        out = tmp_path / "out"
        run_cli("--truth-vcf", vcf, "--sample-id", "solo", "--out", out)
        v = read_verdicts(out)
        assert v["solo"]["verdict"] == "PASS"
