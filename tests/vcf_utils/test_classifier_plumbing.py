"""Regression tests for classifier plumbing (audit M4/M9 + tie fallback + min-alt floor, ticket 05).

Covers:
- M9: caller-rejected (Artifact) records must not count toward caller support
  (``support_callers`` / ``passes_consensus`` / ``N_SUPPORT_CALLERS`` /
  ``PASSES_CONSENSUS``).
- M4: CLI ``--snv_thr/--indel_thr`` reach the rescue classifier instead of
  being silently dropped to DEFAULT_THRESHOLDS.
- Tie fallback: UNIFIED_FILTER_DNA/RNA majority ties resolve to Artifact per
  docs/consensus_vcf_rules.md (no Somatic-favoring tiebreak).
- Min-alt floor: a configurable minimum tumor alt-read floor
  (``consensus_min_alt_support``, default 3) blocks 1-2-alt-read caller PASSes
  from counting as full support votes.

Run with: .venv/bin/python -m pytest tests/vcf_utils/test_classifier_plumbing.py -v
"""

import subprocess
import sys
from pathlib import Path

import pytest
from cyvcf2 import VCF

# Ensure bin/ is on the Python path
_BIN_DIR = Path(__file__).resolve().parent.parent.parent / "bin"
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

from vcf_utils import classification
from vcf_utils.aggregation import (
    aggregate_genotypes,
    aggregate_variants,
    read_variants_from_vcf,
)
from vcf_utils.io_utils import write_union_vcf

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"

# ---------------------------------------------------------------------------
# Synthetic caller VCF fixtures
# ---------------------------------------------------------------------------

# Mutect2 (normal CRAM first -> sample 0 = normal). FILTER=weak_evidence ->
# classified Artifact by classify_mutect2_variant.
MUTECT2_REJECTED_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=weak_evidence,Description="Mutation does not meet likelihood threshold">
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOO8801DN\tCOO8801DT
chr1\t1000\t.\tA\tG\t.\tweak_evidence\t.\tGT:AD:DP:AF\t0/0:80,0:80:0.0\t0/1:60,40:100:0.4
"""

# Mutect2 PASS with good tumor alt support (40).
MUTECT2_PASS_VCF = MUTECT2_REJECTED_VCF.replace("weak_evidence\t", "PASS\t")

# Mutect2 PASS with only 1 tumor alt read.
MUTECT2_LOWALT_VCF = MUTECT2_REJECTED_VCF.replace(
    "weak_evidence\t.\tGT:AD:DP:AF\t0/0:80,0:80:0.0\t0/1:60,40:100:0.4",
    "PASS\t.\tGT:AD:DP:AF\t0/0:80,0:80:0.0\t0/1:99,1:100:0.01",
)

# Strelka (fixed NORMAL,TUMOR order). FILTER=LowDepth with NT=conflict ->
# classified Artifact by classify_strelka_variant.
STRELKA_REJECTED_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowDepth,Description="Low depth">
##contig=<ID=chr1,length=1000000>
##INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal sample">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AU,Number=2,Type=Integer,Description="A allele counts (tier1,tier2)">
##FORMAT=<ID=CU,Number=2,Type=Integer,Description="C allele counts (tier1,tier2)">
##FORMAT=<ID=GU,Number=2,Type=Integer,Description="G allele counts (tier1,tier2)">
##FORMAT=<ID=TU,Number=2,Type=Integer,Description="T allele counts (tier1,tier2)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR
chr1\t1000\t.\tA\tG\t.\tLowDepth\tNT=conflict\tGT:DP:AU:CU:GU:TU\t0/0:50:50,0:0,0:0,0:0,0\t0/1:100:60,0:0,0:40,0:0,0
"""

# Strelka PASS with good tumor alt support (G tier1 = 40).
STRELKA_PASS_VCF = STRELKA_REJECTED_VCF.replace("LowDepth\tNT=conflict", "PASS\tNT=ref")

# Strelka PASS with only 1 tumor alt read (G tier1 = 1).
STRELKA_LOWALT_VCF = STRELKA_REJECTED_VCF.replace(
    "LowDepth\tNT=conflict", "PASS\tNT=ref"
).replace("0/1:100:60,0:0,0:40,0:0,0", "0/1:100:99,0:0,0:1,0:0,0")

# Minimal "consensus" VCFs for the rescue path (FILTER holds biological class).
DNA_CONSENSUS_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=Somatic,Description="Somatic variant">
##contig=<ID=chr1,length=1000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1000\t.\tA\tG\t.\tSomatic\t.
"""

RNA_CONSENSUS_VCF = DNA_CONSENSUS_VCF


def _write_vcf(directory, name, content):
    path = directory / name
    path.write_text(content)
    return str(path)


def _read_pair(tmp_path, mutect2_content, strelka_content):
    m = read_variants_from_vcf(
        _write_vcf(tmp_path, "sample.mutect2.variants.vcf", mutect2_content), "mutect2"
    )
    s = read_variants_from_vcf(
        _write_vcf(tmp_path, "sample.strelka.variants.vcf", strelka_content), "strelka"
    )
    return [("mutect2", m, None), ("strelka", s, None)]


def _run_consensus(tmp_path, mutect2_content, strelka_content, extra_args=None):
    input_dir = tmp_path / "callers"
    input_dir.mkdir(exist_ok=True)
    _write_vcf(input_dir, "sample.mutect2.variants.vcf", mutect2_content)
    _write_vcf(input_dir, "sample.strelka.variants.vcf", strelka_content)
    out_prefix = tmp_path / "out.consensus"
    cmd = [
        sys.executable,
        str(BIN_DIR / "run_consensus_vcf.py"),
        "--input_dir",
        str(input_dir),
        "--out_prefix",
        str(out_prefix),
        "--output_format",
        "vcf",
        "--expected_callers",
        "mutect2,strelka",
    ] + (extra_args or [])
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return str(out_prefix) + ".vcf"


def _info_scalar(variant, key):
    val = variant.INFO.get(key)
    if isinstance(val, (tuple, list)):
        val = val[0] if val else None
    return val


# ---------------------------------------------------------------------------
# M9: caller-rejected records must not count as support
# ---------------------------------------------------------------------------


class TestCallerRejectedSupport:
    def test_both_callers_rejected_give_zero_support(self, tmp_path):
        """Two callers' Artifact records at the same site: support count 0,
        no consensus. Pre-fix both counted -> n_support 2 -> PASSES_CONSENSUS."""
        collections = _read_pair(tmp_path, MUTECT2_REJECTED_VCF, STRELKA_REJECTED_VCF)
        aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)

        assert len(aggregated) == 1
        data = next(iter(aggregated.values()))
        # Detection is preserved; only support votes change
        assert sorted(data["callers"]) == ["mutect2", "strelka"]
        assert data["support_callers"] == set()
        assert data["passes_consensus"] is False

    def test_e2e_rejected_site_has_no_support_and_artifact_filter(self, tmp_path):
        """End-to-end: PASSES_CONSENSUS=NO, N_SUPPORT_CALLERS=0, and the
        majority vote still labels the site Artifact."""
        out = _run_consensus(tmp_path, MUTECT2_REJECTED_VCF, STRELKA_REJECTED_VCF)
        (rec,) = list(VCF(out))
        assert _info_scalar(rec, "PASSES_CONSENSUS") == "NO"
        assert int(_info_scalar(rec, "N_SUPPORT_CALLERS")) == 0
        assert rec.FILTER == "Artifact"

    def test_one_rejected_one_passing_counts_one_support(self, tmp_path):
        collections = _read_pair(tmp_path, MUTECT2_REJECTED_VCF, STRELKA_PASS_VCF)
        aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)
        data = next(iter(aggregated.values()))
        assert data["support_callers"] == {"strelka"}
        assert data["passes_consensus"] is False


# ---------------------------------------------------------------------------
# Min-alt floor: configurable minimum tumor alt reads for a support vote
# ---------------------------------------------------------------------------


class TestMinAltFloor:
    def test_one_alt_read_votes_blocked_by_default_floor(self, tmp_path):
        """Both callers PASS but with a single tumor alt read each: with the
        default floor (3) neither counts as support. Pre-fix both counted."""
        collections = _read_pair(tmp_path, MUTECT2_LOWALT_VCF, STRELKA_LOWALT_VCF)
        aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)
        data = next(iter(aggregated.values()))
        assert data["support_callers"] == set()
        assert data["passes_consensus"] is False

    def test_floor_is_configurable(self, tmp_path):
        collections = _read_pair(tmp_path, MUTECT2_LOWALT_VCF, STRELKA_LOWALT_VCF)
        aggregated = aggregate_variants(
            collections, snv_threshold=2, indel_threshold=2, min_alt_support=1
        )
        data = next(iter(aggregated.values()))
        assert data["support_callers"] == {"mutect2", "strelka"}
        assert data["passes_consensus"] is True

    def test_floor_disabled_with_zero(self, tmp_path):
        collections = _read_pair(tmp_path, MUTECT2_LOWALT_VCF, STRELKA_LOWALT_VCF)
        aggregated = aggregate_variants(
            collections, snv_threshold=2, indel_threshold=2, min_alt_support=0
        )
        data = next(iter(aggregated.values()))
        assert data["passes_consensus"] is True

    def test_records_without_alt_evidence_keep_legacy_behavior(self):
        """No AD evidence (e.g. consensus callers in rescue mode): the floor
        cannot be evaluated, so the record counts as before."""
        record = {
            "CHROM": "chr1",
            "POS": 1000,
            "REF": "A",
            "ALT": "G",
            "is_snv": True,
            "filter_original": "Somatic",
            "filter_normalized": "Somatic",
            "filter_category": "Somatic",
            "quality": None,
            "genotype": {"GT": None, "DP": None, "AD": None, "VAF": None, "GQ": None},
            "id": None,
            "classification": "Somatic",
        }
        vkey = "1:1000:A:G"
        collections = [
            ("DNA_consensus", {vkey: dict(record)}, "DNA"),
            ("RNA_consensus", {vkey: dict(record)}, "RNA"),
        ]
        aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)
        data = next(iter(aggregated.values()))
        assert data["support_callers"] == {"DNA_consensus", "RNA_consensus"}
        assert data["passes_consensus"] is True


# ---------------------------------------------------------------------------
# M4: CLI thresholds must reach the rescue classifier
# ---------------------------------------------------------------------------


class TestRescueThresholdWiring:
    def test_thresholds_reach_rescue_classifier(self, monkeypatch):
        """compute_unified_classification_rescue must build a classifier
        configured with the given CLI thresholds (mirroring the consensus
        path) instead of silently using DEFAULT_THRESHOLDS."""
        captured = {}

        class FakeClassifier:
            def __init__(self, config=None):
                captured["config"] = config

            def classify_rescue_variant(self, variant_data, modality_map):
                return "Somatic"

        monkeypatch.setattr(
            "vcf_utils.variant_classifier_unified.UnifiedVariantClassifier",
            FakeClassifier,
        )

        result = classification.compute_unified_classification_rescue(
            {"callers": []}, {}, snv_threshold=3, indel_threshold=4
        )
        assert result == "Somatic"
        assert captured["config"] == {
            "consensus_snv_threshold": 3,
            "consensus_indel_threshold": 4,
        }

    def test_default_path_uses_real_classifier(self):
        """Without explicit thresholds the rescue classification is unchanged."""
        variant_data = {
            "callers": ["DNA_consensus"],
            "filters_normalized": ["Somatic"],
        }
        assert (
            classification.compute_unified_classification_rescue(variant_data, {})
            == "Somatic"
        )

    def test_write_union_vcf_forwards_thresholds(self, monkeypatch, tmp_path):
        """write_union_vcf must forward its snv_threshold/indel_threshold to
        the rescue classification function."""
        captured = {}
        real = classification.compute_unified_classification_rescue

        def spy(variant_data, modality_map, **kwargs):
            captured.update(kwargs)
            return real(variant_data, modality_map, **kwargs)

        monkeypatch.setattr(
            classification, "compute_unified_classification_rescue", spy
        )

        data = _minimal_aggregated_record(
            callers=["DNA_mutect2", "DNA_strelka"],
            filters=["Somatic", "Somatic"],
        )
        modality_map = {"DNA_mutect2": "DNA", "DNA_strelka": "DNA"}
        template = _template_header(tmp_path)
        out = tmp_path / "rescue.vcf"
        write_union_vcf(
            {"1:1000:A:G": data},
            template,
            "SAMPLE",
            str(out),
            "vcf",
            ["DNA_mutect2", "DNA_strelka"],
            modality_map=modality_map,
            snv_threshold=5,
            indel_threshold=3,
        )
        assert captured.get("snv_threshold") == 5
        assert captured.get("indel_threshold") == 3

    def test_e2e_rescue_cli_thresholds_change_support(self, tmp_path):
        """CLI --snv_thr changes rescue-run behavior (PASSES_CONSENSUS)."""
        dna = _write_vcf(tmp_path, "sample.dna.consensus.vcf", DNA_CONSENSUS_VCF)
        rna = _write_vcf(tmp_path, "sample.rna.consensus.vcf", RNA_CONSENSUS_VCF)

        def run(thr):
            out_prefix = tmp_path / f"rescued.thr{thr}"
            cmd = [
                sys.executable,
                str(BIN_DIR / "run_rescue_vcf.py"),
                "--dna_consensus",
                dna,
                "--rna_consensus",
                rna,
                "--out_prefix",
                str(out_prefix),
                "--output_format",
                "vcf",
                "--snv_thr",
                str(thr),
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, result.stderr
            (rec,) = list(VCF(str(out_prefix) + ".vcf"))
            return _info_scalar(rec, "PASSES_CONSENSUS")

        assert run(2) == "YES"  # DNA_consensus + RNA_consensus = 2 votes
        assert run(3) == "NO"


# ---------------------------------------------------------------------------
# Tie fallback: UNIFIED_FILTER_DNA/RNA ties resolve to Artifact
# ---------------------------------------------------------------------------


def _template_header(tmp_path):
    template_path = _write_vcf(tmp_path, "template.vcf", DNA_CONSENSUS_VCF)
    return VCF(template_path)


def _minimal_aggregated_record(callers, filters, chrom="chr1", pos=1000):
    return {
        "CHROM": chrom,
        "POS": pos,
        "REF": "A",
        "ALT": "G",
        "is_snv": True,
        "callers": callers,
        "filters_original": list(filters),
        "filters_normalized": list(filters),
        "filters_category": list(filters),
        "qualities": [],
        "genotypes": {},
        "ids": [],
        "support_callers": set(callers),
        "passes_consensus": True,
        "gt_aggregated": aggregate_genotypes({}, callers),
    }


class TestTieFallback:
    def test_dna_tie_resolves_to_artifact(self, tmp_path):
        """DNA callers split 1:1 Somatic/Germline -> UNIFIED_FILTER_DNA must be
        Artifact (docs/consensus_vcf_rules.md), not Somatic via priority."""
        data = _minimal_aggregated_record(
            callers=["DNA_mutect2", "DNA_strelka"],
            filters=["Somatic", "Germline"],
        )
        modality_map = {"DNA_mutect2": "DNA", "DNA_strelka": "DNA"}
        out = tmp_path / "tie.vcf"
        write_union_vcf(
            {"1:1000:A:G": data},
            _template_header(tmp_path),
            "SAMPLE",
            str(out),
            "vcf",
            ["DNA_mutect2", "DNA_strelka"],
            modality_map=modality_map,
        )
        (rec,) = list(VCF(str(out)))
        assert _info_scalar(rec, "UNIFIED_FILTER_DNA") == "Artifact"

    def test_rna_tie_resolves_to_artifact(self, tmp_path):
        data = _minimal_aggregated_record(
            callers=["RNA_mutect2", "RNA_strelka"],
            filters=["Somatic", "Reference"],
        )
        modality_map = {"RNA_mutect2": "RNA", "RNA_strelka": "RNA"}
        out = tmp_path / "tie_rna.vcf"
        write_union_vcf(
            {"1:1000:A:G": data},
            _template_header(tmp_path),
            "SAMPLE",
            str(out),
            "vcf",
            ["RNA_mutect2", "RNA_strelka"],
            modality_map=modality_map,
        )
        (rec,) = list(VCF(str(out)))
        assert _info_scalar(rec, "UNIFIED_FILTER_RNA") == "Artifact"

    def test_clear_majority_unchanged(self, tmp_path):
        data = _minimal_aggregated_record(
            callers=["DNA_mutect2", "DNA_strelka", "DNA_deepsomatic"],
            filters=["Somatic", "Somatic", "Germline"],
        )
        modality_map = {c: "DNA" for c in data["callers"]}
        out = tmp_path / "majority.vcf"
        write_union_vcf(
            {"1:1000:A:G": data},
            _template_header(tmp_path),
            "SAMPLE",
            str(out),
            "vcf",
            list(data["callers"]),
            modality_map=modality_map,
        )
        (rec,) = list(VCF(str(out)))
        assert _info_scalar(rec, "UNIFIED_FILTER_DNA") == "Somatic"
