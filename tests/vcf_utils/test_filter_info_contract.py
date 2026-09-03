"""Regression tests for the FILTER/INFO contract and header fix (audit M8a, ticket 07).

Bug (a): the consensus output header added a sample taken from the first input
VCF's first sample (the NORMAL for Strelka-ordered VCFs) while records carried
an empty FORMAT/sample column (``FORMAT "." COO8801DN "."``) — a misleading,
dangling sample name.

Fix under test:
- Output header carries NO sample column; records carry no FORMAT/sample data.
- Every emitted record carries CLASSIFICATION_RATIONALE INFO explaining why its
  FILTER was assigned (majority vote result, tie->Artifact, insufficient
  callers, rescue promotion), so FILTER is derivable from the record's own INFO.
- FILTER vocabulary is unchanged (Somatic/Germline/Reference/Artifact/
  NoConsensus/RNAedit).

Run with: .venv/bin/python -m pytest tests/vcf_utils/test_filter_info_contract.py -v
"""

import subprocess
import sys
from pathlib import Path

import pytest
import pysam
from cyvcf2 import VCF

# Ensure bin/ is on the Python path
_BIN_DIR = Path(__file__).resolve().parent.parent.parent / "bin"
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

from vcf_utils.aggregation import aggregate_genotypes
from vcf_utils.io_utils import write_union_vcf

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"

# Mutect2 (normal first, tumor second; tumor alt = 40). Three sites:
#   1000: PASS           -> Somatic vote
#   2000: germline       -> Germline vote
#   3000: PASS           -> Somatic vote (single-caller site)
MUTECT2_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=germline,Description="Evidence indicates this site is germline">
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOO8801DN\tCOO8801DT
chr1\t1000\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:AF:GQ\t0/0:80,0:80:0.0:99\t0/1:60,40:100:0.4:99
chr1\t2000\t.\tC\tT\t.\tgermline\t.\tGT:AD:DP:AF:GQ\t0/1:40,40:80:0.5:99\t0/1:50,40:90:0.44:99
chr1\t3000\t.\tG\tA\t.\tPASS\t.\tGT:AD:DP:AF:GQ\t0/0:80,0:80:0.0:99\t0/1:70,40:110:0.36:99
"""

# Strelka (NORMAL,TUMOR fixed order; tumor alt = 40). Sites 1000 and 2000 only,
# so 3000 is a single-caller (NoConsensus) site.
STRELKA_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000000>
##INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal sample">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AU,Number=2,Type=Integer,Description="A allele counts (tier1,tier2)">
##FORMAT=<ID=CU,Number=2,Type=Integer,Description="C allele counts (tier1,tier2)">
##FORMAT=<ID=GU,Number=2,Type=Integer,Description="G allele counts (tier1,tier2)">
##FORMAT=<ID=TU,Number=2,Type=Integer,Description="T allele counts (tier1,tier2)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR
chr1\t1000\t.\tA\tG\t.\tPASS\tNT=ref\tGT:DP:AU:CU:GU:TU\t0/0:80:80,0:0,0:0,0:0,0\t0/1:100:60,0:0,0:40,0:0,0
chr1\t2000\t.\tC\tT\t.\tPASS\tNT=ref\tGT:DP:AU:CU:GU:TU\t0/0:80:80,0:0,0:0,0:0,0\t0/1:90:50,0:0,0:0,0:40,0
"""

INPUT_SAMPLE_NAMES = {"COO8801DN", "COO8801DT", "NORMAL", "TUMOR"}

BIOLOGICAL_FILTERS = {
    "Somatic",
    "Germline",
    "Reference",
    "Artifact",
    "NoConsensus",
    "RNAedit",
}


def _write_vcf(directory, name, content):
    path = directory / name
    path.write_text(content)
    return str(path)


@pytest.fixture(scope="module")
def consensus_vcf(tmp_path_factory):
    """Real consensus VCF produced from synthetic Mutect2 + Strelka caller VCFs."""
    tmp_path = tmp_path_factory.mktemp("consensus")
    input_dir = tmp_path / "callers"
    input_dir.mkdir()
    _write_vcf(input_dir, "sample.mutect2.variants.vcf", MUTECT2_VCF)
    _write_vcf(input_dir, "sample.strelka.variants.vcf", STRELKA_VCF)
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
        "--snv_thr",
        "2",
        "--indel_thr",
        "2",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return str(out_prefix) + ".vcf"


def _info_scalar(variant, key):
    val = variant.INFO.get(key)
    if isinstance(val, (tuple, list)):
        val = val[0] if val else None
    return val


def _rationale_parts(raw):
    """Parse 'rule:majority|class:Somatic|...' into a dict."""
    parts = {}
    for clause in raw.split("|"):
        key, _, value = clause.partition(":")
        parts[key] = value
    return parts


def _records_by_pos(vcf_path):
    return {v.POS: v for v in VCF(vcf_path)}


class TestHeaderSampleColumn:
    """Ticket 07(a): the output header must not carry a misleading sample column."""

    def test_header_has_no_samples(self, consensus_vcf):
        with pysam.VariantFile(consensus_vcf) as vcf:
            assert list(vcf.header.samples) == []

    def test_column_header_line_has_eight_columns(self, consensus_vcf):
        for line in Path(consensus_vcf).read_text().splitlines():
            if line.startswith("#CHROM"):
                fields = line.split("\t")
                assert fields == [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                ]
                break
        else:
            pytest.fail("no #CHROM header line found")

    def test_records_carry_no_sample_data(self, consensus_vcf):
        for line in Path(consensus_vcf).read_text().splitlines():
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            assert len(fields) == 8, f"record has {len(fields)} columns: {line}"

    def test_no_input_sample_names_leak_into_output(self, consensus_vcf):
        """No input sample name appears as an output sample (column header or
        declared sample). Substring-scanning the whole header is not valid:
        field names like FILTERS_NORMALIZED legitimately contain 'NORMAL'."""
        for line in Path(consensus_vcf).read_text().splitlines():
            if line.startswith("#CHROM"):
                column_samples = set(line.split("\t")[9:])
                assert column_samples.isdisjoint(INPUT_SAMPLE_NAMES)
                break
        else:
            pytest.fail("no #CHROM header line found")
        with pysam.VariantFile(consensus_vcf) as vcf:
            assert set(vcf.header.samples).isdisjoint(INPUT_SAMPLE_NAMES)

    def test_pysam_strict_roundtrip(self, consensus_vcf):
        """A strict VCF reader parses header and records without error, and
        every record INFO key is declared in the header."""
        with pysam.VariantFile(consensus_vcf) as vcf:
            records = list(vcf)
            assert len(records) == 3
            for rec in records:
                for key in rec.info:
                    assert key in vcf.header.info, f"undeclared INFO key {key}"


class TestClassificationRationale:
    """Ticket 07(b,c): every record's FILTER is derivable from its own INFO."""

    def test_header_declares_rationale_field(self, consensus_vcf):
        with pysam.VariantFile(consensus_vcf) as vcf:
            entry = vcf.header.info.get("CLASSIFICATION_RATIONALE")
            assert entry is not None
            assert entry.number == 1
            assert entry.type == "String"
            assert entry.description

    def test_somatic_majority(self, consensus_vcf):
        rec = _records_by_pos(consensus_vcf)[1000]
        assert rec.FILTER == "Somatic"
        rationale = _rationale_parts(_info_scalar(rec, "CLASSIFICATION_RATIONALE"))
        assert rationale["class"] == "Somatic"
        assert rationale["rule"] == "majority"
        assert "Somatic" in rationale["votes"]

    def test_tie_resolves_to_artifact(self, consensus_vcf):
        """Germline (mutect2) vs Somatic (strelka) tie -> Artifact."""
        rec = _records_by_pos(consensus_vcf)[2000]
        assert rec.FILTER == "Artifact"
        rationale = _rationale_parts(_info_scalar(rec, "CLASSIFICATION_RATIONALE"))
        assert rationale["class"] == "Artifact"
        assert rationale["rule"] == "majority_tie"
        assert "Germline" in rationale["votes"]
        assert "Somatic" in rationale["votes"]

    def test_insufficient_callers_noconsensus(self, consensus_vcf):
        rec = _records_by_pos(consensus_vcf)[3000]
        assert rec.FILTER == "NoConsensus"
        rationale = _rationale_parts(_info_scalar(rec, "CLASSIFICATION_RATIONALE"))
        assert rationale["class"] == "NoConsensus"
        assert rationale["rule"] == "insufficient_callers"
        assert rationale["callers"] == "1"
        assert rationale["threshold"] == "2"

    def test_rationale_class_matches_filter_everywhere(self, consensus_vcf):
        for rec in VCF(consensus_vcf):
            raw = _info_scalar(rec, "CLASSIFICATION_RATIONALE")
            assert raw, f"record at {rec.CHROM}:{rec.POS} lacks rationale"
            assert _rationale_parts(raw)["class"] == rec.FILTER

    def test_filter_vocabulary_unchanged(self, consensus_vcf):
        """FILTER keeps the biological-class vocabulary (model label contract)."""
        with pysam.VariantFile(consensus_vcf) as vcf:
            for rec in vcf:
                assert rec.filter.keys()[0] in BIOLOGICAL_FILTERS


class TestRescueRationale:
    """Rescue-mode records also carry a FILTER-consistent rationale."""

    def _write_rescue_record(self, tmp_path, callers, filters):
        data = {
            "CHROM": "chr1",
            "POS": 1000,
            "REF": "A",
            "ALT": "G",
            "is_snv": True,
            "callers": list(callers),
            "filters_original": list(filters),
            "filters_normalized": list(filters),
            "filters_category": list(filters),
            "qualities": [],
            "genotypes": {},
            "ids": [],
            "support_callers": set(callers),
            "passes_consensus": True,
            "rescued": False,
            "gt_aggregated": aggregate_genotypes({}, callers),
        }
        modality_map = {
            c: "DNA" if c.startswith("DNA_") else "RNA"
            for c in callers
            if c.startswith(("DNA_", "RNA_"))
        }
        template = _write_vcf(tmp_path, "template.vcf", MUTECT2_VCF)
        out = tmp_path / "rescue.vcf"
        write_union_vcf(
            {"1:1000:A:G": data},
            VCF(template),
            "SAMPLE",
            str(out),
            "vcf",
            list(callers),
            modality_map=modality_map,
        )
        (rec,) = list(VCF(str(out)))
        return rec

    def test_rescue_promotion_rationale(self, tmp_path):
        """No consensus labels; DNA+RNA individual callers agree on Somatic ->
        rescued as Somatic with a promotion rationale (ticket 06 contract)."""
        rec = self._write_rescue_record(
            tmp_path,
            callers=["DNA_mutect2", "RNA_strelka"],
            filters=["Somatic", "Somatic"],
        )
        assert rec.FILTER == "Somatic"
        rationale = _rationale_parts(_info_scalar(rec, "CLASSIFICATION_RATIONALE"))
        assert rationale["class"] == "Somatic"
        assert rationale["rule"] == "rescue_promotion"

    def test_rescue_rationale_matches_filter(self, tmp_path):
        rec = self._write_rescue_record(
            tmp_path,
            callers=["DNA_consensus", "RNA_consensus"],
            filters=["NoConsensus", "Artifact"],
        )
        rationale = _rationale_parts(_info_scalar(rec, "CLASSIFICATION_RATIONALE"))
        assert rationale["class"] == rec.FILTER

    def test_rescue_output_has_no_sample_column(self, tmp_path):
        self._write_rescue_record(
            tmp_path, callers=["DNA_mutect2"], filters=["Somatic"]
        )
        with pysam.VariantFile(str(tmp_path / "rescue.vcf")) as vcf:
            assert list(vcf.header.samples) == []
