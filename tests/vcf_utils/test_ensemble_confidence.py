"""Tests for descriptive caller support and complete caller-panel validation.

ENS_SUPPORT is k supporting callers out of the explicit expected panel n.
Missing, duplicate, and unexpected caller files are errors. The removed Wilson
fields are not label-confidence probabilities and must not reappear.

Run with: .venv/bin/python -m pytest tests/vcf_utils/test_ensemble_confidence.py -v
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

from vcf_utils.io_utils import create_output_header

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"


# --- Writer integration ---------------------------------------------------

# Variant at chr1:1000 is called by mutect2 and strelka only (deepsomatic
# absent at the site) -> ENS_SUPPORT=2/3. Variant at chr1:2000 is called by
# all three -> 3/3. Tumor alt reads are 40, well above the default
# min_alt_support floor of 3.
MUTECT2_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR
chr1\t1000\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:AF\t0/0:80,0:80:0.0\t0/1:60,40:100:0.4
chr1\t2000\t.\tC\tT\t.\tPASS\t.\tGT:AD:DP:AF\t0/0:80,0:80:0.0\t0/1:60,40:100:0.4
"""

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
chr1\t2000\t.\tC\tT\t.\tPASS\tNT=ref\tGT:DP:AU:CU:GU:TU\t0/0:80:80,0:0:0,0:0,0\t0/1:100:60,0:0,0:0,0:40,0
"""

# DeepSomatic: no record at chr1:1000 (absent = non-support vote), supports 2000.
DEEPSOMATIC_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR
chr1\t2000\t.\tC\tT\t.\tPASS\t.\tGT:AD:DP\t0/0:80,0:80\t0/1:60,40:100
"""


def _write_vcf(directory, name, content):
    path = directory / name
    path.write_text(content)
    return str(path)


@pytest.fixture(scope="module")
def consensus_vcf(tmp_path_factory):
    """Consensus VCF from 3 configured callers, one absent at chr1:1000."""
    tmp_path = tmp_path_factory.mktemp("ens_conf")
    input_dir = tmp_path / "callers"
    input_dir.mkdir()
    _write_vcf(input_dir, "sample.mutect2.variants.vcf", MUTECT2_VCF)
    _write_vcf(input_dir, "sample.strelka.variants.vcf", STRELKA_VCF)
    _write_vcf(input_dir, "sample.deepsomatic.variants.vcf", DEEPSOMATIC_VCF)
    out_prefix = tmp_path / "out.consensus"
    cmd = [
        sys.executable,
        str(BIN_DIR / "run_consensus_vcf.py"),
        "--input_dir",
        str(input_dir),
        "--expected_callers",
        "mutect2,strelka,deepsomatic",
        "--out_prefix",
        str(out_prefix),
        "--output_format",
        "vcf",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return str(out_prefix) + ".vcf"


def _records_by_pos(vcf_path):
    return {v.POS: v for v in VCF(vcf_path)}


class TestEnsembleSupportInfo:
    """End-to-end descriptive support over the explicit three-caller panel."""

    def test_support_fraction_counts_absent_caller_as_non_support(self, consensus_vcf):
        records = _records_by_pos(consensus_vcf)
        assert records[1000].INFO.get("ENS_SUPPORT") == "2/3"
        assert records[2000].INFO.get("ENS_SUPPORT") == "3/3"

    def test_header_declares_only_descriptive_support(self, consensus_vcf):
        header = str(VCF(consensus_vcf).raw_header)
        assert "ID=ENS_SUPPORT" in header
        assert "ID=ENS_CONF_LO" not in header
        assert "ID=ENS_CONF_HI" not in header


def _run_panel(tmp_path, files, expected):
    input_dir = tmp_path / "panel"
    input_dir.mkdir()
    for name, content in files:
        _write_vcf(input_dir, name, content)
    out_prefix = tmp_path / "panel.consensus"
    return subprocess.run(
        [
            sys.executable,
            str(BIN_DIR / "run_consensus_vcf.py"),
            "--input_dir", str(input_dir),
            "--expected_callers", expected,
            "--out_prefix", str(out_prefix),
            "--output_format", "vcf",
            "--snv_thr", "2",
            "--indel_thr", "2",
        ],
        capture_output=True,
        text=True,
    )


class TestCallerPanelValidation:
    def test_missing_caller_is_fatal(self, tmp_path):
        result = _run_panel(
            tmp_path,
            [("s.mutect2.variants.vcf", MUTECT2_VCF),
             ("s.strelka.variants.vcf", STRELKA_VCF)],
            "mutect2,strelka,deepsomatic",
        )
        assert result.returncode == 2
        assert "missing=['deepsomatic']" in result.stderr

    def test_unexpected_caller_is_fatal(self, tmp_path):
        result = _run_panel(
            tmp_path,
            [("s.mutect2.variants.vcf", MUTECT2_VCF),
             ("s.strelka.variants.vcf", STRELKA_VCF),
             ("s.deepsomatic.variants.vcf", DEEPSOMATIC_VCF)],
            "mutect2,strelka",
        )
        assert result.returncode == 2
        assert "unexpected=['deepsomatic']" in result.stderr

    def test_duplicate_caller_file_is_fatal(self, tmp_path):
        result = _run_panel(
            tmp_path,
            [("s.mutect2.variants.vcf", MUTECT2_VCF),
             ("other.mutect2.variants.vcf", MUTECT2_VCF),
             ("s.strelka.variants.vcf", STRELKA_VCF)],
            "mutect2,strelka",
        )
        assert result.returncode == 2
        assert "duplicate VCFs discovered for caller mutect2" in result.stderr


class TestRescueHeaderUnchanged:
    """Rescue mode remains free of consensus-only support fields."""

    def _header_info_ids(self, tmp_path, include_rescue_fields):
        template_path = _write_vcf(tmp_path, "template.vcf", MUTECT2_VCF)
        template = VCF(template_path)
        header = create_output_header(
            template, "sample", include_rescue_fields=include_rescue_fields
        )
        return set(header.info.keys())

    def test_consensus_header_has_descriptive_support(self, tmp_path):
        info_ids = self._header_info_ids(tmp_path, include_rescue_fields=False)
        assert "ENS_SUPPORT" in info_ids
        assert not {"ENS_CONF_LO", "ENS_CONF_HI"} & info_ids

    def test_rescue_header_has_no_ens_fields(self, tmp_path):
        info_ids = self._header_info_ids(tmp_path, include_rescue_fields=True)
        assert not {"ENS_SUPPORT", "ENS_CONF_LO", "ENS_CONF_HI"} & info_ids
