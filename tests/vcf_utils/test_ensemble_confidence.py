"""Tests for per-variant ensemble confidence annotation in consensus VCF output.

Semantics under test: ENS_SUPPORT / ENS_CONF_LO / ENS_CONF_HI carry the 95%
Wilson score confidence interval on the caller support fraction k/n, where
k = supporting callers (non-Artifact record clearing the min alt-read floor)
and n = ALL callers configured for the invocation (a caller absent at the
site counts as a non-support vote). Consensus mode only — rescue output is
unchanged.

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

from vcf_utils.ensemble_confidence import wilson_interval
from vcf_utils.io_utils import create_output_header

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"


class TestWilsonInterval:
    """Wilson score math (z = 1.96, 95% interval)."""

    def test_full_support(self):
        lo, hi = wilson_interval(3, 3)
        assert lo == pytest.approx(0.4385, abs=1e-3)
        assert hi == pytest.approx(1.0, abs=1e-9)

    def test_partial_support(self):
        lo, hi = wilson_interval(2, 3)
        assert lo == pytest.approx(0.2077, abs=1e-3)
        assert hi == pytest.approx(0.9385, abs=1e-3)

    def test_zero_support_lower_bound_is_zero(self):
        lo, hi = wilson_interval(0, 3)
        assert lo == 0.0
        assert hi == pytest.approx(0.5615, abs=1e-3)

    def test_single_support(self):
        lo, hi = wilson_interval(1, 3)
        assert lo == pytest.approx(0.0615, abs=1e-3)
        assert hi == pytest.approx(0.7923, abs=1e-3)

    def test_no_callers_returns_none(self):
        assert wilson_interval(0, 0) is None
        assert wilson_interval(2, 0) is None

    def test_bounds_clamped_to_unit_interval(self):
        for k in range(0, 6):
            lo, hi = wilson_interval(k, 5)
            assert 0.0 <= lo <= hi <= 1.0


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


class TestEnsembleConfidenceInfo:
    """End-to-end: run_consensus_vcf.py annotates ENS_* on every record."""

    def test_support_fraction_counts_absent_caller_as_non_support(self, consensus_vcf):
        records = _records_by_pos(consensus_vcf)
        assert records[1000].INFO.get("ENS_SUPPORT") == "2/3"
        assert records[2000].INFO.get("ENS_SUPPORT") == "3/3"

    def test_interval_bounds_match_wilson(self, consensus_vcf):
        records = _records_by_pos(consensus_vcf)
        for pos, (k, n) in [(1000, (2, 3)), (2000, (3, 3))]:
            lo, hi = wilson_interval(k, n)
            assert float(records[pos].INFO.get("ENS_CONF_LO")) == pytest.approx(
                round(lo, 4), abs=1e-6
            )
            assert float(records[pos].INFO.get("ENS_CONF_HI")) == pytest.approx(
                round(hi, 4), abs=1e-6
            )

    def test_header_declares_ens_fields(self, consensus_vcf):
        header = str(VCF(consensus_vcf).raw_header)
        for field in ("ENS_SUPPORT", "ENS_CONF_LO", "ENS_CONF_HI"):
            assert f"ID={field}" in header


class TestRescueHeaderUnchanged:
    """Rescue mode (include_rescue_fields=True) must not gain ENS_* fields."""

    def _header_info_ids(self, tmp_path, include_rescue_fields):
        template_path = _write_vcf(tmp_path, "template.vcf", MUTECT2_VCF)
        template = VCF(template_path)
        header = create_output_header(
            template, "sample", include_rescue_fields=include_rescue_fields
        )
        return set(header.info.keys())

    def test_consensus_header_has_ens_fields(self, tmp_path):
        info_ids = self._header_info_ids(tmp_path, include_rescue_fields=False)
        assert {"ENS_SUPPORT", "ENS_CONF_LO", "ENS_CONF_HI"} <= info_ids

    def test_rescue_header_has_no_ens_fields(self, tmp_path):
        info_ids = self._header_info_ids(tmp_path, include_rescue_fields=True)
        assert not {"ENS_SUPPORT", "ENS_CONF_LO", "ENS_CONF_HI"} & info_ids
