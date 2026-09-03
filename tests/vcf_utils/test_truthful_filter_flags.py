"""Regression tests for truthful RaVeX filter flags (audit finding C2, ticket 03).

Bug: the filtering stage read ``variant.format("AD")`` on consensus VCFs that
have no FORMAT column, so ``min_alt_reads`` fired on 100% of records; and
``vc_filter`` fired on 100% because the consensus FILTER column holds a
biological class (Somatic/Germline/...), never PASS.

Fix under test:
- Consensus VCFs carry per-caller tumor alt-count INFO (ALT_COUNT_BY_CALLER /
  ALT_COUNT_MAX), extracted from the tumor sample (ticket 02).
- ``min_alt_reads`` sources alt support from that INFO (max tumor alt count
  across callers), falling back to FORMAT/AD only for VCFs that have one.
- Biological-class FILTER values are exempt from the ``vc_filter`` check.

Run with: .venv/bin/python -m pytest tests/vcf_utils/test_truthful_filter_flags.py -v
"""

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import pysam
from cyvcf2 import VCF

# Ensure bin/ is on the Python path
_BIN_DIR = Path(__file__).resolve().parent.parent.parent / "bin"
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

from filter_vcf import apply_filters
from vcf_utils.unified_filters import apply_ravex_filters

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"

# Mutect2 is invoked with the normal CRAM first (sample 0 = normal); tumor alt = 40.
MUTECT2_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOO8801DN\tCOO8801DT
chr1\t1000\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:AF:GQ\t0/0:80,0:80:0.0:99\t0/1:60,40:100:0.4:99
"""

# Strelka writes samples in fixed NORMAL,TUMOR order with tiered base counts;
# tumor alt (G tier1) = 40.
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
"""

# Consensus-style VCF: NO FORMAT column (as written by write_union_vcf).
# Records exercise: 2-caller Somatic with good alt support, low-alt Somatic,
# biological-class FILTERs that must not trip vc_filter, and a genuine caller
# failure filter that still must.
CONSENSUS_STYLE_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=Somatic,Description="Somatic variant">
##FILTER=<ID=Artifact,Description="Artifact">
##FILTER=<ID=NoConsensus,Description="No consensus">
##FILTER=<ID=RNAedit,Description="RNA editing event">
##FILTER=<ID=LowDepth,Description="Low depth (caller rejection)">
##INFO=<ID=ALT_COUNT_MAX,Number=1,Type=Integer,Description="Maximum tumor alt-read count across callers">
##INFO=<ID=ALT_COUNT_BY_CALLER,Number=.,Type=String,Description="Tumor alt-read count per caller">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1000\t.\tA\tG\t.\tSomatic\tALT_COUNT_MAX=40;ALT_COUNT_BY_CALLER=mutect2:40|strelka:40
chr1\t2000\t.\tC\tT\t.\tSomatic\tALT_COUNT_MAX=1;ALT_COUNT_BY_CALLER=mutect2:1|strelka:1
chr1\t3000\t.\tG\tA\t.\tArtifact\tALT_COUNT_MAX=25
chr1\t4000\t.\tT\tC\t.\tNoConsensus\tALT_COUNT_BY_CALLER=mutect2:10|strelka:12
chr1\t5000\t.\tA\tT\t.\tLowDepth\tALT_COUNT_MAX=30
chr1\t6000\t.\tG\tC\t.\tRNAedit\tALT_COUNT_MAX=8
"""

# Same shape but with no alt-count evidence at all (legacy consensus VCFs):
# min_alt_reads still fires (no evidence), vc_filter must not.
CONSENSUS_NO_EVIDENCE_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=Somatic,Description="Somatic variant">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1000\t.\tA\tG\t.\tSomatic\t.
"""

# Raw caller VCF (has FORMAT/AD, no consensus INFO): AD fallback must be used.
RAW_CALLER_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOO8801DT
chr1\t1000\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:50,1
chr1\t2000\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/1:50,10
"""


def _write_vcf(directory, name, content):
    path = directory / name
    path.write_text(content)
    return str(path)


# args namespace mirroring filter_vcf.py / filter_rescue_vcf.py CLI defaults
ARGS = SimpleNamespace(
    whitelist=None,
    blacklist=None,
    gnomad_thr=0.001,
    min_alt_reads=2,
    filters=["PASS"],
    filter_multiallelic=False,
)


def _cyvcf2_flags(vcf_path):
    vcf_in = VCF(vcf_path)
    return [(v.POS, flags) for v, _, flags in apply_filters(vcf_in, ARGS, None)]


def _pysam_flags(vcf_path):
    with pysam.VariantFile(vcf_path) as vcf:
        return [(r.pos, apply_ravex_filters(r, ARGS, use_cyvcf2=False)) for r in vcf]


@pytest.fixture
def consensus_vcf(tmp_path):
    """Real consensus VCF produced from synthetic Mutect2 + Strelka caller VCFs."""
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


def _parse_by_caller(raw):
    return dict(entry.split(":", 1) for entry in raw.split("|"))


class TestConsensusAltCountInfo:
    """Ticket 03(a): consensus output carries per-caller tumor alt counts."""

    def test_alt_count_info_written(self, consensus_vcf):
        records = list(VCF(consensus_vcf))
        assert len(records) == 1
        rec = records[0]

        alt_by_caller = _parse_by_caller(_info_scalar(rec, "ALT_COUNT_BY_CALLER"))
        assert alt_by_caller == {"mutect2": "40", "strelka": "40"}
        assert int(_info_scalar(rec, "ALT_COUNT_MAX")) == 40

    def test_alt_count_info_header_declared(self, consensus_vcf):
        header = str(VCF(consensus_vcf).raw_header)
        assert 'ID=ALT_COUNT_BY_CALLER' in header
        assert 'ID=ALT_COUNT_MAX' in header


class TestFilterVcfConsensusInput:
    """Ticket 03(b,c): the original 100%-firing bug, reproduced end-to-end."""

    def test_somatic_consensus_record_has_no_bogus_flags(self, consensus_vcf):
        """A 2-caller Somatic record (tumor alt 40) must not accumulate
        min_alt_reads/vc_filter flags. Pre-fix this record got BOTH flags."""
        flags_by_pos = dict(_cyvcf2_flags(consensus_vcf))
        flags = flags_by_pos[1000]
        assert "min_alt_reads" not in flags
        assert "vc_filter" not in flags


class TestConsensusStyleFlags:
    """Flag semantics on consensus-style VCFs (no FORMAT column)."""

    @pytest.fixture
    def flags_cyvcf2(self, tmp_path):
        return dict(_cyvcf2_flags(_write_vcf(tmp_path, "c.vcf", CONSENSUS_STYLE_VCF)))

    @pytest.fixture
    def flags_pysam(self, tmp_path):
        return dict(_pysam_flags(_write_vcf(tmp_path, "c.vcf", CONSENSUS_STYLE_VCF)))

    @pytest.mark.parametrize("reader", ["flags_cyvcf2", "flags_pysam"])
    def test_high_alt_somatic_passes(self, reader, request):
        assert request.getfixturevalue(reader)[1000] == []

    @pytest.mark.parametrize("reader", ["flags_cyvcf2", "flags_pysam"])
    def test_low_alt_record_flagged(self, reader, request):
        flags = request.getfixturevalue(reader)[2000]
        assert "min_alt_reads" in flags
        assert "vc_filter" not in flags

    @pytest.mark.parametrize("reader", ["flags_cyvcf2", "flags_pysam"])
    def test_biological_class_filters_exempt_from_vc_filter(self, reader, request):
        flags = request.getfixturevalue(reader)
        # Artifact / NoConsensus / RNAedit are biological classes, not caller rejections
        assert flags[3000] == []
        assert flags[4000] == []  # ALT_COUNT_BY_CALLER alone (max 12) suffices
        assert flags[6000] == []

    def test_genuine_caller_filter_still_flagged(self, flags_cyvcf2):
        """A true caller-rejection FILTER (LowDepth) still trips vc_filter."""
        flags = flags_cyvcf2[5000]
        assert "vc_filter" in flags
        assert "min_alt_reads" not in flags

    def test_no_evidence_still_flags_min_alt_reads_but_not_vc_filter(self, tmp_path):
        """Legacy consensus VCFs without alt-count INFO and without FORMAT keep
        the min_alt_reads flag (no evidence of support), but vc_filter must not
        fire on the biological-class FILTER. Pre-fix BOTH flags fired here."""
        flags_by_pos = dict(
            _cyvcf2_flags(_write_vcf(tmp_path, "legacy.vcf", CONSENSUS_NO_EVIDENCE_VCF))
        )
        flags = flags_by_pos[1000]
        assert "min_alt_reads" in flags
        assert "vc_filter" not in flags


class TestFormatAdFallback:
    """FORMAT/AD is only a fallback for VCFs that genuinely have it."""

    def test_format_ad_fallback(self, tmp_path):
        flags_by_pos = dict(_cyvcf2_flags(_write_vcf(tmp_path, "raw.vcf", RAW_CALLER_VCF)))
        assert "min_alt_reads" in flags_by_pos[1000]  # AD alt = 1 < 2
        assert flags_by_pos[2000] == []  # AD alt = 10

    def test_format_ad_fallback_pysam(self, tmp_path):
        flags_by_pos = dict(_pysam_flags(_write_vcf(tmp_path, "raw.vcf", RAW_CALLER_VCF)))
        assert "min_alt_reads" in flags_by_pos[1000]
        assert flags_by_pos[2000] == []
