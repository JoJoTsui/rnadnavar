"""
Regression tests for tumor-sample-aware genotype extraction (audit finding C1).

Bug: extract_genotype_info() in bin/vcf_utils/aggregation.py always read
sample index 0 (the NORMAL for Mutect2 and Strelka in this pipeline), and the
Strelka path picked DP and VAF from *different* samples (max-read-count row).
Consensus VCFs therefore carried normal-sample GT/AD/DP/VAF, e.g. a somatic
site with tumor VAF 0.40 emitted CONSENSUS_GT=0/0 and
VAF_BY_CALLER=mutect2:0.0000.

These tests build synthetic paired caller VCFs (normal first, tumor VAF 0.40)
and assert that per-caller GT/DP/AD/VAF and the mean aggregates come from the
TUMOR sample of each caller's VCF.
"""

import subprocess
import sys
from pathlib import Path

import pytest
from cyvcf2 import VCF

from vcf_utils.aggregation import (
    read_variants_from_vcf,
    resolve_tumor_sample_index,
)

BIN_DIR = Path(__file__).resolve().parents[2] / "bin"

# Mutect2 in this pipeline is invoked with the normal CRAM first, so the
# normal is sample 0. Sample names come from BAM read groups and do NOT
# contain "tumor"/"normal" substrings (e.g. COO8801DN / COO8801DT).
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

# Strelka writes samples in fixed NORMAL,TUMOR order and uses non-standard
# per-sample tiered base counts (AU/CU/GU/TU) instead of AD/AF.
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

# DeepSomatic names samples <id>_tumor / <id>_normal and lists the tumor first.
DEEPSOMATIC_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fraction">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTEST_tumor\tTEST_normal
chr1\t1000\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:VAF:GQ\t0/1:60,40:100:0.4:99\t0/0:80,0:80:0.0:99
"""

# Mutect2 paired VCF whose header declares ##normal_sample AND whose sample
# order is TUMOR-first, deliberately violating the pipeline's normal-first
# ordering convention. The header is ground truth and must win.
MUTECT2_NORMAL_HEADER_TUMOR_FIRST_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000000>
##normal_sample=COO8801DN
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOO8801DT\tCOO8801DN
chr1\t1000\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:AF:GQ\t0/1:60,40:100:0.4:99\t0/0:80,0:80:0.0:99
"""

# Tumor-only (single-sample) Mutect2 run: the only sample must be used.
MUTECT2_TUMOR_ONLY_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOO8801DT
chr1\t1000\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:AF\t0/1:60,40:100:0.4
"""

TUMOR_GT = "0/1"
TUMOR_DP = 100
TUMOR_AD = "60,40"
TUMOR_VAF = 0.4


def _write_vcf(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(content)
    return str(path)


@pytest.fixture
def mutect2_vcf(tmp_path):
    return _write_vcf(tmp_path, "sample.mutect2.variants.vcf", MUTECT2_VCF)


@pytest.fixture
def strelka_vcf(tmp_path):
    return _write_vcf(tmp_path, "sample.strelka.variants.vcf", STRELKA_VCF)


@pytest.fixture
def deepsomatic_vcf(tmp_path):
    return _write_vcf(tmp_path, "sample.deepsomatic.variants.vcf", DEEPSOMATIC_VCF)


def _read_single_genotype(vcf_path, caller):
    variants = read_variants_from_vcf(vcf_path, caller)
    assert len(variants) == 1
    return next(iter(variants.values()))["genotype"]


def _assert_tumor_genotype(gt):
    assert gt["GT"] == TUMOR_GT, f"GT should come from the tumor sample: {gt}"
    assert gt["DP"] == TUMOR_DP, f"DP should come from the tumor sample: {gt}"
    assert gt["AD"] == TUMOR_AD, f"AD should come from the tumor sample: {gt}"
    assert gt["VAF"] == pytest.approx(TUMOR_VAF), (
        f"VAF should come from the tumor sample: {gt}"
    )


class TestResolveTumorSampleIndex:
    def test_strelka_fixed_order_by_name(self):
        assert resolve_tumor_sample_index(["NORMAL", "TUMOR"], "strelka") == 1

    def test_deepsomatic_names(self):
        assert (
            resolve_tumor_sample_index(["TEST_tumor", "TEST_normal"], "deepsomatic")
            == 0
        )

    def test_mutect2_normal_not_first_in_name_but_first_in_order(self):
        # Pipeline convention: normal CRAM is passed first, read-group names
        # carry no tumor/normal hint -> tumor is the last sample.
        assert resolve_tumor_sample_index(["COO8801DN", "COO8801DT"], "mutect2") == 1

    def test_named_normal_picks_the_other_sample(self):
        assert resolve_tumor_sample_index(["TUMOUR", "NORMAL"], "strelka") == 0
        assert resolve_tumor_sample_index(["NORMAL", "CASE1"], "mutect2") == 1

    def test_single_sample(self):
        assert resolve_tumor_sample_index(["COO8801DT"], "mutect2") == 0

    def test_unknown_caller_falls_back_to_zero(self):
        assert resolve_tumor_sample_index(["A", "B"], "unknowncaller") == 0

    def test_normal_sample_header_overrides_ordering_convention(self):
        # Mutect2 header declares the normal; the samples are TUMOR-first,
        # violating the normal-first ordering convention. The header wins.
        assert (
            resolve_tumor_sample_index(
                ["COO8801DT", "COO8801DN"], "mutect2", normal_sample="COO8801DN"
            )
            == 0
        )

    def test_normal_sample_header_absent_or_unknown_is_ignored(self):
        # No header value -> legacy resolution; unknown name -> legacy resolution
        assert resolve_tumor_sample_index(["COO8801DN", "COO8801DT"], "mutect2") == 1
        assert (
            resolve_tumor_sample_index(
                ["COO8801DN", "COO8801DT"], "mutect2", normal_sample="NOT_IN_VCF"
            )
            == 1
        )


class TestExtractGenotypeInfo:
    def test_mutect2_tumor_sample(self, mutect2_vcf):
        _assert_tumor_genotype(_read_single_genotype(mutect2_vcf, "mutect2"))

    def test_strelka_tumor_sample(self, strelka_vcf):
        # DP and VAF must come from the SAME (tumor) sample, not DP from the
        # normal (index 0) and VAF from the max-read-count row.
        _assert_tumor_genotype(_read_single_genotype(strelka_vcf, "strelka"))

    def test_deepsomatic_tumor_sample(self, deepsomatic_vcf):
        _assert_tumor_genotype(_read_single_genotype(deepsomatic_vcf, "deepsomatic"))

    def test_single_sample_tumor_only(self, tmp_path):
        vcf = _write_vcf(tmp_path, "tumor_only.mutect2.variants.vcf",
                         MUTECT2_TUMOR_ONLY_VCF)
        _assert_tumor_genotype(_read_single_genotype(vcf, "mutect2"))

    def test_mutect2_normal_sample_header_tumor_first(self, tmp_path):
        # ##normal_sample header resolves the tumor even when the sample
        # order violates the normal-first Mutect2 convention.
        vcf = _write_vcf(
            tmp_path,
            "sample.mutect2.tumor_first.variants.vcf",
            MUTECT2_NORMAL_HEADER_TUMOR_FIRST_VCF,
        )
        _assert_tumor_genotype(_read_single_genotype(vcf, "mutect2"))


def _info_scalar(variant, key):
    val = variant.INFO.get(key)
    if isinstance(val, (tuple, list)):
        val = val[0] if val else None
    return val


def _parse_by_caller(raw):
    return dict(entry.split(":", 1) for entry in raw.split("|"))


class TestConsensusCliAggregates:
    """End-to-end regression for the exact symptom reported in audit C1."""

    def test_consensus_vcf_carries_tumor_genotypes(
        self, tmp_path, mutect2_vcf, strelka_vcf
    ):
        out_prefix = tmp_path / "out.consensus"
        cmd = [
            sys.executable,
            str(BIN_DIR / "run_consensus_vcf.py"),
            "--input_dir",
            str(tmp_path),
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

        records = list(VCF(str(out_prefix) + ".vcf"))
        assert len(records) == 1
        rec = records[0]

        assert _info_scalar(rec, "CONSENSUS_GT") == TUMOR_GT

        gt_by_caller = _parse_by_caller(_info_scalar(rec, "GT_BY_CALLER"))
        assert gt_by_caller == {"mutect2": TUMOR_GT, "strelka": TUMOR_GT}

        dp_by_caller = _parse_by_caller(_info_scalar(rec, "DP_BY_CALLER"))
        assert dp_by_caller == {
            "mutect2": str(TUMOR_DP),
            "strelka": str(TUMOR_DP),
        }

        vaf_by_caller = _parse_by_caller(_info_scalar(rec, "VAF_BY_CALLER"))
        assert vaf_by_caller == {"mutect2": "0.4000", "strelka": "0.4000"}

        assert float(_info_scalar(rec, "VAF_MEAN")) == pytest.approx(TUMOR_VAF)
        assert float(_info_scalar(rec, "DP_MEAN")) == pytest.approx(TUMOR_DP)
