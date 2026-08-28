"""Regression tests for gnomAD AF INFO lookup (audit finding C3).

The gnomAD annotator (bin/vcf_utils/gnomad_annotator.py) writes
``INFO/GNOMAD_AF`` (uppercase, copied from gnomAD's Number=A AF field).
The RaVeX gnomAD guard probed only lowercase/mixed-case spellings, so the
``gnomad`` filter flag never fired. These tests run annotated VCF fixtures
through the rescue-filter path and assert the guard fires for every
historical field spelling, case-insensitively.

Run with: .venv/bin/python -m pytest tests/vcf_utils/test_gnomad_af_lookup.py -v
"""

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

from vcf_utils.unified_filters import apply_ravex_filters, get_gnomad_af
from filter_vcf import get_gnomad_af_cyvcf2

# Spellings seen in the wild: GNOMAD_AF is what gnomad_annotator.py writes;
# the rest are historical aliases probed by the filter.
GNOMAD_AF_FIELD_NAMES = [
    "GNOMAD_AF",
    "gnomAD_AF",
    "AF_gnomad",
    "gnomad_AF",
    "MAX_AF",
    "GnomAD_af",  # arbitrary case — matching must be case-insensitive
]

VCF_TEMPLATE = """\
##fileformat=VCFv4.2
##contig=<ID=chr1>
##INFO=<ID={field},Number=A,Type=Float,Description="gnomAD allele frequency">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tSomatic\t{field}=0.5
chr1\t200\t.\tC\tT\t.\tSomatic\t{field}=0.00001
"""

VCF_NO_GNOMAD = """\
##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tSomatic\t.
"""


def _write_vcf(tmp_path, content):
    vcf_path = tmp_path / "annotated.vcf"
    vcf_path.write_text(content)
    return str(vcf_path)


def _pysam_records(vcf_path):
    with pysam.VariantFile(vcf_path) as vcf:
        return list(vcf)


# args namespace mirroring filter_rescue_vcf.py CLI defaults relevant here;
# min_alt_reads=0 isolates the gnomAD rule from the alt-read rule.
ARGS = SimpleNamespace(gnomad_thr=0.001, min_alt_reads=0)


@pytest.mark.parametrize("field", GNOMAD_AF_FIELD_NAMES)
def test_common_polymorphism_is_flagged_pysam(tmp_path, field):
    """A record annotated with AF=0.5 must be flagged 'gnomad' by the rescue filter."""
    vcf_path = _write_vcf(tmp_path, VCF_TEMPLATE.format(field=field))
    high_af, low_af = _pysam_records(vcf_path)

    high_flags = apply_ravex_filters(high_af, ARGS, use_cyvcf2=False)
    low_flags = apply_ravex_filters(low_af, ARGS, use_cyvcf2=False)

    assert "gnomad" in high_flags
    assert "gnomad" not in low_flags


@pytest.mark.parametrize("field", GNOMAD_AF_FIELD_NAMES)
def test_get_gnomad_af_pysam(tmp_path, field):
    """get_gnomad_af returns the annotated value regardless of field spelling."""
    vcf_path = _write_vcf(tmp_path, VCF_TEMPLATE.format(field=field))
    high_af, low_af = _pysam_records(vcf_path)

    assert get_gnomad_af(high_af, use_cyvcf2=False) == pytest.approx(0.5)
    assert get_gnomad_af(low_af, use_cyvcf2=False) == pytest.approx(0.00001)


@pytest.mark.parametrize("field", GNOMAD_AF_FIELD_NAMES)
def test_get_gnomad_af_cyvcf2(tmp_path, field):
    """filter_vcf.py's cyvcf2 lookup has the same case-insensitive behaviour."""
    vcf_path = _write_vcf(tmp_path, VCF_TEMPLATE.format(field=field))
    records = list(VCF(vcf_path))

    assert get_gnomad_af_cyvcf2(records[0]) == pytest.approx(0.5)
    assert get_gnomad_af_cyvcf2(records[1]) == pytest.approx(0.00001)


def test_missing_gnomad_field_returns_zero(tmp_path):
    """Records without any gnomAD annotation yield AF 0.0 and no flag."""
    vcf_path = _write_vcf(tmp_path, VCF_NO_GNOMAD)
    (record,) = _pysam_records(vcf_path)

    assert get_gnomad_af(record, use_cyvcf2=False) == 0.0
    assert "gnomad" not in apply_ravex_filters(record, ARGS, use_cyvcf2=False)
