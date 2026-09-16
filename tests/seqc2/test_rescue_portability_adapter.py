from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "examples/seqc2/scripts"))
from audit_rescue_gate_portability import caller_evidence


@pytest.mark.parametrize("reverse", [False, True])
def test_repeated_low_depth_pass_does_not_become_eligible(tmp_path, reverse):
    samples = [("case-17", "0/1:0,2"), ("control-42", "0/0:40,0")]
    if reverse:
        samples.reverse()
    header = (
        '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
        '##normal_sample=control-42\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
        + '\t'.join(name for name, _ in samples) + '\n'
    )
    row = ('chr1\t10\t.\tT\tG\t50\tPASS\t.\tGT:AD\t'
           + '\t'.join(value for _, value in samples) + '\n')
    source = tmp_path / "arbitrary.vcf"
    source.write_text(header + row + row)
    key = ("chr1", 10, "T", "G")
    observed, eligible = caller_evidence(source, "mutect2", {key}, False)
    assert observed == {key}
    assert eligible == set()
    # The same positive evidence can nominate in DNA without being an RNA vote.
    assert caller_evidence(source, "mutect2", {key}, True) == {key}
    source.write_text(header + row.replace("0,2", "0,3"))
    assert caller_evidence(source, "mutect2", {key}, False) == ({key}, {key})
