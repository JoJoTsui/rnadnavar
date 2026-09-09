"""Rescue source independence and reported/output outcome reconciliation."""
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin/run_rescue_vcf.py"


@pytest.fixture
def inputs(tmp_path):
    header = '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
    columns = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    consensus = (header + '##FILTER=<ID=NoConsensus,Description="Insufficient support">\n'
                 + columns + '\nchr1\t100\t.\tA\tG\t50\tNoConsensus\t.\n')
    caller = (header
              + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              + '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
              + '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
              + '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">\n'
              + columns + '\tFORMAT\tTUMOR\n'
              + 'chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP:AD:AF\t0/1:20:12,8:0.4\n')
    for name, content in (("dna.vcf", consensus), ("rna.vcf", consensus),
                          ("dna.mutect2.vcf", caller), ("rna.mutect2.vcf", caller),
                          ("other.mutect2.vcf", caller)):
        (tmp_path / name).write_text(content)
    return tmp_path


def run_rescue(directory, *extra):
    return subprocess.run([
        sys.executable, str(SCRIPT), '--dna_consensus', str(directory / 'dna.vcf'),
        '--rna_consensus', str(directory / 'rna.vcf'), '--out_prefix', str(directory / 'result'),
        '--output_format', 'vcf', *map(str, extra),
    ], capture_output=True, text=True)


@pytest.mark.parametrize('second', ['dna.mutect2.vcf', 'other.mutect2.vcf'])
def test_duplicate_caller_identity_rejected_before_output(inputs, second):
    result = run_rescue(inputs, '--dna_vcf', inputs / 'dna.mutect2.vcf',
                        '--dna_vcf', inputs / second)
    assert result.returncode == 2
    assert 'duplicate DNA caller identity' in result.stderr
    assert not (inputs / 'result.vcf').exists()


def test_missing_caller_rejected_before_output(inputs):
    result = run_rescue(inputs, '--rna_vcf', inputs / 'missing.strelka.vcf')
    assert result.returncode == 2
    assert 'RNA caller VCF not found' in result.stderr
    assert not (inputs / 'result.vcf').exists()


def test_promotion_statistics_match_written_record(inputs):
    result = run_rescue(inputs, '--dna_vcf', inputs / 'dna.mutect2.vcf',
                        '--rna_vcf', inputs / 'rna.mutect2.vcf')
    assert result.returncode == 0, result.stderr
    with pysam.VariantFile(inputs / 'result.vcf') as reader:
        records = list(reader)
        assert len(records) == 1
        assert list(records[0].filter) == ['Somatic']
        for field in ('RESCUE_PROMOTED', 'RESCUED'):
            value = records[0].info[field]
            assert (value[0] if isinstance(value, tuple) else value) == 'YES'
    assert 'Rescued variants (cross-modality support): 1' in result.stdout
