"""Regression checks for complete, evidence-preserving candidate exports."""
import importlib.util
from pathlib import Path
import shutil
import pyarrow.parquet as pq
import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("refined_export", ROOT/"examples/seq2neo/scripts/build_refined_rerun_artifacts.py")
export = importlib.util.module_from_spec(spec)
spec.loader.exec_module(export)

@pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools required")
def test_three_classes_keep_evidence_and_full_alt(tmp_path):
    vcf = tmp_path/"input.vcf"
    vcf.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
        + ''.join(f'##FILTER=<ID={x},Description="test">\n' for x in (*export.LABELS,"NoConsensus"))
        + '##INFO=<ID=RESCUED,Number=1,Type=String,Description="test">\n'
        + '##INFO=<ID=GNOMAD_AF,Number=1,Type=Float,Description="test">\n'
        + '##INFO=<ID=NORMAL_DP_BY_CALLER,Number=1,Type=String,Description="test">\n'
        + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        + 'chr1\t10\t.\tA\tC\t.\tSomatic\tRESCUED=YES;GNOMAD_AF=0.02\n'
        + 'chr1\t20\t.\tA\tC,G\t.\tGermline\tNORMAL_DP_BY_CALLER=mutect2:40\n'
        + 'chr1\t30\t.\tA\tAT\t.\tReference\tRESCUED=NO\n'
        + 'chr1\t40\t.\tA\tT\t.\tNoConsensus\t.\n')
    out = tmp_path/"variants.parquet"
    with pq.ParquetWriter(out, export.SCHEMA) as writer:
        counts, _, reasons = export.stream_rows(vcf, "sample", writer, batch_size=1)
    rows = pq.read_table(out).to_pylist()
    assert counts == dict(Somatic=1, Germline=1, Reference=1)
    assert rows[0]["RESCUED"] == "YES"
    assert rows[0]["review_reason"] == "annotation_conflict:common_population_af"
    assert rows[1]["ALT"] == "C,G"
    assert rows[1]["NORMAL_DP_BY_CALLER"] == "mutect2:40"
    assert rows[2]["GNOMAD_AF"] is None
    assert all(row["training_eligible"] is False for row in rows)

def test_verification_conflict_and_missing_evidence():
    assert export.review_reason("Somatic", {"DNA_VERIFICATION":"rejected"}) == "dna_verification_conflict"
    assert export.review_reason("Reference", {}) == "inherited_or_legacy_negative_requires_paired_validation"

