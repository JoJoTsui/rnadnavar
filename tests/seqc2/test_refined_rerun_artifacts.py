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
        + '##INFO=<ID=NEGATIVE_EVIDENCE_STATUS,Number=1,Type=String,Description="test">\n'
        + '##INFO=<ID=NEGATIVE_NORMAL_COUNTS,Number=1,Type=String,Description="test">\n'
        + '##INFO=<ID=THREE_CLASS_POLICY,Number=1,Type=String,Description="test">\n'
        + '##INFO=<ID=THREE_CLASS_NATIVE_RATIONALE,Number=1,Type=String,Description="test">\n'
        + '##INFO=<ID=TRAINING_ELIGIBLE,Number=1,Type=String,Description="test">\n'
        + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        + 'chr1\t10\t.\tA\tC\t.\tSomatic\tRESCUED=YES;GNOMAD_AF=0.02;THREE_CLASS_POLICY=separated_three_class_v2;THREE_CLASS_NATIVE_RATIONALE=native_trace;TRAINING_ELIGIBLE=NO\n'
        + 'chr1\t20\t.\tA\tC,G\t.\tGermline\tNORMAL_DP_BY_CALLER=mutect2:40\n'
        + 'chr1\t30\t.\tA\tAT\t.\tReference\tRESCUED=NO;NEGATIVE_EVIDENCE_STATUS=SUPPORTED;NEGATIVE_NORMAL_COUNTS=299:0:0:299\n'
        + 'chr1\t40\t.\tA\tT\t.\tNoConsensus\t.\n')
    out = tmp_path/"variants.parquet"
    with pq.ParquetWriter(out, export.SCHEMA) as writer:
        counts, _, reasons = export.stream_rows(vcf, "sample", writer, batch_size=1)
    rows = pq.read_table(out).to_pylist()
    assert counts == dict(Somatic=1, Germline=1, Reference=1)
    assert rows[0]["RESCUED"] == "YES"
    assert rows[0]["review_reason"] == "annotation_conflict:common_population_af"
    assert rows[0]['THREE_CLASS_POLICY']=='separated_three_class_v2'
    assert rows[0]['THREE_CLASS_NATIVE_RATIONALE']=='native_trace'
    assert rows[0]['TRAINING_ELIGIBLE']=='NO'
    assert rows[1]["ALT"] == "C,G"
    assert rows[1]["NORMAL_DP_BY_CALLER"] == "mutect2:40"
    assert rows[2]["GNOMAD_AF"] is None
    assert rows[2]["NEGATIVE_NORMAL_COUNTS"] == "299:0:0:299"
    assert rows[2]["review_reason"] == "negative_read_supported_biological_approval_pending"
    assert all(row["training_eligible"] is False for row in rows)

def test_verification_conflict_and_missing_evidence():
    assert export.review_reason("Somatic", {"DNA_VERIFICATION":"rejected"}) == "dna_verification_conflict"
    assert export.review_reason("Reference", {}) == "inherited_or_legacy_negative_requires_paired_validation"
    assert export.review_reason('Reference',{'THREE_CLASS_POLICY':'separated_three_class_v2'})=='native_negative_requires_paired_validation'
    assert export.review_reason('Somatic',{'THREE_CLASS_REVIEW_REASON':'native_negative_conflict'})=='three_class_review:native_negative_conflict'


def test_new_policy_export_never_falls_back_to_old_baseline(tmp_path):
    state={'output':str(tmp_path),'identity':{'policy':'separated_three_class_v2'}}
    with pytest.raises(ValueError,match='refuse baseline fallback'):export.candidate_paths(state)
    state['final_artifacts']={k:str(tmp_path/f'three_class.{k}.vcf.gz') for k in ('consensus','rescue')}
    rescue,consensus=export.candidate_paths(state)
    assert rescue.name=='three_class.rescue.vcf.gz'
    assert consensus.name=='three_class.consensus.vcf.gz'
    state['final_artifacts']['rescue']=str(tmp_path/'refined.rescue.vcf.gz')
    with pytest.raises(ValueError,match='old baseline'):export.candidate_paths(state)
    assert export.candidate_paths({'output':str(tmp_path)})==(tmp_path/'refined.rescue.vcf.gz',tmp_path/'refined.vcf.gz')


def test_new_counts_never_keep_old_output_paths_or_hashes(tmp_path):
    consensus=tmp_path/'three_class.consensus.vcf.gz'; rescue=tmp_path/'three_class.rescue.vcf.gz'
    old=dict(new_consensus='/old/consensus.vcf.gz',new_rescue='/old/rescue.vcf.gz',
             new_consensus_sha256='old-c',new_rescue_sha256='old-r',new_rescue_records=999,
             original_rescue='/original/rescue.vcf.gz',original_sha256='original')
    state=dict(identity={'policy':'separated_three_class_v2'},outputs={str(consensus):'new-c',str(rescue):'new-r'})
    result=export.sample_summary(old,state,consensus,rescue,{'Somatic':1},{'Somatic':1,'Reference':2},
                                 {'Somatic':1,'Reference':2},{},{})
    assert result['new_consensus']==str(consensus) and result['new_consensus_sha256']=='new-c'
    assert result['new_rescue']==str(rescue) and result['new_rescue_sha256']=='new-r'
    assert result['new_rescue_records']==3
    assert result['original_rescue']==old['original_rescue'] and result['original_sha256']=='original'
