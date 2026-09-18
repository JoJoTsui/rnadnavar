import sys
from pathlib import Path
import pysam
import pytest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'bin'))
from apply_three_class_labels import decide_label, run, digest
from assess_negative_label_evidence import evidence_index, native_nomination


def vcf(path, rows, native=False):
    text='##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
    text+='##INFO=<ID=CLASSIFICATION_RATIONALE,Number=1,Type=String,Description="trace">\n'
    text+='##INFO=<ID=GNOMAD_AF,Number=1,Type=Float,Description="af">\n'
    text+=''.join(f'##FILTER=<ID={x},Description="label">\n' for x in ['Somatic','Germline','Reference','NoConsensus'])
    text+='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    for pos,label in rows:
        trace='three_class_policy:native_three_class_v1' if native else 'rule:seqc2_refined_v2'
        text+=f'chr1\t{pos}\t.\tA\tT\t.\t{label}\tCLASSIFICATION_RATIONALE={trace};GNOMAD_AF=0.1\n'
    path.write_text(text)


def test_keep_somatic_even_with_native_negative_but_no_new_somatic():
    assert decide_label('Somatic','Reference')[0]=='Somatic'
    assert decide_label('Somatic','Germline')[0]=='Somatic'
    assert decide_label('NoConsensus','Somatic')[0]=='NoConsensus'
    assert decide_label('Reference',None)[0]=='NoConsensus'


def test_complete_union_provenance_and_somatic_parity(tmp_path):
    b,n=tmp_path/'b.vcf',tmp_path/'n.vcf'
    vcf(b,[(1,'Somatic'),(2,'Reference'),(3,'Germline')])
    vcf(n,[(1,'Reference'),(2,'Germline'),(4,'Reference')],True)
    hashes=digest(b),digest(n)
    report=run(b,n,tmp_path/'out',*hashes,'realignment')
    assert report['somatic_membership_mismatches']==0
    assert report['sources_unchanged']
    with pysam.VariantFile(report['output']) as f:
        rows=list(f)
    assert [next(iter(r.filter)) for r in rows]==['Somatic','Germline','NoConsensus','Reference']
    assert all(r.info['TRAINING_ELIGIBLE']=='NO' for r in rows)
    assert 'native_negative_conflict' in rows[0].info['THREE_CLASS_REVIEW_REASON']
    assert native_nomination(rows[1].info)
    assert (digest(b),digest(n))==hashes
    with pytest.raises(FileExistsError):run(b,n,tmp_path/'out',*hashes,'realignment')


def test_reject_unnominated_negative_and_wrong_checksum(tmp_path):
    b,n=tmp_path/'b.vcf',tmp_path/'n.vcf'
    vcf(b,[(1,'Somatic')]);vcf(n,[(2,'Reference')])
    with pytest.raises(ValueError,match='checksum'):run(b,n,tmp_path/'bad','wrong',digest(n),'consensus')
    with pytest.raises(ValueError,match='provenance'):run(b,n,tmp_path/'bad2',digest(b),digest(n),'consensus')


def test_evidence_must_attest_primary_independent_read_filters():
    report=dict(vcf_sha256='a',thresholds=dict(MAPQ=20,BQ=20,BAQ=True),results=[])
    with pytest.raises(ValueError,match='primary'):evidence_index(report,'a')
