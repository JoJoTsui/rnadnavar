import importlib.util
from pathlib import Path

import pytest
import json
import subprocess
import sys
import pysam

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('negative_gate', ROOT/'bin/assess_negative_label_evidence.py')
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


def reads(ref, alt=0, other=0):
    return dict(ref=ref,alt=alt,other=other,depth=ref+alt+other)


def test_reference_detection_limit_and_no_summed_depth():
    assert gate.required_zero_alt_depth(.01,.95)==299
    assert gate.required_zero_alt_depth(.05,.95)==59
    assert gate.assess('Reference','A','T',reads(299),reads(299))[0]=='SUPPORTED'
    assert gate.assess('Reference','A','T',reads(298),reads(999))[0]=='WITHHELD'
    assert gate.assess('Reference','A','T',reads(500),reads(500,1))[0]=='WITHHELD'
    assert gate.assess('Reference','A','T',reads(39,2),reads(41,2))[0]=='WITHHELD'


@pytest.mark.parametrize('bad',[None,{},reads(0),dict(depth=20,ref=20,alt=None,other=0),dict(depth=20,ref=18,alt=0,other=0)])
def test_missing_invalid_reads_never_become_reference(bad):
    assert gate.assess('Reference','A','T',bad,reads(500))[0]=='WITHHELD'


def test_germline_and_ambiguous_context_are_separate():
    assert gate.assess('Germline','A','T',reads(20,20),reads(20,20))[0]=='SUPPORTED'
    assert gate.assess('Germline','G','GT',reads(63,57),reads(1,86))[0]=='WITHHELD'
    assert gate.assess('Germline','A','T',reads(20,20),reads(0,40))[0]=='WITHHELD'
    assert gate.assess('Germline','A','T',reads(60),reads(20,20))[0]=='CONFLICT'
    assert gate.assess('Somatic','A','T',None,None)[0]=='NOT_APPLICABLE'
    for context in ('LOH','loh','somatic_on_germline','unknown',{'invalid':'schema'}):
        assert gate.assess('Germline','A','T',reads(20,20),reads(20,20),context_review=context)[0]=='WITHHELD'


@pytest.mark.parametrize('value',[0,1,float('nan'),float('inf'),-.1])
def test_bad_detection_limit(value):
    with pytest.raises(ValueError): gate.required_zero_alt_depth(value,.95)


def test_evidence_binding_and_duplicates():
    r=dict(vcf_sha256='abc',thresholds=dict(MAPQ=20,BQ=20,BAQ=True,flag_filter=0xF04,ignore_overlaps=True,ignore_orphans=True),results=[])
    with pytest.raises(ValueError):gate.evidence_index(r,'other')
    r['results']=[dict(site=['chr1',1,'A','T'])]*2
    with pytest.raises(ValueError):gate.evidence_index(r,'abc')


def test_cli_preserves_labels_but_never_grants_training_approval(tmp_path):
    vcf=tmp_path/'input.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
        '##FILTER=<ID=Reference,Description="candidate">\n'
        '##INFO=<ID=CLASSIFICATION_RATIONALE,Number=1,Type=String,Description="rule">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        'chr1\t10\t.\tA\tT\t.\tReference\tCLASSIFICATION_RATIONALE=three_class_policy:native_three_class_v1\n'
        'chr1\t20\t.\tA\tT\t.\tReference\tCLASSIFICATION_RATIONALE=three_class_policy:native_three_class_v1\n')
    before=gate.digest(vcf)
    evidence=tmp_path/'evidence.json'
    evidence.write_text(json.dumps(dict(vcf_sha256=before,thresholds=dict(MAPQ=20,BQ=20,BAQ=True,flag_filter=0xF04,ignore_overlaps=True,ignore_orphans=True),results=[
        dict(site=['chr1',10,'A','T'],**{'class':'Reference'},normal=reads(299),tumor=reads(299))])))
    cmd=[sys.executable,str(ROOT/'bin/assess_negative_label_evidence.py'),'--vcf',str(vcf),
         '--bam-evidence',str(evidence),'--outdir',str(tmp_path/'out')]
    subprocess.run(cmd,check=True,capture_output=True,text=True)
    with pysam.VariantFile(str(tmp_path/'out/candidate.evidence.vcf.gz')) as f: records=list(f)
    assert [set(r.filter) for r in records]==[{'Reference'},{'Reference'}]
    assert [r.info['NEGATIVE_EVIDENCE_STATUS'] for r in records]==['SUPPORTED','WITHHELD']
    assert all(r.info['TRAINING_ELIGIBLE']=='NO' for r in records)
    assert records[0].info['NEGATIVE_NORMAL_COUNTS']=='299:0:0:299'
    assert 'NEGATIVE_NORMAL_COUNTS' not in records[1].info
    assert gate.digest(vcf)==before
    assert subprocess.run(cmd,capture_output=True).returncode != 0
