import csv
import importlib.util
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('export_verification',ROOT/'examples/seq2neo/scripts/verify_three_class_export.py')
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


@pytest.fixture
def exported(tmp_path):
    path=tmp_path/'candidates.parquet'
    rows=[dict(sample_id='one',FILTER=cls,variant_type='SNP',POS=i+1,training_eligible=False,
               TRAINING_ELIGIBLE='NO',THREE_CLASS_POLICY='separated_three_class_v2')
          for i,cls in enumerate(('Somatic','Germline','Reference'))]
    pq.write_table(pa.Table.from_pylist(rows),path)
    summary=dict(parquet=str(path),parquet_sha256=m.sha(path),samples=1,candidate_rows=3,
                 class_counts={r['FILTER']:1 for r in rows},class_type_counts={r['FILTER']+':SNP':1 for r in rows},
                 samples_detail=[dict(sample_id='one',new_rescue='new.rescue',new_consensus='new.consensus',
                                      exported_counts={r['FILTER']:1 for r in rows})])
    summary_path=tmp_path/'summary.json';summary_path.write_text(json.dumps(summary))
    manifest=tmp_path/'manifest.tsv'
    row=dict(sample_id='one',training_label_vcf='',label_qc_verdict='NOT_APPROVED',candidate_policy='separated_three_class_v2',
             rescue_vcf_path='new.rescue',consensus_vcf_path='new.consensus',variant_parquet_path=str(path))
    with manifest.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(row),delimiter='\t');writer.writeheader();writer.writerow(row)
    return summary_path,manifest,tmp_path/'verification.json'


def test_every_class_and_fresh_destination(exported):
    result=m.verify(*exported)
    assert result['rows']==3 and result['samples']==1 and not result['training_approved']
    with pytest.raises(ValueError,match='fresh'):m.verify(*exported)


@pytest.mark.parametrize('column,value', [('training_eligible',True),('TRAINING_ELIGIBLE','YES'),
                                         ('THREE_CLASS_POLICY','wrong'),('POS',0),('FILTER','NoConsensus')])
def test_rejects_invalid_export_rows(exported,column,value):
    summary=json.loads(exported[0].read_text());path=Path(summary['parquet'])
    rows=pq.read_table(path).to_pylist();rows[0][column]=value
    pq.write_table(pa.Table.from_pylist(rows),path)
    summary['parquet_sha256']=m.sha(path);exported[0].write_text(json.dumps(summary))
    with pytest.raises(ValueError):m.verify(*exported)
    assert not exported[2].exists()


def test_rejects_count_mismatch(exported):
    summary=json.loads(exported[0].read_text());summary['class_counts']['Reference']=2
    exported[0].write_text(json.dumps(summary))
    with pytest.raises(ValueError,match='class-count'):m.verify(*exported)


def test_rejects_checksum_mismatch(exported):
    summary=json.loads(exported[0].read_text());summary['parquet_sha256']='0'*64
    exported[0].write_text(json.dumps(summary))
    with pytest.raises(ValueError,match='checksum'):m.verify(*exported)
