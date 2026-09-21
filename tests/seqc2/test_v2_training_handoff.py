import importlib.util
from pathlib import Path
import sys
import pyarrow as pa

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'examples/seq2neo/scripts'))
spec = importlib.util.spec_from_file_location('handoff',ROOT/'examples/seq2neo/scripts/prepare_v2_training_handoff.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_quarantine_samples_and_conflicts_not_native_negatives():
    table=pa.table({'sample_id':['PRJNA298376_4255','PRJNA298330_4032','PRJNA298376_4081']+['PRJNA298376_4007']*4,
                    'FILTER':['Reference']*3+['Reference','Germline','Somatic','Somatic'],
                    'review_reason':['native_negative_requires_paired_validation']*5+
                    ['annotation_conflict:common_population_af','somatic_candidate_requires_label_qc']})
    assert module.selection(table).to_pylist()==[False,False,False,True,True,False,True]
    assert 'PRJNA298376_4007' in module.THREE
    assert not module.THREE & module.EXCLUDED
