import json, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def test_policy_selector_requires_all_slices(tmp_path):
 m={'required_slices':['WES-LL:SNV'],'partitions':{'development':'dev','holdout':'holdout'},'baseline':{'WES-LL:SNV':{'precision':.9,'recall':.8}},'candidates':[{'policy':{'x':1},'rows':[{'pair':'WES-LL','variant_type':'SNV','precision':.91,'recall':.81}]}]}
 i=tmp_path/'m.json'; o=tmp_path/'o.json'; i.write_text(json.dumps(m)); subprocess.run([sys.executable,str(ROOT/'examples/seqc2/scripts/select_development_policy.py'),'--metrics',str(i),'--out',str(o)],check=True); assert json.loads(o.read_text())['status']=='qualified'
def test_heldout_refuses_without_frozen_policy(tmp_path):
 f=tmp_path/'f.json'; m=tmp_path/'m.json'; o=tmp_path/'o.json'; f.write_text(json.dumps({'status':'no_qualifying_policy'})); m.write_text(json.dumps({'rows':[{'pair':'WES-LL','variant_type':'SNV'}]})); subprocess.run([sys.executable,str(ROOT/'examples/seqc2/scripts/evaluate_heldout_policy.py'),'--frozen',str(f),'--metrics',str(m),'--out',str(o)],check=True); assert json.loads(o.read_text())['status']=='incomplete_no_frozen_policy'


def run_selector(tmp_path, metrics):
 i=tmp_path/'m.json'; o=tmp_path/'o.json'; i.write_text(json.dumps(metrics))
 return subprocess.run([sys.executable,str(ROOT/'examples/seqc2/scripts/select_development_policy.py'),'--metrics',str(i),'--out',str(o)],capture_output=True,text=True), o

def test_policy_selector_rejects_empty_baseline(tmp_path):
 result, _ = run_selector(tmp_path, {'required_slices':['WES-LL:SNV'],'partitions':{'development':'dev','holdout':'holdout'},'baseline':{},'candidates':[{'policy':{},'rows':[]}]})
 assert result.returncode != 0 and 'baseline' in result.stderr.lower()

def test_policy_selector_rejects_duplicate_and_nonfinite_slices(tmp_path):
 base={'required_slices':['WES-LL:SNV'],'partitions':{'development':'dev','holdout':'holdout'},'baseline':{'WES-LL:SNV':{'precision':.9,'recall':.8}},'candidates':[{'policy':{},'rows':[{'pair':'WES-LL','variant_type':'SNV','precision':.91,'recall':.81},{'pair':'WES-LL','variant_type':'SNV','precision':.91,'recall':.81}]}]}
 result, _ = run_selector(tmp_path, base)
 assert result.returncode != 0 and 'duplicate' in result.stderr.lower()

def test_policy_selector_rejects_missing_required_slice(tmp_path):
 base={'required_slices':['WES-LL:SNV','WES-IL:SNV'],'partitions':{'development':'dev','holdout':'holdout'},'baseline':{'WES-LL:SNV':{'precision':.9,'recall':.8}},'candidates':[{'policy':{},'rows':[{'pair':'WES-LL','variant_type':'SNV','precision':.91,'recall':.81}]}]}
 result, _ = run_selector(tmp_path, base)
 assert result.returncode != 0 and 'required' in result.stderr.lower()


def test_policy_selector_rejects_nonfinite_metric(tmp_path):
 base={'required_slices':['WES-LL:SNV'],'partitions':{'development':'dev','holdout':'holdout'},'baseline':{'WES-LL:SNV':{'precision':.9,'recall':.8}},'candidates':[{'policy':{},'rows':[{'pair':'WES-LL','variant_type':'SNV','precision':float('inf'),'recall':.81}]}]}
 result, _ = run_selector(tmp_path, base)
 assert result.returncode != 0 and 'finite' in result.stderr.lower()


def test_heldout_recomputes_gates_instead_of_trusting_flags(tmp_path):
 frozen={'status':'qualified','frozen_policy':{'x':1},'required_slices':['WES-LL:SNV'],'partitions':{'development':'dev','holdout':'holdout'},'min_delta':0.0,'baseline':{'WES-LL:SNV':{'precision':.9,'recall':.8}}}
 held={'required_slices':['WES-LL:SNV'],'rows':[{'pair':'WES-LL','variant_type':'SNV','precision':.89,'recall':.81,'precision_gate':True,'recall_gate':True}]}
 f=tmp_path/'f.json'; m=tmp_path/'m.json'; o=tmp_path/'o.json'; f.write_text(json.dumps(frozen)); m.write_text(json.dumps(held))
 subprocess.run([sys.executable,str(ROOT/'examples/seqc2/scripts/evaluate_heldout_policy.py'),'--frozen',str(f),'--metrics',str(m),'--out',str(o)],check=True)
 result=json.loads(o.read_text()); assert result['status']=='complete'; assert result['rows'][0]['accepted'] is False; assert result['rows'][0]['precision_gate'] is False

def test_heldout_rejects_unbound_or_incomplete_slices(tmp_path):
 frozen={'status':'qualified','frozen_policy':{'x':1},'required_slices':['WES-LL:SNV'],'partitions':{'development':'dev','holdout':'holdout'},'min_delta':0.0,'baseline':{'WES-LL:SNV':{'precision':.9,'recall':.8}}}
 held={'required_slices':['WES-IL:SNV'],'rows':[]}
 f=tmp_path/'f.json'; m=tmp_path/'m.json'; o=tmp_path/'o.json'; f.write_text(json.dumps(frozen)); m.write_text(json.dumps(held))
 result=subprocess.run([sys.executable,str(ROOT/'examples/seqc2/scripts/evaluate_heldout_policy.py'),'--frozen',str(f),'--metrics',str(m),'--out',str(o)],capture_output=True,text=True)
 assert result.returncode != 0 and 'slice' in result.stderr.lower()
