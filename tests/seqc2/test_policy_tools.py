import json, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def test_policy_selector_requires_all_slices(tmp_path):
 m={'baseline':{'WES-LL:SNV':{'precision':.9,'recall':.8}},'candidates':[{'policy':{'x':1},'rows':[{'pair':'WES-LL','variant_type':'SNV','precision':.91,'recall':.81}]}]}
 i=tmp_path/'m.json'; o=tmp_path/'o.json'; i.write_text(json.dumps(m)); subprocess.run([sys.executable,str(ROOT/'examples/seqc2/scripts/select_development_policy.py'),'--metrics',str(i),'--out',str(o)],check=True); assert json.loads(o.read_text())['status']=='qualified'
def test_heldout_refuses_without_frozen_policy(tmp_path):
 f=tmp_path/'f.json'; m=tmp_path/'m.json'; o=tmp_path/'o.json'; f.write_text(json.dumps({'status':'no_qualifying_policy'})); m.write_text(json.dumps({'rows':[{'pair':'WES-LL','variant_type':'SNV'}]})); subprocess.run([sys.executable,str(ROOT/'examples/seqc2/scripts/evaluate_heldout_policy.py'),'--frozen',str(f),'--metrics',str(m),'--out',str(o)],check=True); assert json.loads(o.read_text())['status']=='incomplete_no_frozen_policy'
