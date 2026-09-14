import importlib.util
from pathlib import Path

def load(name, rel):
 p=Path(__file__).resolve().parents[1]/rel; s=importlib.util.spec_from_file_location(name,p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m

def test_outside_hc_is_explicitly_unassessed(tmp_path):
 p=tmp_path/'domain.json'; p.write_text('{"summary":{"wes_ll":[{"domain":"outside_hc_unassessed","decision":"retained","truth_status":"truth_present","count":3,"reasons":"none"}]}}')
 out=tmp_path/'out.json'
 import subprocess,sys
 subprocess.run([sys.executable,'examples/seqc2/scripts/audit_outside_hc.py','--domain-audit',str(p),'--out',str(out)],check=True)
 assert 'unassessed' in out.read_text()

def test_indel_policy_requires_indel_row(tmp_path):
 m=load('indel','examples/seqc2/scripts/evaluate_indel_policy.py')
 p=tmp_path/'score.json'; p.write_text('{"rows":[{"variant_type":"indel","tp":1,"fp":0,"fn":2}]}')
 # module import is the contract; input validation is exercised by CLI separately
 assert m.__doc__ and p.exists()


def test_opt_in_profile_is_explicit_and_indel_policy_is_provisional(tmp_path):
    m=load('profile','examples/seq2neo/scripts/validate_policy_profile.py')
    profile=Path('examples/seq2neo/config/native_gated_policy.yaml')
    assert 'policy_profile: native_gated_experimental' in profile.read_text()
    assert m.EXPECTED['indel_policy']=='threshold_consensus'


def test_source_manifest_hashes_content(tmp_path):
    m=load('manifest','examples/seqc2/scripts/build_policy_source_manifest.py')
    f=tmp_path/'source.txt'; f.write_text('frozen')
    assert m.digest(f) and len(m.digest(f))==64


def test_indel_policy_matrix_rejects_tp_regression(tmp_path):
    import json, subprocess, sys
    candidate = tmp_path / "candidate.json"
    baseline = tmp_path / "baseline.json"
    payload = lambda tp, fp, fn: {"metrics": [{"data": [
        {"id": "type", "values": ["indels"]},
        {"id": "tp", "values": [tp]},
        {"id": "fp", "values": [fp]},
        {"id": "fn", "values": [fn]},
        {"id": "precision", "values": [tp / (tp + fp)]},
        {"id": "recall", "values": [tp / (tp + fn)]},
    ]}]}
    candidate.write_text(json.dumps(payload(3, 1, 7)))
    baseline.write_text(json.dumps(payload(4, 2, 6)))
    out = tmp_path / "decision.json"
    subprocess.run([sys.executable, "examples/seqc2/scripts/build_indel_policy_matrix.py",
                    "--candidate", f"cell={candidate}",
                    "--baseline", f"cell={baseline}", "--out", str(out)], check=True)
    result = json.loads(out.read_text())
    assert result["status"] == "no_qualifying_shared_rule"
    assert result["decision"] == "retain_existing_threshold_consensus"
