from pathlib import Path
import importlib.util

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('reconsensus', ROOT / 'examples/seq2neo/scripts/run_reconsensus_rerun.py')
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def test_default_policy_is_native_and_gated():
    assert MOD.DEFAULTS['native_evidence_snv'] is True
    assert MOD.DEFAULTS['rescue_promotion_enabled'] is True
    assert MOD.DEFAULTS['rescue_min_dna_callers'] == 1
    assert MOD.DEFAULTS['rescue_min_rna_callers'] == 2


def test_build_command_emits_policy_flags(tmp_path):
    cfg = dict(MOD.DEFAULTS, main_nf='main.nf', rdv_conf='shared.config')
    cmd = MOD.build_command(cfg, tmp_path / 'sample.csv', tmp_path / 'out')
    joined = ' '.join(cmd)
    assert '--native_evidence_snv true' in joined
    assert '--rescue_promotion_enabled true' in joined
    assert '--rescue_min_dna_callers 1' in joined
    assert '--rescue_min_rna_callers 2' in joined
    assert '--rescue_veto dna' in joined
    assert '--tools consensus,rescue,filtering,vep' in joined


def test_manifest_paths_are_preferred_and_outcomes_classified(tmp_path):
    paths = {}
    for name in MOD.MANIFEST_VCF_FIELDS:
        path = tmp_path / f'{name}.vcf.gz'
        path.write_bytes(b'placeholder')
        (tmp_path / f'{name}.vcf.gz.tbi').write_bytes(b'index')
        paths[name] = str(path)
    row = {'base_output_dir': str(tmp_path / 'legacy'), 'dir_name': 'x', 'vcf_prefix': 'x'}
    row.update({field: paths[name] for name, field in MOD.MANIFEST_VCF_FIELDS.items()})
    found = MOD.locate_caller_vcfs(row)
    assert set(found) == set(MOD.MANIFEST_VCF_FIELDS)
    assert MOD.classify_inputs(found) == ('READY', [])
    dna_only = {k: v for k, v in found.items() if k.startswith('dna_')}
    assert MOD.classify_inputs(dna_only)[0] == 'DNA_ONLY'
    partial = dict(dna_only)
    partial.pop('dna_strelka')
    assert MOD.classify_inputs(partial)[0] == 'PARTIAL'


def test_dna_only_completion_omits_rescue_artifacts(tmp_path, monkeypatch):
    out = tmp_path / 'out'; out.mkdir()
    (out / 'pipeline_info').mkdir()
    (out / 'pipeline_info' / 'execution_trace.txt').write_text('trace')
    (out / 'consensus').mkdir(); (out / 'consensus' / 'x.vcf.gz').write_bytes(b'x')
    (out / 'filtered').mkdir(); (out / 'filtered' / 'x.filtered.vcf.gz').write_bytes(b'x')
    (out / 'filtered' / 'x.filtered.vcf.stripped.vcf.gz').write_bytes(b'x')
    monkeypatch.setattr(MOD, 'vcf_gz_sane', lambda _: True)
    ok, reason = MOD.evaluate_completion(out, dict(MOD.DEFAULTS, require_rescue=False))
    assert ok, reason
    ok, reason = MOD.evaluate_completion(out, dict(MOD.DEFAULTS, require_rescue=True))
    assert not ok and 'completion artifacts' in reason


def test_partial_rna_panel_is_not_dna_only(tmp_path):
    found = {}
    for name in ("dna_deepsomatic", "dna_mutect2", "dna_strelka", "rna_deepsomatic"):
        path = tmp_path / f"{name}.vcf.gz"
        path.write_bytes(b"x")
        (tmp_path / f"{name}.vcf.gz.tbi").write_bytes(b"i")
        found[name] = path
    status, missing = MOD.classify_inputs(found)
    assert status == "PARTIAL"
    assert "rna_mutect2" in missing
