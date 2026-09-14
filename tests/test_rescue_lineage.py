import importlib.util
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"examples/seq2neo/scripts/validate_rescue_lineage.py"
spec=importlib.util.spec_from_file_location("lineage", SCRIPT)
lineage=importlib.util.module_from_spec(spec)
spec.loader.exec_module(lineage)

def test_final_realign_filtered_path_passes():
    rows=lineage.validate_rows([{
        "sample_id":"S1",
        "rescue_vcf_path":"output/vcf_realignment/rescue/S1_RT_realign.filtered.vcf.stripped.vcf.gz",
    }])
    assert rows[0]["status"]=="pass"

def test_intermediate_first_round_path_is_rejected():
    rows=lineage.validate_rows([{
        "sample_id":"S1", "rescue_vcf_path":"output/rescue/S1.rescued.vcf.gz",
    }])
    assert rows[0]["status"]=="error"
    assert "intermediate rescue artifact" in rows[0]["reasons"]
    assert "realignment lineage is not identifiable" in rows[0]["reasons"]

def test_first_round_can_be_explicitly_allowed_but_must_be_final_filtered():
    rows=lineage.validate_rows([{
        "sample_id":"S1", "rescue_vcf_path":"output/rescue/S1.filtered.vcf.gz",
    }], require_realign=False)
    assert rows[0]["status"]=="pass"
