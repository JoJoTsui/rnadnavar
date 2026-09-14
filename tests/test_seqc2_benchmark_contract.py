import csv
import gzip
import importlib.util
import json
from pathlib import Path


def load_script(name):
    path = Path(__file__).parents[1] / "examples" / "seqc2" / "scripts" / name
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_vcf(path, records):
    with gzip.open(path, "wt") as out:
        out.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt, filt in records:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t.\n")


def test_stage_transition_counts_true_and_false(tmp_path):
    mod = load_script("compare_stage_transitions.py")
    before = tmp_path / "before.vcf.gz"
    after = tmp_path / "after.vcf.gz"
    truth = tmp_path / "truth.vcf.gz"
    write_vcf(before, [("chr1", 1, "A", "C", "PASS"), ("chr1", 2, "G", "T", "PASS")])
    write_vcf(after, [("chr1", 1, "A", "C", "PASS"), ("chr1", 3, "C", "G", "Somatic")])
    write_vcf(truth, [("chr1", 2, "G", "T", "PASS"), ("chr1", 3, "C", "G", "PASS")])
    assert mod.rows(before) == {("chr1", 1, "A", "C"): {"filter": "PASS", "info": "."},
                                ("chr1", 2, "G", "T"): {"filter": "PASS", "info": "."}}
    added, removed = set(mod.rows(after)) - set(mod.rows(before)), set(mod.rows(before)) - set(mod.rows(after))
    assert (len(added), len(removed), len(set(mod.rows(truth)) & added), len(set(mod.rows(truth)) & removed)) == (1, 1, 1, 1)


def test_aggregate_uses_semantic_type_labels(tmp_path):
    mod = load_script("aggregate_benchmark.py")
    metrics = {"metrics": [{"data": [
        {"id": "type", "values": ["records", "indels", "SNVs"]},
        {"id": "tp", "values": [10, 3, 7]}, {"id": "fp", "values": [1, 0, 2]},
        {"id": "fn", "values": [2, 1, 4]},
    ]}]}
    path = tmp_path / "m.json"
    path.write_text(json.dumps(metrics))
    parsed = mod.parse_metrics_json(path)
    assert list(parsed) == ["snp", "indel", "records"]
    assert parsed["snp"]["tp"] == 7
    assert parsed["indel"]["fp"] == 0
    assert parsed["records"]["fn"] == 2


def test_policy_f1_is_derived_and_zero_safe():
    mod = load_script("policy_metrics.py")
    assert mod.f1({"precision": 0.8, "recall": 0.5}, "x") == 0.6153846153846154
    assert mod.f1({"precision": 0.0, "recall": 0.0}, "x") == 0.0
    assert mod.f1({"f1": 0.7, "precision": 0.8, "recall": 0.5}, "x") == 0.7
    assert mod.f1({"tp": 5, "fp": 0, "fn": 5}, "x") == 2 / 3


def test_matrix_gate_requires_all_cells_and_comparators(tmp_path):
    mod = load_script("validate_policy_matrix.py")
    assert mod.load
