import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[2] / "examples" / "seqc2" / "scripts" / "aggregate_benchmark.py"
spec = importlib.util.spec_from_file_location("aggregate_benchmark", SCRIPT)
aggregate_benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aggregate_benchmark)


def metrics_doc(labels):
    return {
        "metrics": [
            {
                "data": [
                    {"id": "type", "values": labels},
                    {"id": "tp", "values": [10, 20, 30]},
                    {"id": "fp", "values": [1, 2, 3]},
                    {"id": "fn", "values": [4, 5, 6]},
                ]
            }
        ]
    }


def test_metrics_are_keyed_by_semantic_type_not_position(tmp_path):
    path = tmp_path / "query.metrics.json"
    path.write_text(json.dumps(metrics_doc(["records", "SNVs", "indels"])))
    parsed = aggregate_benchmark.parse_metrics_json(path)
    assert parsed["snp"]["tp"] == 20
    assert parsed["indel"]["tp"] == 30
    assert parsed["records"]["tp"] == 10


def test_missing_or_duplicate_type_labels_fail_closed(tmp_path):
    path = tmp_path / "query.metrics.json"
    path.write_text(json.dumps(metrics_doc(["SNVs", "SNVs", "records"])))
    with pytest.raises(ValueError, match="invalid or duplicate"):
        aggregate_benchmark.parse_metrics_json(path)
