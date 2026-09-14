import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "examples/seqc2/scripts/validate_ticket_coverage.py"
spec = importlib.util.spec_from_file_location("ticket_coverage", SCRIPT)
coverage = importlib.util.module_from_spec(spec)
spec.loader.exec_module(coverage)


def test_matrix_preserves_known_gates_without_failure():
    report = {
        "status": "pass_with_known_gates",
        "checks": [
            {"name": name, "status": "pass"}
            for name in (
                "benchmark:wes_ll/ukb", "benchmark:wes_ll/medexome",
                "benchmark:wgs_il/ukb", "benchmark:wgs_il/medexome",
                "read-evidence:available", "default:native-flag",
                "policy:experimental-not-promoted", "policy:no-production-nomination-wiring",
                "process-reachability:available", "focused-tests",
                "ingress-contract:available", "rerun:checksum-output",
                "cohort-rerun:checksum-output",
            )
        ],
    }
    matrix = coverage.build_matrix(report)
    assert [row["ticket"] for row in matrix] == ["01", "02", "03", "04", "05", "06", "07"]
    assert matrix[-1]["status"] == "pass_with_known_gates"


def test_missing_evidence_is_inconclusive():
    matrix = coverage.build_matrix({"status": "pass_with_known_gates", "checks": []})
    assert matrix[0]["status"] == "inconclusive"
    assert matrix[6]["status"] == "pass_with_known_gates"
