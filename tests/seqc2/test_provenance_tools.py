import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "examples/seqc2/scripts/build_artifact_provenance.py"


def test_provenance_manifest_hashes_declared_inputs_and_preserves_metadata(tmp_path):
    artifact = tmp_path / "final.rescue.vcf.gz"; reference = tmp_path / "reference.fa"
    artifact.write_bytes(b"artifact"); reference.write_bytes(b"reference")
    metadata = tmp_path / "metadata.json"
    metadata.write_text(json.dumps({"sample": "P1", "stage": "final_second_rescue", "selector": "Somatic", "effective_args": {"min_alt_support": 0}, "tools": {"deepsomatic": "unknown"}}))
    out = tmp_path / "provenance.json"
    result = subprocess.run([sys.executable, str(SCRIPT), "--artifact", str(artifact), "--reference", str(reference), "--metadata", str(metadata), "--command", "run_rescue_vcf.py --min_alt_support 0", "--out", str(out)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    data = json.loads(out.read_text())
    assert data["artifact"]["sha256"]
    assert data["inputs"]["reference"]["sha256"]
    assert data["stage"] == "final_second_rescue"
    assert data["effective_args"]["min_alt_support"] == 0
    assert data["tools"]["deepsomatic"] == "unknown"
    assert data["commands"] == ["run_rescue_vcf.py --min_alt_support 0"]


def test_provenance_manifest_fails_on_missing_artifact(tmp_path):
    out = tmp_path / "provenance.json"
    result = subprocess.run([sys.executable, str(SCRIPT), "--artifact", str(tmp_path / "missing.vcf.gz"), "--out", str(out)], capture_output=True, text=True)
    assert result.returncode != 0
    assert "artifact" in result.stderr.lower()
    assert not out.exists()
