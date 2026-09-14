"""Guard parity with the historical experiment, not production policy choices."""
import importlib.util
from pathlib import Path
import sys

SCRIPTS = Path(__file__).resolve().parents[1] / "examples/seqc2/scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("historical_replay", SCRIPTS / "replay_historical_native_gate.py")
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


def row(filt="PASS", qual="30", evidence="."):
    return ["chr1", "10", ".", "A", "C", qual, filt, evidence]


def test_historical_native_retains_pass_indels_but_does_not_promote_indels():
    accepted = ("chr1", 10, "A", "AT")
    rejected = ("chr1", 20, "C", "CT")
    ds = {accepted: row(), rejected: row("RefCall")}
    m2 = {rejected: row(evidence="TLOD=20;GERMQ=80")}
    assert replay.historical_native(ds, m2) == {accepted}


def test_historical_veto_also_removes_deepsomatic_pass():
    key = ("chr1", 10, "A", "C")
    ds = {key: row()}
    m2 = {key: row("contamination;orientation;weak_evidence", evidence="TLOD=20;GERMQ=80")}
    assert replay.historical_native(ds, m2) == set()


def test_additions_require_deepsomatic_presence_and_positive_quality():
    key = ("chr1", 10, "A", "C")
    m2 = {key: row(evidence="TLOD=12;GERMQ=60")}
    assert replay.historical_native({}, m2) == set()
    assert replay.historical_native({key: row("RefCall", "0")}, m2) == set()
    assert replay.historical_native({key: row("RefCall", "1")}, m2) == {key}


def test_gate_preserves_baseline_and_requires_somatic_snv_with_both_counts():
    baseline = ("chr1", 10, "A", "AT")
    good = ("chr1", 20, "A", "C")
    missing = ("chr1", 30, "A", "G")
    artifact = ("chr1", 40, "C", "G")
    indel = ("chr1", 50, "C", "CT")
    evidence = "N_DNA_CALLERS_SUPPORT=1;N_RNA_CALLERS_SOMATIC=2"
    rescue = {good: row("Somatic", evidence=evidence),
              missing: row("Somatic", evidence="N_RNA_CALLERS_SOMATIC=2"),
              artifact: row("Artifact", evidence=evidence),
              indel: row("Somatic", evidence=evidence)}
    assert replay.gated({baseline}, rescue) == {baseline, good}
