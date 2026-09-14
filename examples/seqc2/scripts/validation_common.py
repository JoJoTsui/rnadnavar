"""Shared status and report primitives for validation-only scripts."""
from pathlib import Path
import json

PASS = "pass"
ERROR = "error"
BLOCKED = "blocked"
INCONCLUSIVE = "inconclusive"
KNOWN_GATES = "pass_with_known_gates"

def check(name, passed, detail="", *, failure=ERROR):
    return {"name": name, "status": PASS if passed else failure, "detail": detail}

def write_report(outdir, filename, report):
    outdir=Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    (outdir/filename).write_text(json.dumps(report, indent=2)+"\n")

def overall(checks, *, known_gates=False):
    if any(item["status"] in {ERROR, BLOCKED} for item in checks): return BLOCKED
    if known_gates: return KNOWN_GATES
    return PASS
