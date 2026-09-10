#!/usr/bin/env python3
"""
run_pipeline.py
───────────────
Run the rnadnavar pipeline for the SEQC2 DNA-only benchmark: loads
config/runner.yaml, builds the micromamba/nextflow command, executes it,
and records run state for resume.

Modeled on examples/seq2neo/scripts/run_batch_from_json.py, minus the
JSON-manifest/CSV-generation machinery (seqc2 uses one static samplesheet).

Usage:
  python3 scripts/run_pipeline.py --config config/runner.yaml
  python3 scripts/run_pipeline.py --config config/runner.yaml --dry-run
  python3 scripts/run_pipeline.py --config config/runner.yaml --outdir /path/to/out
"""

import argparse
import json
import os
import hashlib
import subprocess
import sys
from datetime import datetime
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("PyYAML required:  pip install pyyaml")

EXAMPLE_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = EXAMPLE_DIR.parent.parent

DEFAULTS = {
    "main_nf": "",           # empty -> REPO_ROOT/main.nf
    "rdv_conf": "",          # empty -> EXAMPLE_DIR/seqc2.shared.config
    "nxf_conda_cachedir": "",
    "nxf_conda_usemamba": "true",
    "micromamba_env": "nextflow",
    "https_proxy": "",
    "seqc2_root": "",        # empty -> EXAMPLE_DIR
    "input_csv": "csv/seqc2_wes_ll.csv",
    "outdir": "output/seqc2.wes.ll",
    "state_file": "runs/run_state.json",
    "dry_run": False,
    "resume": True,
    "offline": True,
    "completion_artifacts": [
        "**/pipeline_info/execution_trace*.txt",
        "consensus/**/*.consensus.vcf.gz",
    ],
    "preflight_validator": [],
    "completion_validator": [],
    "provenance_file": "",
}


def load_config(config_path, cli_overrides: dict) -> dict:
    cfg = dict(DEFAULTS)
    if config_path and Path(config_path).exists():
        with open(config_path) as f:
            file_cfg = yaml.safe_load(f) or {}
        cfg.update({k: v for k, v in file_cfg.items() if v is not None})
    cfg.update({k: v for k, v in cli_overrides.items() if v is not None})
    if not cfg["seqc2_root"]:
        # Hybrid profiles live below examples/seqc2/hybrid; anchor their
        # relative inputs/configs to that directory so invocation CWD is free.
        config_parent = Path(config_path).resolve().parent if config_path else EXAMPLE_DIR
        cfg["seqc2_root"] = str(config_parent if config_parent.name == "hybrid" else EXAMPLE_DIR)
    cfg["main_nf"] = cfg["main_nf"] or str(REPO_ROOT / "main.nf")
    if cfg["rdv_conf"]:
        rdv = Path(cfg["rdv_conf"])
        cfg["rdv_conf"] = str(rdv if rdv.is_absolute() else Path(cfg["seqc2_root"]) / rdv)
    else:
        cfg["rdv_conf"] = str(EXAMPLE_DIR / "seqc2.shared.config")
    return cfg


def resolve(cfg: dict, key: str) -> Path:
    """Resolve path relative to seqc2_root when not absolute."""
    p = Path(cfg[key])
    return p if p.is_absolute() else Path(cfg["seqc2_root"]) / p


def build_command(cfg: dict, input_csv: Path, outdir: Path) -> list:
    cmd = ["micromamba", "run", "-n", cfg["micromamba_env"]]
    cmd += ["nextflow", "run", cfg["main_nf"], "-c", cfg["rdv_conf"]]
    cmd += ["--input", str(input_csv)]
    cmd += ["--outdir", str(outdir)]
    if cfg.get("step"):
        cmd += ["--step", str(cfg["step"])]
    if cfg.get("tools"):
        cmd += ["--tools", str(cfg["tools"])]
    if cfg.get("offline"):
        cmd.append("-offline")
    if cfg.get("resume"):
        cmd.append("-resume")
    cmd.append("-with-conda")
    return cmd


def build_env(cfg: dict) -> dict:
    env = os.environ.copy()
    if cfg.get("https_proxy"):
        for var in ("HTTPS_PROXY", "https_proxy", "HTTP_PROXY", "http_proxy"):
            env[var] = cfg["https_proxy"]
    if cfg.get("nxf_conda_cachedir"):
        env["NXF_CONDA_CACHEDIR"] = cfg["nxf_conda_cachedir"]
    if cfg.get("nxf_conda_usemamba"):
        env["NXF_CONDA_USEMAMBA"] = str(cfg["nxf_conda_usemamba"])
    return env


def run_validator(command: list, outdir: Path, input_csv: Path | None = None) -> int:
    """Run an optional contract validator without exposing its output on success."""
    if not command:
        return 0
    substitutions = {"outdir": str(outdir)}
    if input_csv is not None:
        substitutions["input_csv"] = str(input_csv)
    cmd = [str(part).format(**substitutions) for part in command]
    return subprocess.run(cmd, text=True).returncode


def is_complete(outdir: Path, artifacts: list, cfg: dict | None = None) -> bool:
    """Require published artifacts and an optional final-artifact validator."""
    if not outdir.exists():
        return False
    if not all(list(outdir.glob(p)) for p in artifacts):
        return False
    return run_validator((cfg or {}).get("completion_validator", []), outdir) == 0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_provenance(cfg: dict, input_csv: Path, outdir: Path) -> None:
    if not cfg.get("provenance_file"):
        return
    trace = sorted((outdir / "pipeline_info").glob("execution_trace*.txt"), key=lambda p: p.stat().st_mtime)
    cached = 0
    if trace:
        cached = sum(1 for line in trace[-1].read_text().splitlines()[1:] if "\tCACHED\t" in line)
    data = {
        "input": {"path": str(input_csv), "sha256": sha256(input_csv)},
        "configuration": {"path": cfg["rdv_conf"], "sha256": sha256(Path(cfg["rdv_conf"]))},
        "code": {"main_nf": cfg["main_nf"], "sha256": sha256(Path(cfg["main_nf"]))},
        "outdir": str(outdir), "trace": str(trace[-1]) if trace else None,
        "cached_tasks": cached,
        "completion_artifacts": {pattern: [str(p) for p in outdir.glob(pattern)] for pattern in cfg["completion_artifacts"]},
    }
    path = resolve(cfg, "provenance_file")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Run the SEQC2 DNA-only pipeline run (single samplesheet)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--config", default=str(EXAMPLE_DIR / "config" / "runner.yaml"))
    ap.add_argument("--input", default=None, dest="input_csv",
                    help="Override samplesheet path")
    ap.add_argument("--outdir", default=None, help="Override nextflow outdir")
    ap.add_argument("--dry-run", action="store_true", default=None, dest="dry_run")
    ap.add_argument("--no-resume", action="store_false", default=None, dest="resume")
    args = ap.parse_args()

    cli_overrides = {k: v for k, v in vars(args).items() if k != "config"}
    cfg = load_config(args.config, cli_overrides)

    input_csv = resolve(cfg, "input_csv")
    outdir = resolve(cfg, "outdir")
    state_path = resolve(cfg, "state_file")
    dry_run = bool(cfg.get("dry_run"))

    for label, p in [("main_nf", Path(cfg["main_nf"])),
                     ("rdv_conf", Path(cfg["rdv_conf"])),
                     ("input_csv", input_csv)]:
        if not p.exists():
            sys.exit(f"ERROR: {label} not found: {p}")

    key = input_csv.stem
    state = json.loads(state_path.read_text()) if state_path.exists() else {}

    if is_complete(outdir, cfg["completion_artifacts"], cfg):
        print(f"[SKIP]  {key}  already complete (artifacts present in {outdir})")
        return

    cmd = build_command(cfg, input_csv, outdir)
    print(f"[{'DRY' if dry_run else 'RUN'}]  {key}")
    print(f"       conf: {cfg['rdv_conf']}")
    print(f"       csv : {input_csv}")
    print(f"       out : {outdir}")
    print(f"       cmd : {' '.join(cmd)}")

    if dry_run:
        return

    if run_validator(cfg.get("preflight_validator", []), outdir, input_csv) != 0:
        sys.exit("ERROR: preflight validation failed; no Nextflow tasks were launched")

    state[key] = {"status": "running", "started": datetime.now().isoformat()}
    state_path.parent.mkdir(parents=True, exist_ok=True)
    state_path.write_text(json.dumps(state, indent=2))

    try:
        subprocess.run(cmd, env=build_env(cfg), check=True)
        validator_status = run_validator(cfg.get("completion_validator", []), outdir)
        if validator_status == 3:
            state[key] = {"status": "zero_candidates", "finished": datetime.now().isoformat()}
            write_provenance(cfg, input_csv, outdir)
            print(f"[ZERO]  {key}  no candidates reached realignment")
            sys.exit(3)
        if is_complete(outdir, cfg["completion_artifacts"], cfg):
            state[key] = {"status": "succeeded",
                          "finished": datetime.now().isoformat()}
            write_provenance(cfg, input_csv, outdir)
            print(f"[OK]    {key}")
        else:
            state[key] = {"status": "failed",
                          "reason": "completion artifacts missing after run"}
            print(f"[FAIL]  {key}  completion artifacts missing after run")
            sys.exit(1)
    except subprocess.CalledProcessError as e:
        state[key] = {"status": "failed", "reason": str(e)}
        print(f"[FAIL]  {key}  {e}")
        sys.exit(e.returncode or 1)
    finally:
        state_path.write_text(json.dumps(state, indent=2))


if __name__ == "__main__":
    main()
