#!/usr/bin/env python3
"""
run_batch_from_json.py
──────────────────────
Batch runner: reads merged.json, generates per-sample nextflow CSV inputs,
executes the rnadnavar pipeline, and tracks run state for resume.

Config (YAML, all values overridable via CLI):
  See config/runner.example.yaml for all options.

Filters (all optional, combinable):
  --project   PRJNA298376
  --disease   "Colon cancer"     (substring, case-insensitive)
  --patient   4007
  --set       1                  (partition set 1-4)
  --status    standard           (standard | extra | all)

Usage examples:
  # Dry-run all eligible samples
  python3 scripts/run_batch_from_json.py --config config/runner.yaml --dry-run

  # Run set1 (colorectal) standard samples only
  python3 scripts/run_batch_from_json.py --config config/runner.yaml --set 1 --status standard

  # Run a single patient
  python3 scripts/run_batch_from_json.py --config config/runner.yaml \\
      --project PRJNA298376 --patient 3812

  # Resume (re-run same command — finished samples are auto-skipped)
  python3 scripts/run_batch_from_json.py --config config/runner.yaml --set 2
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("PyYAML required:  pip install pyyaml")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.common import (
    MOD_STATUS_CODE,
    REQUIRED_MODALITIES,
    is_eligible,
    sample_key,
)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

DEFAULTS = {
    "main_nf": "",
    "rdv_conf": "",
    "nxf_conda_cachedir": "",
    "nxf_conda_usemamba": "true",
    "micromamba_env": "nextflow",
    "https_proxy": "",
    "seq2neo_root": str(Path(__file__).resolve().parent.parent),
    "merged_json": "data/processed/merged.json",
    "csv_dir": "runs/csv",
    "outdir_base": "output",
    "state_file": "runs/run_state.json",
    "lane": "LX",
    "dry_run": False,
    "resume": True,
    "offline": True,
    "max_parallel": 1,
    "max_retries": 1,
    # Generic completion artifacts: all must exist.
    "completion_artifacts": [
        "**/pipeline_info/execution_trace*.txt",
    ],
    # Rescue artifacts: at least one must exist for successful completion.
    "rescue_success_patterns": [
        "rescue/**/*.filtered.vcf.gz",
        "vcf_realignment/rescue/**/*.filtered.vcf.gz",
        "rescue/**/*.rescued.vcf.gz",
        "vcf_realignment/rescue/**/*.rescued.vcf.gz",
    ],
    # Treat failed RNA-branch tasks in trace as failed sample, even if nextflow exits 0.
    "fail_on_failed_trace": True,
    "trace_file_pattern": "**/pipeline_info/execution_trace*.txt",
    "failed_trace_statuses": ["FAILED"],
    "failed_trace_process_regex": (
        "STAR_ALIGN|FASTQ_ALIGN_STAR|RNA_REALIGNMENT_WORKFLOW|SECOND_RESCUE_WORKFLOW"
        "|VCF_RESCUE|RNA_FILTERING|FILTER_RNA_MUTATIONS"
    ),
}


def load_config(config_path, cli_overrides: dict) -> dict:
    cfg = dict(DEFAULTS)
    if config_path and Path(config_path).exists():
        with open(config_path) as f:
            file_cfg = yaml.safe_load(f) or {}
        cfg.update({k: v for k, v in file_cfg.items() if v is not None})
    cfg.update({k: v for k, v in cli_overrides.items() if v is not None})
    return cfg


def resolve(cfg: dict, key: str) -> Path:
    """Resolve path relative to seq2neo_root when not absolute."""
    p = Path(cfg[key])
    return p if p.is_absolute() else Path(cfg["seq2neo_root"]) / p


# ---------------------------------------------------------------------------
# State management
# ---------------------------------------------------------------------------


def load_state(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def save_state(path: Path, state: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, indent=2))


# ---------------------------------------------------------------------------
# Completion check (strict)
# ---------------------------------------------------------------------------


def is_complete(outdir: Path, artifacts: list) -> bool:
    """outdir must exist AND every artifact glob must match ≥1 file."""
    if not outdir.exists():
        return False
    return all(list(outdir.glob(p)) for p in artifacts)


def has_any_match(outdir: Path, patterns: list) -> bool:
    """outdir must exist AND at least one glob pattern must match ≥1 file."""
    if not outdir.exists():
        return False
    return any(list(outdir.glob(p)) for p in patterns)


def latest_trace_file(outdir: Path, trace_pattern: str) -> Path | None:
    """Return the newest trace file under outdir matching the configured pattern."""
    traces = list(outdir.glob(trace_pattern))
    if not traces:
        return None
    traces.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return traces[0]


def trace_has_failed_rna_process(outdir: Path, cfg: dict) -> tuple[bool, str]:
    """
    Parse the newest execution trace and flag FAILED rows in RNA-critical processes.
    Returns: (has_failure, reason)
    """
    trace_file = latest_trace_file(outdir, cfg["trace_file_pattern"])
    if not trace_file:
        return False, ""

    failed_status = {s.upper() for s in cfg.get("failed_trace_statuses", ["FAILED"])}
    proc_re = re.compile(
        cfg.get("failed_trace_process_regex", "STAR_ALIGN"), re.IGNORECASE
    )
    failures = []

    with trace_file.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = set(reader.fieldnames or [])
        # Standard Nextflow trace has "status" and either "process" or "name"
        status_key = "status" if "status" in fields else None
        process_key = "process" if "process" in fields else None
        name_key = "name" if "name" in fields else None
        exit_key = "exit" if "exit" in fields else None

        if not status_key or (not process_key and not name_key):
            return False, ""

        for row in reader:
            status = (row.get(status_key) or "").upper()
            if status not in failed_status:
                continue
            proc_text = " ".join(
                [
                    row.get(process_key, "") if process_key else "",
                    row.get(name_key, "") if name_key else "",
                ]
            )
            if proc_re.search(proc_text):
                exit_code = row.get(exit_key, "") if exit_key else ""
                if exit_code:
                    failures.append(f"{proc_text} (exit {exit_code})")
                else:
                    failures.append(proc_text)

    if failures:
        preview = "; ".join(failures[:2])
        if len(failures) > 2:
            preview += f"; +{len(failures) - 2} more"
        return True, f"RNA branch failed in trace: {preview}"
    return False, ""


def evaluate_completion(outdir: Path, has_rna: bool, cfg: dict) -> tuple[bool, str]:
    """Evaluate sample completion using strict rescue and trace-aware criteria."""
    if not is_complete(outdir, cfg["completion_artifacts"]):
        return False, "completion artifacts missing"

    if not has_any_match(outdir, cfg["rescue_success_patterns"]):
        return False, "rescue artifacts missing"

    if has_rna and cfg.get("fail_on_failed_trace", True):
        failed_rna, reason = trace_has_failed_rna_process(outdir, cfg)
        if failed_rna:
            return False, reason

    return True, ""


# ---------------------------------------------------------------------------
# CSV generation
# ---------------------------------------------------------------------------


def write_sample_csv(sample: dict, csv_path: Path, lane: str):
    """
    Generate nextflow input CSV for one sample.
    CSV columns: patient, status, sample, lane, fastq_1, fastq_2
    status codes: DN=0, DT=1, RT=2
    For 'extra' samples only the first pair of each modality is used.
    """
    project_id = sample["_project_id"]
    patient_id = sample["patient_id"]
    nf_patient = f"{project_id}_{patient_id}"

    rows = ["patient,status,sample,lane,fastq_1,fastq_2"]
    for mod in REQUIRED_MODALITIES:
        mod_data = sample["modalities"].get(mod)
        if not mod_data or not mod_data["pairs"]:
            continue
        pair = mod_data["pairs"][0]  # always first pair
        rows.append(
            f"{nf_patient},{MOD_STATUS_CODE[mod]},{nf_patient}{mod},"
            f"{lane},{pair['r1']},{pair['r2']}"
        )

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.write_text("\n".join(rows) + "\n")


# ---------------------------------------------------------------------------
# Nextflow command builder
# ---------------------------------------------------------------------------


def build_command(cfg: dict, csv_path: Path, outdir: Path) -> list:
    cmd = ["micromamba", "run", "-n", cfg["micromamba_env"]]
    cmd += ["nextflow", "run"]
    if cfg.get("main_nf"):
        cmd.append(cfg["main_nf"])
    if cfg.get("rdv_conf"):
        cmd += ["-c", cfg["rdv_conf"]]
    cmd += ["--input", str(csv_path)]
    cmd += ["--outdir", str(outdir)]
    if cfg.get("offline"):
        cmd.append("-offline")
    if cfg.get("resume"):
        cmd.append("-resume")
    cmd.append("-with-conda")
    return cmd


def build_env(cfg: dict) -> dict:
    env = os.environ.copy()
    if cfg.get("https_proxy"):
        env["HTTPS_PROXY"] = cfg["https_proxy"]
        env["https_proxy"] = cfg["https_proxy"]
    if cfg.get("nxf_conda_cachedir"):
        env["NXF_CONDA_CACHEDIR"] = cfg["nxf_conda_cachedir"]
    if cfg.get("nxf_conda_usemamba"):
        env["NXF_CONDA_USEMAMBA"] = str(cfg["nxf_conda_usemamba"])
    return env


def format_shell_command(cfg: dict, cmd: list) -> str:
    """Format the full shell-equivalent command including env var prefixes."""
    env_prefix = []
    if cfg.get("https_proxy"):
        env_prefix.append(f'HTTPS_PROXY="{cfg["https_proxy"]}"')
    if cfg.get("nxf_conda_cachedir"):
        env_prefix.append(f'NXF_CONDA_CACHEDIR="{cfg["nxf_conda_cachedir"]}"')
    if cfg.get("nxf_conda_usemamba"):
        env_prefix.append(f"NXF_CONDA_USEMAMBA={cfg['nxf_conda_usemamba']}")
    # join env vars on one line, then the command as a single line
    env_str = " ".join(env_prefix)
    cmd_str = " ".join(cmd)
    return f"{env_str} \\\n  {cmd_str}" if env_str else cmd_str


def validate_config(cfg: dict):
    """Abort with a clear message if required pipeline paths are not set."""
    missing = []
    if not cfg.get("main_nf"):
        missing.append("main_nf")
    if not cfg.get("rdv_conf"):
        missing.append("rdv_conf")
    if missing:
        sys.exit(
            f"ERROR: required config values not set: {missing}\n"
            f"Edit config/runner.yaml and set:\n"
            + "\n".join(f"  {k}: /path/to/..." for k in missing)
        )


# ---------------------------------------------------------------------------
# Sample filtering
# ---------------------------------------------------------------------------


def filter_samples(all_samples: list, args) -> list:
    """Apply project/disease/patient/set/status filters. Only eligible samples pass."""
    out = []
    for s in all_samples:
        if not is_eligible(s):
            continue
        if args.project and s["_project_id"] != args.project:
            continue
        if args.disease and args.disease.lower() not in s["disease"].lower():
            continue
        if args.patient and str(s["patient_id"]) != str(args.patient):
            continue
        if args.set is not None and s.get("partition_set") != args.set:
            continue
        if args.status and args.status != "all" and s["status"] != args.status:
            continue
        out.append(s)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(
        description="Batch nextflow runner from merged JSON",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--config", default=None, help="Path to runner.yaml")
    # filters
    ap.add_argument("--project", default=None, help="Filter by project ID")
    ap.add_argument("--disease", default=None, help="Filter by disease (substring)")
    ap.add_argument("--patient", default=None, help="Filter by patient ID")
    ap.add_argument(
        "--set", default=None, type=int, help="Filter by partition set (1-4)"
    )
    ap.add_argument(
        "--status", default="all", help="standard | extra | all  (default: all)"
    )
    # config overrides
    ap.add_argument("--main-nf", default=None, dest="main_nf")
    ap.add_argument("--rdv-conf", default=None, dest="rdv_conf")
    ap.add_argument("--outdir-base", default=None, dest="outdir_base")
    ap.add_argument("--seq2neo", default=None, dest="seq2neo_root")
    ap.add_argument("--dry-run", action="store_true", default=None, dest="dry_run")
    ap.add_argument("--no-resume", action="store_false", default=None, dest="resume")
    args = ap.parse_args()

    # build config: file → CLI overrides
    filter_keys = {"project", "disease", "patient", "set", "status", "config"}
    cli_overrides = {k: v for k, v in vars(args).items() if k not in filter_keys}
    cfg = load_config(args.config, cli_overrides)

    # load merged JSON
    merged_path = resolve(cfg, "merged_json")
    if not merged_path.exists():
        sys.exit(
            f"merged.json not found: {merged_path}\n"
            f"Run parse_projects_to_json.py first."
        )

    data = json.loads(merged_path.read_text())
    all_samples = [s for proj in data["projects"] for s in proj["samples"]]

    selected = filter_samples(all_samples, args)
    if not selected:
        print("No samples match the given filters.")
        return

    state_path = resolve(cfg, "state_file")
    state = load_state(state_path)
    csv_dir = resolve(cfg, "csv_dir")
    outdir_base = resolve(cfg, "outdir_base")
    dry_run = bool(cfg.get("dry_run"))

    print(f"\n{'[DRY RUN] ' if dry_run else ''}Selected {len(selected)} sample(s)\n")

    counts = {"succeeded": 0, "skipped": 0, "failed": 0}

    for s in selected:
        proj = s["_project_id"]
        pid = s["patient_id"]
        key = sample_key(proj, pid)
        nf_patient = f"{proj}_{pid}"
        outdir = outdir_base / nf_patient
        csv_path = csv_dir / f"{nf_patient}.csv"
        has_rna = bool((s.get("modalities", {}).get("RT") or {}).get("pairs"))

        # ── completion check ──────────────────────────────────────────────
        complete_ok, complete_reason = evaluate_completion(outdir, has_rna, cfg)
        if complete_ok:
            print(f"[SKIP]  {key}  already complete")
            state[key] = {"status": "succeeded", "reason": "artifacts present"}
            counts["skipped"] += 1
            continue
        if outdir.exists() and complete_reason:
            print(f"[NOTE]  {key}  not complete: {complete_reason}")

        prev = state.get(key, {})
        retries = prev.get("retries", 0)
        if prev.get("status") == "failed" and retries >= cfg["max_retries"]:
            print(f"[SKIP]  {key}  max retries ({cfg['max_retries']}) reached")
            counts["skipped"] += 1
            continue

        # ── generate CSV ──────────────────────────────────────────────────
        write_sample_csv(s, csv_path, cfg["lane"])

        cmd = build_command(cfg, csv_path, outdir)
        env = build_env(cfg)

        tag = f"set{s.get('partition_set')} | {s['disease']} | {s['status']}"
        print(f"[{'DRY' if dry_run else 'RUN'}]  {key}  {tag}")
        print(f"       csv : {csv_path}")
        print(f"       out : {outdir}")
        print(f"       cmd : {format_shell_command(cfg, cmd)}")

        if dry_run:
            continue

        # ── validate config before first real execution ───────────────────
        validate_config(cfg)

        # ── execute ───────────────────────────────────────────────────────
        state[key] = {
            "status": "running",
            "started": datetime.now().isoformat(),
            "retries": retries,
        }
        save_state(state_path, state)

        try:
            subprocess.run(cmd, env=env, check=True)
            complete_ok, complete_reason = evaluate_completion(outdir, has_rna, cfg)
            if complete_ok:
                state[key] = {
                    "status": "succeeded",
                    "finished": datetime.now().isoformat(),
                }
                print(f"[OK]    {key}")
                counts["succeeded"] += 1
            else:
                state[key] = {
                    "status": "failed",
                    "reason": complete_reason or "artifacts missing after run",
                    "retries": retries + 1,
                }
                print(
                    f"[FAIL]  {key}  {complete_reason or 'artifacts missing after run'}"
                )
                counts["failed"] += 1
        except subprocess.CalledProcessError as e:
            state[key] = {"status": "failed", "reason": str(e), "retries": retries + 1}
            print(f"[FAIL]  {key}  {e}")
            counts["failed"] += 1

        save_state(state_path, state)

    print(
        f"\n=== Done  succeeded={counts['succeeded']}  "
        f"skipped={counts['skipped']}  failed={counts['failed']} ==="
    )


if __name__ == "__main__":
    main()
