#!/usr/bin/env python3
"""
run_reconsensus_cohort.py
─────────────────────────
Parallel cohort launcher for the re-consensus rerun: partitions the remaining
(not-yet-succeeded) cohort samples into N disjoint groups and launches N
detached run_reconsensus_rerun.py driver processes (nohup-style), each running
its group sequentially with the cohort config (no VEP — see
config/rerun_cohort.yaml).

Concurrency-safety design:
  * Samples are partitioned disjointly, so per-sample artifacts
    (runs/rerun_csv/<sid>.csv, runs/rerun_checksums/<sid>.json,
    output_reconsensus/<sid>/) never collide.
  * The driver state file is NOT safe for concurrent writers (it is loaded
    once at startup and written back wholesale — last writer wins). Each
    group therefore gets its own state file via the driver's --state-file
    flag: runs/cohort_state/group<i>.json. Already-succeeded samples (from
    runs/rerun_state.json) are excluded before partitioning.
  * Each group driver runs with its own CWD (runs/cohort_work/group<i>/), so
    the nextflow work/ and .nextflow/ session dirs are per-group — concurrent
    nextflow instances must not share a launch dir.

Usage:
  # Show the partition and the exact commands; launch nothing
  python3 scripts/run_reconsensus_cohort.py --dry-run

  # Launch 6 detached group drivers (default groups: 6)
  python3 scripts/run_reconsensus_cohort.py

  # Re-launch (skips groups whose driver is still alive; completed samples
  # are auto-skipped inside each group by the driver itself)
  python3 scripts/run_reconsensus_cohort.py
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import run_reconsensus_rerun as rr  # noqa: E402

SEQ2NEO_ROOT = SCRIPT_DIR.parent
DEFAULT_CONFIG = SEQ2NEO_ROOT / "config" / "rerun_cohort.yaml"

STATE_DIR = SEQ2NEO_ROOT / "runs" / "cohort_state"
LOG_DIR = SEQ2NEO_ROOT / "runs" / "cohort_logs"
WORK_DIR = SEQ2NEO_ROOT / "runs" / "cohort_work"


def load_done_samples(state_path: Path) -> set:
    """sample_ids marked succeeded in the (sequential) driver state file."""
    if not state_path.exists():
        return set()
    state = json.loads(state_path.read_text())
    return {sid for sid, rec in state.items() if rec.get("status") == "succeeded"}


def partition(sample_ids: list, n_groups: int) -> list:
    """Round-robin partition (sorted ids → deterministic); balances sets,
    diseases, and expected runtimes across groups."""
    groups = [[] for _ in range(n_groups)]
    for i, sid in enumerate(sorted(sample_ids)):
        groups[i % n_groups].append(sid)
    return groups


def pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except (ProcessLookupError, OverflowError, ValueError):
        return False
    except PermissionError:
        return True
    return True


def main():
    ap = argparse.ArgumentParser(
        description="Parallel launcher for the re-consensus cohort rerun",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--config", default=str(DEFAULT_CONFIG),
                    help="Cohort config (default: config/rerun_cohort.yaml)")
    ap.add_argument("--groups", type=int, default=6,
                    help="Number of parallel group drivers (default: 6)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print the partition and commands; launch nothing")
    ap.add_argument("--force", action="store_true",
                    help="Launch even if a group's driver PID is still alive")
    args = ap.parse_args()

    cfg = rr.load_config(args.config, {})
    rows = rr.load_manifest(rr.resolve(cfg, "manifest_tsv"))
    done = load_done_samples(rr.resolve(cfg, "state_file"))
    remaining = [r["sample_id"] for r in rows if r["sample_id"] not in done]
    groups = partition(remaining, args.groups)

    print(f"manifest: {len(rows)} samples | already succeeded: {len(done)} "
          f"({', '.join(sorted(done)) or 'none'}) | to run: {len(remaining)}")
    print(f"config: {args.config}")
    print(f"groups: {args.groups} "
          f"(sizes: {', '.join(str(len(g)) for g in groups)})\n")

    launched = []
    for gi, sids in enumerate(groups, start=1):
        log_path = LOG_DIR / f"group{gi}.log"
        pid_path = LOG_DIR / f"group{gi}.pid"
        state_path = STATE_DIR / f"group{gi}.json"
        cwd = WORK_DIR / f"group{gi}"
        cmd = [
            sys.executable, str(SCRIPT_DIR / "run_reconsensus_rerun.py"),
            "--config", str(Path(args.config).resolve()),
            "--sample", ",".join(sids),
            "--state-file", str(state_path),
        ]
        print(f"group{gi}: {len(sids)} samples: {','.join(sids)}")
        print(f"  log   : {log_path}")
        print(f"  state : {state_path}")
        print(f"  cwd   : {cwd}  (per-group nextflow work/ + .nextflow/)")
        print(f"  cmd   : {' '.join(cmd)}")

        if pid_path.exists() and not args.force:
            pid = int(pid_path.read_text().strip())
            if pid_alive(pid):
                print(f"  SKIP  : driver already running (pid {pid}); "
                      f"use --force to override\n")
                continue

        if args.dry_run:
            print("  [dry-run] not launched\n")
            continue

        cwd.mkdir(parents=True, exist_ok=True)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        state_path.parent.mkdir(parents=True, exist_ok=True)
        log_fh = open(log_path, "a")
        proc = subprocess.Popen(
            cmd, stdout=log_fh, stderr=subprocess.STDOUT,
            cwd=str(cwd), start_new_session=True,  # detached (nohup-style)
        )
        pid_path.write_text(str(proc.pid) + "\n")
        launched.append((gi, proc.pid))
        print(f"  LAUNCHED pid {proc.pid}\n")

    if launched:
        print("Launched groups: "
              + ", ".join(f"group{gi}=pid{pid}" for gi, pid in launched))
        print(f"Monitor: tail -f {LOG_DIR}/group<N>.log ; "
              f"cat {STATE_DIR}/group<N>.json")


if __name__ == "__main__":
    main()
