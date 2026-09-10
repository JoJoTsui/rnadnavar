#!/usr/bin/env python3
"""
run_reconsensus_rerun.py
────────────────────────
Re-consensus rerun driver: regenerates consensus + rescue VCFs for the seq2neo
cohort from EXISTING per-caller VCFs (Mutect2/Strelka2/DeepSomatic, DNA branch
and realigned RNA branch), writing to a NEW output location.

Hard properties (ticket 08 / audit finding M5):
  1. CONSENSUS+RESCUE ONLY — the generated samplesheet references caller VCFs
     only; nextflow is invoked with --step consensus and tools without any
     caller names, so FASTQ→BAM alignment and per-caller variant calling are
     structurally unreachable.
  2. NORMALIZATION GUARANTEED — the pipeline now routes caller VCFs through
     VCF_NORMALIZE (vt decompose + bcftools norm) when entering at the
     consensus step (see subworkflows/local/vcf_normalize/main.nf).
  3. READ-ONLY INPUTS + NEW OUTDIR — input VCFs are only opened for reading
     (md5); a checksum manifest recorded before each run is verified after it.
     Outputs land under a separate rerun outdir, never inside the source roots.

Sample list: data/processed/sample_manifest.tsv (build_sample_manifest.py).
Raw caller VCF layout per sample (base_output_dir/dir_name):
  variant_calling/{caller}/{p}DT_vs_{p}DN/{p}DT_vs_{p}DN<suffix>
  vcf_realignment/variant_calling/{caller}/{p}RT_realign_vs_{p}DN/...<suffix>
  suffix: mutect2=.mutect2.filtered.vcf.gz  strelka=.strelka.variants.vcf.gz
          deepsomatic=.deepsomatic.vcf.gz

Usage examples:
  # Dry-run the whole cohort (prints CSVs, commands, checksums; writes nothing)
  python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --dry-run

  # Run one sample
  python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml \\
      --sample PRJNA298330_4032

  # Run set 2
  python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --set 2

  # Re-verify recorded input checksums without running anything
  python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --verify-only

  # Resume (re-run same command — completed samples are auto-skipped)
  python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml

  # Parallel cohort run (6 detached group drivers, disjoint samples/state):
  python3 scripts/run_reconsensus_cohort.py --config config/rerun_cohort.yaml

  # One parallel group by hand: comma-separated --sample + a per-group state
  # file (the state file is loaded once and written back wholesale, so it is
  # NOT safe for concurrent processes sharing one file).
  python3 scripts/run_reconsensus_rerun.py --config config/rerun_cohort.yaml \\
      --sample ID1,ID2,ID3 --state-file runs/cohort_state/group1.json
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None  # only required when a YAML config file is actually loaded

# ---------------------------------------------------------------------------
# Caller VCF discovery
# ---------------------------------------------------------------------------

# Raw (pre-normalization) caller VCF file suffixes, matching what the original
# runs fed into VCF_NORMALIZE (see csv/variantcalled.csv in existing outdirs).
CALLER_VCF_SUFFIX = {
    "mutect2": ".mutect2.filtered.vcf.gz",  # FilterMutectCalls-filtered VCF
    "strelka": ".strelka.variants.vcf.gz",  # MergeVcfs SNV+indel VCF
    "deepsomatic": ".deepsomatic.vcf.gz",
}

# modality -> (status code, pair-id template, subdirectory under the sample dir)
MANIFEST_VCF_FIELDS = {
    "dna_mutect2": "caller_dna_mutect2",
    "dna_deepsomatic": "caller_dna_deepsomatic",
    "dna_strelka": "caller_dna_strelka",
    "rna_mutect2": "caller_rna_mutect2",
    "rna_deepsomatic": "caller_rna_deepsomatic",
    "rna_strelka": "caller_rna_strelka",
}

MODALITY_CONFIG = {
    "dna": {
        "status": 1,
        "pair": "{p}DT_vs_{p}DN",
        "subdir": "variant_calling",
    },
    "rna": {
        "status": 2,
        "pair": "{p}RT_realign_vs_{p}DN",
        "subdir": "vcf_realignment/variant_calling",
    },
}

CSV_HEADER = "patient,sample,status,variantcaller,vcf"


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
    "manifest_tsv": "data/processed/sample_manifest.tsv",
    "csv_dir": "runs/rerun_csv",
    "outdir_base": "output_reconsensus",
    "state_file": "runs/rerun_state.json",
    "checksum_dir": "runs/rerun_checksums",
    "step": "consensus",
    # No caller names here on purpose: with this tool set no variant-calling
    # process can be triggered (caller modules gate on their name in --tools).
    "tools": "consensus,rescue,filtering,vep",
    "native_evidence_snv": True,
    "rescue_promotion_enabled": True,
    "rescue_min_dna_callers": 1,
    "rescue_min_rna_callers": 2,
    "rescue_veto": "dna",
    "extra_nextflow_args": [],
    "dry_run": False,
    "resume": True,
    "offline": True,
    "max_retries": 2,
    # Generic completion artifacts: all must exist. Strict final-deliverable
    # check: the true label artifacts are the FILTERED VCFs (standard +
    # stripped) for the consensus branches and the rescue branch; the
    # intermediate *.rescued.vcf.gz alone is NOT completion (audit: stale-8).
    "completion_artifacts": [
        "**/pipeline_info/execution_trace*.txt",
        "filtered/**/*.filtered.vcf.gz",
        "filtered/**/*.filtered.vcf.stripped.vcf.gz",
        "rescue/**/*.filtered.vcf.gz",
        "rescue/**/*.filtered.vcf.stripped.vcf.gz",
    ],
    # Consensus artifacts: at least one must exist.
    "consensus_success_patterns": [
        "consensus/**/*.vcf.gz",
    ],
    # Rescue artifacts: at least one must exist for successful completion.
    # No *.rescued.vcf.gz fallback: that intermediate is produced mid-chain,
    # so its presence does not imply the run reached the label artifacts.
    "rescue_success_patterns": [
        "rescue/**/*.filtered.vcf.gz",
    ],
}


def load_config(config_path, cli_overrides: dict) -> dict:
    cfg = dict(DEFAULTS)
    if config_path and Path(config_path).exists():
        if yaml is None:
            sys.exit("PyYAML required for --config:  pip install pyyaml")
        with open(config_path) as f:
            file_cfg = yaml.safe_load(f) or {}
        cfg.update({k: v for k, v in file_cfg.items() if v is not None})
    cfg.update({k: v for k, v in cli_overrides.items() if v is not None})
    # Anchor relative config roots to the config location, not invocation CWD.
    # In the checked-in configs, `seq2neo_root: .` means the parent of config/.
    if config_path and cfg.get("seq2neo_root"):
        root = Path(str(cfg["seq2neo_root"]))
        if not root.is_absolute():
            config_dir = Path(config_path).resolve().parent
            cfg["seq2neo_root"] = str((config_dir.parent if str(root) == "." else config_dir / root).resolve())
    return cfg


def resolve(cfg: dict, key: str) -> Path:
    """Resolve path relative to seq2neo_root when not absolute."""
    p = Path(cfg[key])
    return p if p.is_absolute() else Path(cfg["seq2neo_root"]) / p


# ---------------------------------------------------------------------------
# Manifest loading and input discovery
# ---------------------------------------------------------------------------


def load_manifest(path: Path) -> list:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def locate_caller_vcfs(row: dict, manifest_dir: Path | None = None) -> dict:
    """Locate caller VCFs, preferring the manifest's audited explicit paths.

    The manifest is authoritative for seq2neo because its normalized VCF paths
    may not match the legacy raw-output layout. Legacy discovery remains a
    compatibility fallback for older manifests.
    """
    found = {}
    declared = set()
    for input_name, field in MANIFEST_VCF_FIELDS.items():
        value = str(row.get(field, "")).strip()
        if value:
            declared.add(input_name)
            candidate = Path(value)
            if not candidate.is_absolute() and manifest_dir is not None:
                candidate = manifest_dir / candidate
            if candidate.is_file():
                found[input_name] = candidate

    base = Path(row["base_output_dir"])
    if not base.is_absolute() and manifest_dir is not None:
        base = manifest_dir / base
    base = base / row["dir_name"]
    prefix = row["vcf_prefix"]
    for modality, mcfg in MODALITY_CONFIG.items():
        pair_id = mcfg["pair"].format(p=prefix)
        for caller, suffix in CALLER_VCF_SUFFIX.items():
            input_name = f"{modality}_{caller}"
            if input_name in found or input_name in declared:
                continue
            vcf = base / mcfg["subdir"] / caller / pair_id / f"{pair_id}{suffix}"
            if vcf.is_file():
                found[input_name] = vcf
    return found


def expected_inputs() -> list:
    return [
        f"{modality}_{caller}"
        for modality in MODALITY_CONFIG
        for caller in CALLER_VCF_SUFFIX
    ]


def patient_column(row: dict) -> str:
    """Value for the CSV patient column.

    nf-schema coerces a numeric-only CSV field to an integer, which fails
    the schema's patient-is-string validation. For cohort samples whose
    vcf_prefix is bare digits (e.g. "4278"), use sample_id instead
    (identical to vcf_prefix for all non-numeric prefixes).
    """
    if row["vcf_prefix"].isdigit():
        return row["sample_id"]
    return row["vcf_prefix"]


def write_sample_csv(row: dict, inputs: dict, csv_path: Path):
    """Generate a VCF-only samplesheet for the available modality panel."""
    patient = patient_column(row)
    prefix = row["vcf_prefix"]
    rows = [CSV_HEADER]
    for modality, mcfg in MODALITY_CONFIG.items():
        pair_id = mcfg["pair"].format(p=prefix)
        for caller in CALLER_VCF_SUFFIX:
            input_name = f"{modality}_{caller}"
            if input_name in inputs:
                rows.append(f"{patient},{pair_id},{mcfg['status']},{caller},{inputs[input_name]}")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.write_text("\n".join(rows) + "\n")


def render_sample_csv(row: dict, inputs: dict) -> str:
    """Render the VCF-only samplesheet without writing it."""
    patient = patient_column(row)
    prefix = row["vcf_prefix"]
    lines = [CSV_HEADER]
    for modality, mcfg in MODALITY_CONFIG.items():
        pair_id = mcfg["pair"].format(p=prefix)
        for caller in CALLER_VCF_SUFFIX:
            input_name = f"{modality}_{caller}"
            if input_name in inputs:
                lines.append(f"{patient},{pair_id},{mcfg['status']},{caller},{inputs[input_name]}")
    return "\n".join(lines) + "\n"


def classify_inputs(inputs: dict) -> tuple[str, list[str]]:
    """Return (READY|DNA_ONLY|PARTIAL, missing evidence names).

    Every caller VCF must have a tabix or CSI index before a run is started;
    otherwise the sample is marked partial rather than allowing a late
    Nextflow failure.
    """
    usable = set()
    missing_index = []
    for name, path in inputs.items():
        if any(Path(str(path) + ext).is_file() for ext in (".tbi", ".csi")):
            usable.add(name)
        else:
            missing_index.append(f"{name} (.tbi/.csi)")
    dna = {f"dna_{caller}" for caller in CALLER_VCF_SUFFIX}
    rna = {f"rna_{caller}" for caller in CALLER_VCF_SUFFIX}
    missing_dna = sorted(dna - usable)
    missing_rna = sorted(rna - usable)
    if missing_dna:
        return "PARTIAL", missing_dna + missing_rna
    if missing_rna:
        # DNA-only is valid only when the complete RNA panel is absent.
        # A partially present RNA panel is incomplete evidence and must not
        # silently skip rescue.
        all_rna = {f"rna_{caller}" for caller in CALLER_VCF_SUFFIX}
        if usable.isdisjoint(all_rna):
            return "DNA_ONLY", missing_rna
        return "PARTIAL", missing_rna
    if missing_index:
        return "PARTIAL", sorted(missing_index)
    return "READY", []


# ---------------------------------------------------------------------------
# Checksum guard (inputs are opened read-only, for hashing only)
# ---------------------------------------------------------------------------


def md5sum(path: Path, chunk_size: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def input_files(inputs: dict) -> list:
    """All files the run reads: caller VCFs plus their tabix indices if present."""
    files = []
    for vcf in inputs.values():
        files.append(vcf)
        for idx_ext in (".tbi", ".csi"):
            idx = Path(str(vcf) + idx_ext)
            if idx.is_file():
                files.append(idx)
    return files


def record_checksums(inputs: dict, manifest_path: Path) -> dict:
    """Record md5 + size for every input file. Returns the manifest dict."""
    entries = {}
    for f in sorted(input_files(inputs)):
        entries[str(f)] = {"md5": md5sum(f), "size": f.stat().st_size}
    manifest = {
        "recorded": datetime.now().isoformat(),
        "files": entries,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2))
    return manifest


def verify_checksums(manifest_path: Path) -> tuple:
    """
    Re-hash every file in a recorded manifest.
    Returns (ok, problems) where problems is a list of human-readable strings.
    """
    if not manifest_path.exists():
        return False, [f"checksum manifest missing: {manifest_path}"]
    manifest = json.loads(manifest_path.read_text())
    problems = []
    for path_str, expected in manifest["files"].items():
        f = Path(path_str)
        if not f.is_file():
            problems.append(f"missing input: {path_str}")
            continue
        actual = md5sum(f)
        if actual != expected["md5"]:
            problems.append(f"checksum mismatch: {path_str}")
    return (not problems), problems


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

# Standard BGZF end-of-file marker block (28 bytes). Published VCFs are
# bgzf-compressed (all carry .tbi); a missing EOF block means the writer was
# killed mid-stream (the cohort's dominant failure mode was OOM truncation).
BGZF_EOF_BLOCK = bytes.fromhex(
    "1f8b08040000000000ff" "0600" "4243" "0200" "1b00" "0300"
    "00000000" "00000000"
)


def vcf_gz_sane(path: Path, min_records: int = 1) -> bool:
    """Cheap integrity check for a bgzf VCF: plausible size, BGZF EOF marker,
    readable header (#CHROM) and at least min_records data lines."""
    try:
        if path.stat().st_size < 1024:
            return False
        with open(path, "rb") as fh:
            fh.seek(-28, os.SEEK_END)
            if fh.read(28) != BGZF_EOF_BLOCK:
                return False
        saw_header = False
        records = 0
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    saw_header = True
                elif not line.startswith("#"):
                    records += 1
                    if saw_header and records >= min_records:
                        return True
        return saw_header and records >= min_records
    except (OSError, EOFError, UnicodeDecodeError):
        return False


def has_any_match(outdir: Path, patterns: list) -> bool:
    """outdir must exist AND at least one glob pattern must match ≥1 file."""
    if not outdir.exists():
        return False
    return any(list(outdir.glob(p)) for p in patterns)


def has_all_matches(outdir: Path, patterns: list) -> bool:
    """outdir must exist AND every glob pattern must match ≥1 file."""
    if not outdir.exists():
        return False
    return all(list(outdir.glob(p)) for p in patterns)


def evaluate_completion(outdir: Path, cfg: dict) -> tuple:
    """Completion = generic artifacts + consensus VCF(s) + rescue VCF(s),
    with an integrity check on every completion-artifact VCF so that
    truncated/corrupt outputs from killed runs are not waved through."""
    completion_patterns = list(cfg["completion_artifacts"])
    if not cfg.get("require_rescue", True):
        completion_patterns = [p for p in completion_patterns if not p.startswith("rescue/")]
    if not has_all_matches(outdir, completion_patterns):
        return False, "completion artifacts missing"
    if not has_any_match(outdir, cfg["consensus_success_patterns"]):
        return False, "consensus artifacts missing"
    if cfg.get("require_rescue", True) and not has_any_match(outdir, cfg["rescue_success_patterns"]):
        return False, "rescue artifacts missing"
    for pattern in completion_patterns:
        for f in outdir.glob(pattern):
            if f.name.endswith(".vcf.gz") and not vcf_gz_sane(f):
                return False, f"corrupt or truncated VCF: {f}"
    return True, ""


# ---------------------------------------------------------------------------
# Nextflow command builder
# ---------------------------------------------------------------------------


def config_path(cfg: dict, key: str) -> str:
    value = Path(str(cfg[key]))
    if value.is_absolute():
        return str(value)
    return str((Path(cfg["seq2neo_root"]) / value).resolve())


def build_command(cfg: dict, csv_path: Path, outdir: Path) -> list:
    cmd = ["micromamba", "run", "-n", cfg["micromamba_env"]]
    cmd += ["nextflow", "run"]
    if cfg.get("main_nf"):
        cmd.append(config_path(cfg, "main_nf"))
    if cfg.get("rdv_conf"):
        cmd += ["-c", config_path(cfg, "rdv_conf")]
    cmd += ["--input", str(csv_path)]
    cmd += ["--outdir", str(outdir)]
    cmd += ["--step", cfg["step"]]
    cmd += ["--tools", cfg["tools"]]
    # Policy flags are explicit in the generated CLI so dry-runs are auditable.
    cmd += ["--native_evidence_snv", str(bool(cfg.get("native_evidence_snv", True))).lower()]
    cmd += ["--rescue_promotion_enabled", str(bool(cfg.get("rescue_promotion_enabled", True))).lower()]
    cmd += ["--rescue_min_dna_callers", str(cfg.get("rescue_min_dna_callers", 1))]
    cmd += ["--rescue_min_rna_callers", str(cfg.get("rescue_min_rna_callers", 2))]
    cmd += ["--rescue_veto", str(cfg.get("rescue_veto", "dna"))]
    for extra in cfg.get("extra_nextflow_args") or []:
        cmd.append(str(extra))
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
            f"Edit config/rerun.yaml and set:\n"
            + "\n".join(f"  {k}: /path/to/..." for k in missing)
        )
    # Structural guard: rerun tools are an explicit VCF-only allowlist.
    allowed = {"consensus", "rescue", "filtering", "vep"}
    requested = {t.strip() for t in str(cfg.get("tools", "")).split(",") if t.strip()}
    unsupported = requested - allowed
    if unsupported:
        sys.exit(
            f"ERROR: rerun tools outside VCF-only allowlist: {sorted(unsupported)}.\n"
            "Allowed tools are consensus,rescue,filtering,vep."
        )
    callers = {"sage", "strelka", "mutect2", "deepsomatic", "manta"}
    named = callers & requested
    if named:
        sys.exit(
            f"ERROR: rerun tools must not name variant callers (found: {sorted(named)}).\n"
            "The rerun path is consensus+rescue only; naming a caller would make\n"
            "variant-calling processes reachable."
        )
    if cfg.get("step") != "consensus":
        sys.exit(
            f"ERROR: rerun step must be 'consensus' (got: {cfg.get('step')}).\n"
            "Other steps can trigger alignment/preprocessing/variant calling."
        )


def check_outdir_separation(rows: list, outdir_base: Path):
    """Refuse to run if the rerun outdir lives inside any source output root."""
    outdir_resolved = outdir_base.resolve()
    source_roots = {str(Path(r["base_output_dir"]).resolve()) for r in rows}
    for root in source_roots:
        if outdir_resolved == Path(root) or str(outdir_resolved).startswith(root + os.sep):
            sys.exit(
                f"ERROR: rerun outdir {outdir_resolved} is inside source output root {root}.\n"
                "Rerun outputs must go to a new, separate location."
            )


# ---------------------------------------------------------------------------
# Sample filtering
# ---------------------------------------------------------------------------


def filter_samples(rows: list, args) -> list:
    # --sample accepts a comma-separated list of sample_ids (used by the
    # parallel cohort launcher to hand one disjoint group to each process).
    wanted = set(args.sample.split(",")) if args.sample else None
    out = []
    for r in rows:
        if wanted and r["sample_id"] not in wanted:
            continue
        if args.set is not None and str(r.get("set_number")) != str(args.set):
            continue
        if args.status and args.status != "all" and r.get("status") != args.status:
            continue
        out.append(r)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(
        description="Consensus+rescue-only rerun driver over existing caller VCFs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--config", default=None, help="Path to rerun.yaml")
    # filters
    ap.add_argument("--sample", default=None,
                    help="Run one sample_id, or a comma-separated list of sample_ids")
    ap.add_argument("--set", default=None, type=int, help="Filter by partition set (1-4)")
    ap.add_argument(
        "--status", default="all", help="standard | extra | all  (default: all)"
    )
    # modes
    ap.add_argument(
        "--verify-only",
        action="store_true",
        help="Only verify recorded input checksums; do not run anything",
    )
    # config overrides
    ap.add_argument("--main-nf", default=None, dest="main_nf")
    ap.add_argument("--rdv-conf", default=None, dest="rdv_conf")
    ap.add_argument("--outdir-base", default=None, dest="outdir_base")
    ap.add_argument("--state-file", default=None, dest="state_file",
                    help="Override state_file (per-group state for parallel runs; "
                         "the state file is not safe for concurrent writers)")
    ap.add_argument("--seq2neo", default=None, dest="seq2neo_root")
    ap.add_argument("--manifest", default=None, dest="manifest_tsv")
    ap.add_argument("--dry-run", action="store_true", default=None, dest="dry_run")
    ap.add_argument("--no-resume", action="store_false", default=None, dest="resume")
    args = ap.parse_args()

    filter_keys = {"sample", "set", "status", "verify_only", "config"}
    cli_overrides = {k: v for k, v in vars(args).items() if k not in filter_keys}
    cfg = load_config(args.config, cli_overrides)

    manifest_path = resolve(cfg, "manifest_tsv")
    if not manifest_path.exists():
        sys.exit(
            f"sample manifest not found: {manifest_path}\n"
            f"Run build_sample_manifest.py first."
        )
    rows = load_manifest(manifest_path)

    selected = filter_samples(rows, args)
    if not selected:
        print("No samples match the given filters.")
        return

    state_path = resolve(cfg, "state_file")
    checksum_dir = resolve(cfg, "checksum_dir")
    csv_dir = resolve(cfg, "csv_dir")
    outdir_base = resolve(cfg, "outdir_base")
    dry_run = bool(cfg.get("dry_run"))

    # ── verify-only mode ────────────────────────────────────────────────────
    if args.verify_only:
        n_ok = n_bad = n_missing = 0
        for r in selected:
            sid = r["sample_id"]
            manifest_file = checksum_dir / f"{sid}.input_checksums.json"
            if not manifest_file.exists():
                print(f"[NONE]  {sid}  no checksum manifest recorded")
                n_missing += 1
                continue
            ok, problems = verify_checksums(manifest_file)
            if ok:
                print(f"[OK]    {sid}  inputs untouched")
                n_ok += 1
            else:
                print(f"[FAIL]  {sid}  {'; '.join(problems)}")
                n_bad += 1
        print(
            f"\n=== Verify-only  ok={n_ok}  failed={n_bad}  no-manifest={n_missing} ==="
        )
        sys.exit(1 if n_bad else 0)

    if not dry_run:
        validate_config(cfg)
        check_outdir_separation(selected, outdir_base)

    state = load_state(state_path)

    print(f"\n{'[DRY RUN] ' if dry_run else ''}Selected {len(selected)} sample(s)\n")

    counts = {"succeeded": 0, "skipped": 0, "failed": 0, "dna_only": 0, "partial": 0}

    for r in selected:
        sid = r["sample_id"]
        key = sid
        outdir = outdir_base / sid
        csv_path = csv_dir / f"{sid}.csv"
        checksum_manifest = checksum_dir / f"{sid}.input_checksums.json"

        # ── locate inputs and classify before completion checks ─────────────
        inputs = locate_caller_vcfs(r, manifest_path.parent)
        input_status, missing = classify_inputs(inputs)
        completion_cfg = dict(cfg, require_rescue=input_status == "READY")
        complete_ok, complete_reason = evaluate_completion(outdir, completion_cfg)
        if complete_ok:
            print(f"[SKIP]  {key}  already complete ({input_status})")
            state[key] = {"status": "succeeded", "reason": "artifacts present", "input_status": input_status}
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
        if input_status == "PARTIAL":
            print(f"[PARTIAL]  {key}  missing caller VCFs: {missing}")
            if not dry_run:
                state[key] = {"status": "partial", "reason": f"missing caller VCFs: {missing}"}
                save_state(state_path, state)
            counts["partial"] += 1
            continue
        run_cfg = dict(cfg)
        if input_status == "DNA_ONLY":
            print(f"[DNA_ONLY] {key}  RNA caller VCFs unavailable; rescue will be skipped")
            counts["dna_only"] += 1
            run_cfg["tools"] = ",".join(t for t in str(cfg["tools"]).split(",") if t != "rescue")

        cmd = build_command(run_cfg, csv_path, outdir)

        print(f"[{'DRY' if dry_run else 'RUN'}]  {key}  set{r.get('set_number')} | {r.get('disease')}")
        print(f"       csv : {csv_path}")
        print(f"       out : {outdir}")
        print(f"       cmd : {format_shell_command(run_cfg, cmd)}")
        if dry_run:
            print("       --- samplesheet ---")
            for line in render_sample_csv(r, inputs).splitlines():
                print(f"         {line}")
            print("       --- input md5 (read-only) ---")
            for f in sorted(input_files(inputs)):
                print(f"         {md5sum(f)}  {f}")
            continue

        # ── generate CSV + record input checksums (before) ──────────────────
        write_sample_csv(r, inputs, csv_path)
        record_checksums(inputs, checksum_manifest)

        # ── execute ─────────────────────────────────────────────────────────
        state[key] = {
            "status": "running",
            "started": datetime.now().isoformat(),
            "retries": retries,
        }
        save_state(state_path, state)

        try:
            subprocess.run(cmd, env=build_env(cfg), check=True)
            complete_ok, complete_reason = evaluate_completion(outdir, completion_cfg)
            checksum_ok, checksum_problems = verify_checksums(checksum_manifest)
            if not checksum_ok:
                state[key] = {
                    "status": "failed",
                    "reason": "input checksums changed during run: "
                    + "; ".join(checksum_problems),
                    "retries": cfg["max_retries"],  # never auto-retry a guard trip
                }
                print(f"[FAIL]  {key}  checksum guard tripped: {checksum_problems}")
                counts["failed"] += 1
            elif complete_ok:
                state[key] = {
                    "status": "succeeded",
                    "finished": datetime.now().isoformat(),
                    "inputs_verified": True,
                }
                print(f"[OK]    {key}  (input checksums verified unchanged)")
                counts["succeeded"] += 1
            else:
                state[key] = {
                    "status": "failed",
                    "reason": complete_reason or "artifacts missing after run",
                    "retries": retries + 1,
                }
                print(f"[FAIL]  {key}  {complete_reason}")
                counts["failed"] += 1
        except subprocess.CalledProcessError as e:
            state[key] = {"status": "failed", "reason": str(e), "retries": retries + 1}
            print(f"[FAIL]  {key}  {e}")
            counts["failed"] += 1

        save_state(state_path, state)

    print(
        f"\n=== Done  succeeded={counts['succeeded']} skipped={counts['skipped']} "
        f"dna_only={counts['dna_only']} partial={counts['partial']} failed={counts['failed']} ==="
    )


if __name__ == "__main__":
    main()
