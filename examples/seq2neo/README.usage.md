# seq2neo — Usage Guide

## Split-root layout

Scripts self-locate using `$BASH_SOURCE` — no hardcoded repo path anywhere.
Before first use, set three values in `config/runner.yaml`:

| Key | What to set |
|-----|-------------|
| `main_nf` | absolute path to `rnadnavar/main.nf` on this machine |
| `rdv_conf` | absolute path to `examples/seq2neo/seq2neo.shared.config` on this machine |
| `seq2neo_root` | where data/output/state will be written |

The repo can live at any path. Scripts always resolve themselves:
```bash
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
```

### Repo tree (versioned, read-only at runtime)

```
$REPO/examples/seq2neo/
├── PRJNA298376.txt                  # raw manifests
├── PRJNA298330.txt
├── PRJNA298310.txt
├── PRJNA298330.disease.tsv          # per-patient disease annotation
├── seq2neo.shared.config            # nextflow pipeline config
├── config/
│   └── runner.yaml                  # execution config (seq2neo_root → DATA)
├── scripts/
│   ├── lib/
│   │   ├── __init__.py
│   │   └── common.py                # shared logic (disease rules, classification, partitioning)
│   ├── parse_projects_to_json.py    # step 1: parse manifests → merged JSON
│   └── run_batch_from_json.py       # step 2: batch nextflow runner
├── parse.sh                         # convenience: run parse step
├── dry_run.sh                       # convenience: preview any filter
├── run_set1.sh                      # Colorectal cancer
├── run_set2.sh                      # Colon cancer
├── run_set3.sh                      # Ampullary / Bile Duct / Cholangiocarcinoma / Esophageal / Melanoma
├── run_set4.sh                      # Gastric / Lung / Pancreatic / Rectal
└── run_single.sh                    # single patient
```

### Workspace tree (generated at runtime, not in repo)

```
$DATA/
├── data/processed/
│   ├── merged.json                  # unified sample database
│   └── set{1-4}_samples.tsv         # partition membership lists
├── runs/
│   ├── csv/                         # per-sample nextflow input CSVs
│   └── run_state.json               # resume state
└── output/
    └── PRJNA298376_4060/            # per-sample nextflow output
        └── ...
```

---

## Sample status

| Status     | Meaning                                                   | Eligible? |
|------------|-----------------------------------------------------------|-----------|
| standard   | Exactly 1 DN + 1 DT + 1 RT pair                          | Yes       |
| extra      | All 3 modalities present, ≥1 has >1 pair                  | Yes — first pair used |
| incomplete | Missing DN, DT, or RT entirely                            | No        |

---

## Partition sets

| Set  | Disease(s)                                                                    |
|------|-------------------------------------------------------------------------------|
| set1 | Colorectal cancer (strict — colon/rectal are NOT included here)               |
| set2 | Colon cancer                                                                  |
| set3 | Ampullary / Bile Duct / Cholangiocarcinoma / Esophageal / Melanoma            |
| set4 | Gastric / Lung / Pancreatic / Rectal                                          |

Sets 2–4 are disease-exclusive (each disease appears in exactly one set).
Size balance is a secondary objective.

---

## Configuration reference

### config/runner.yaml — all options

```yaml
# ── Nextflow pipeline paths ───────────────────────────────────────────────
main_nf:  /path/to/rnadnavar/main.nf
#   Absolute path to the rnadnavar pipeline entry point.
#   To change the nextflow version or repo location, update this path.

rdv_conf: /path/to/rnadnavar/examples/seq2neo/seq2neo.shared.config
#   Absolute path to the shared nextflow config for this project.
#   Contains all pipeline params (genome, tools, reference files, etc.).
#   Can be overridden at runtime with --rdv-conf.

# ── Conda / Nextflow environment ──────────────────────────────────────────
nxf_conda_cachedir: /path/to/nf_conda_envs
#   Where nextflow stores/reuses conda environments (NXF_CONDA_CACHEDIR).
#   Set this to a shared, persistent directory to avoid re-creating envs.
#   This is the nextflow working/cache directory for conda environments.

nxf_conda_usemamba: "true"
#   Use mamba instead of conda for faster env creation (NXF_CONDA_USEMAMBA).

micromamba_env: nextflow
#   The micromamba environment that has nextflow installed.
#   Used as: micromamba run -n <micromamba_env> nextflow run ...

# ── Network proxy ─────────────────────────────────────────────────────────
https_proxy: "http://10.233.17.241:3128"
#   HTTP/HTTPS proxy for nextflow and conda downloads.
#   Leave empty ("") to disable.

# ── Data root ─────────────────────────────────────────────────────────────
seq2neo_root: /path/to/seq2neo_data
#   Root directory for all runtime data. All relative paths below are
#   resolved relative to this directory.

merged_json: data/processed/merged.json
#   Input: parsed sample database produced by parse.sh (step 1).

csv_dir: runs/csv
#   Where per-sample nextflow input CSVs are written before each run.

outdir_base: output
#   Base directory for per-sample nextflow output.
#   Each sample writes to: <outdir_base>/<project>_<patient>/

state_file: runs/run_state.json
#   JSON file tracking run status (succeeded/failed/running) per sample.
#   Used for resume logic — do not delete unless you want to reset all state.

# ── Execution options ─────────────────────────────────────────────────────
lane: LX
#   Lane identifier written into the nextflow input CSV.

dry_run: false
#   If true, print commands without executing. Same as --dry-run flag.

resume: true
#   Pass -resume to nextflow, enabling Nextflow's built-in task caching.

offline: true
#   Pass -offline to nextflow, preventing remote resource fetching.

max_parallel: 1
#   Number of samples to run concurrently. 1 = sequential (safe default).

# ── Retry policy ──────────────────────────────────────────────────────────
max_retries: 1
#   How many times a failed sample is retried before being permanently skipped.

# ── Completion check ──────────────────────────────────────────────────────
completion_artifacts:
  - "**/*.filtered.vcf.gz"
  - "**/pipeline_info/execution_trace*.txt"
#   All globs must match ≥1 file under outdir/<project>_<patient>/ for a
#   sample to be considered successfully finished. Adjust if your pipeline
#   produces different output files.
```

### seq2neo.shared.config — pipeline parameters

This file is the nextflow config passed via `-c`. It sets all pipeline-level
parameters. Key sections:

```
Resource limits     — cpus, memory, time per process
Reference genome    — fasta, fasta_fai, genome (GRCh38)
Aligners            — bwa index, star_index, hisat2_index, aligner
Annotation DBs      — dbsnp, pon, germline_resource, vep_cache, gtf
RNA editing         — rediportal_vcf, min_rna_support
COSMIC / gnomAD     — cosmic_database, gnomad_database, thresholds
Pipeline tools      — tools = deepsomatic,mutect2,strelka,vep,norm,...
Step                — step = mapping (start from raw FASTQ)
```

To change reference paths or enable/disable tools, edit this file directly.
To use a different config file entirely, update `rdv_conf` in `runner.yaml`.

---

## How to change the nextflow path

The nextflow executable is invoked via micromamba:

```
micromamba run -n <micromamba_env> nextflow run <main_nf> ...
```

There are two independent paths to configure:

**1. The pipeline script (`main_nf`)** — which `main.nf` to run:
```yaml
# config/runner.yaml
main_nf: /new/path/to/rnadnavar/main.nf
```
Or override at runtime without editing the file:
```bash
bash run_set2.sh --main-nf /new/path/to/rnadnavar/main.nf
```

**2. The nextflow executable** — which nextflow installation to use:
```yaml
# config/runner.yaml
micromamba_env: nextflow   # micromamba environment containing nextflow
```
If nextflow lives in a different conda/micromamba environment, change `micromamba_env`.
If you use a different launcher (e.g., plain `nextflow` on PATH), you would need to
edit `build_command()` in `scripts/run_batch_from_json.py`.

---

## How to change the nextflow working directory (cache dir)

Nextflow uses two kinds of cache/work directories:

**1. Conda environment cache (`nxf_conda_cachedir`)** — where conda envs are stored:
```yaml
# config/runner.yaml
nxf_conda_cachedir: /shared/path/nf_conda_envs
```
This maps to the `NXF_CONDA_CACHEDIR` environment variable. Set it to a
persistent, shared location to avoid re-creating environments on every run.

**2. Nextflow work directory** — where task intermediate files are cached.
Nextflow defaults this to `./work` relative to where the command is run.
To change it, add `-work-dir` to the nextflow invocation in `seq2neo.shared.config`
or pass it via the pipeline params. Alternatively, set `NXF_WORK` in your environment:
```bash
export NXF_WORK=/scratch/nextflow_work
bash run_set2.sh
```
Or set it in `runner.yaml` under a custom env block if you extend the runner script.

---

## Step 1 — Parse manifests

Edit `config/runner.yaml` first (see above), then:

```bash
# From anywhere — script self-locates
bash /your/path/rnadnavar/examples/seq2neo/parse.sh

# Or if you're already in the examples/seq2neo directory:
bash parse.sh
```

Reads manifest files from `seq2neo_root` (configured in `runner.yaml`):
- `PRJNA298376.txt` — tree-style manifest
- `PRJNA298330.txt` — pipe-style manifest
- `PRJNA298310.txt` — pipe-style manifest (all Melanoma)
- `PRJNA298330.disease.tsv` — per-patient disease annotation

Writes to `$seq2neo_root/data/processed/`:
- `merged.json` — unified sample database with modality paths, status, partition set
- `set1_samples.tsv` through `set4_samples.tsv` — partition membership lists

Prints a validation report: status counts, set counts, disease exclusivity, incomplete list.

---

## Step 2 — Config

Edit `config/runner.yaml` — only three fields need to be set before first use:

```yaml
main_nf:      /your/path/rnadnavar/main.nf
rdv_conf:     /your/path/rnadnavar/examples/seq2neo/seq2neo.shared.config
seq2neo_root: /your/data/output/root
```

Everything else has sensible defaults. See the full configuration reference above.

---

## Step 3 — Dry-run (preview, no execution)

```bash
SEQ2NEO_SCRIPTS=/your/path/rnadnavar/examples/seq2neo

# All eligible samples
bash $SEQ2NEO_SCRIPTS/dry_run.sh

# Set1 only (Colorectal cancer)
bash $SEQ2NEO_SCRIPTS/dry_run.sh --set 1

# Set1, standard samples only
bash $SEQ2NEO_SCRIPTS/dry_run.sh --set 1 --status standard

# One specific patient
bash $SEQ2NEO_SCRIPTS/dry_run.sh --project PRJNA298376 --patient 4060

# Filter by disease substring
bash $SEQ2NEO_SCRIPTS/dry_run.sh --disease "Pancreatic"
```

Dry-run prints the exact nextflow command that would be executed for each sample,
including the CSV path, output directory, and all environment variables.
Nothing is written to disk and no nextflow process is started.

---

## Step 4 — Run

```bash
SEQ2NEO_SCRIPTS=/your/path/rnadnavar/examples/seq2neo

# Run by set
bash $SEQ2NEO_SCRIPTS/run_set1.sh   # Colorectal cancer
bash $SEQ2NEO_SCRIPTS/run_set2.sh   # Colon cancer
bash $SEQ2NEO_SCRIPTS/run_set3.sh   # Ampullary/Bile Duct/Cholangiocarcinoma/Esophageal/Melanoma
bash $SEQ2NEO_SCRIPTS/run_set4.sh   # Gastric/Lung/Pancreatic/Rectal

# Run a single patient
PROJECT=PRJNA298376 PATIENT=4060 bash $SEQ2NEO_SCRIPTS/run_single.sh

# CLI override without editing runner.yaml
bash $SEQ2NEO_SCRIPTS/run_set2.sh \
    --main-nf /other/path/main.nf \
    --outdir-base /scratch/output
```

### CLI flags (all run scripts)

All run scripts pass extra arguments through to `run_batch_from_json.py`.
Available overrides:

| Flag | Overrides | Description |
|------|-----------|-------------|
| `--main-nf PATH` | `main_nf` | Path to `main.nf` |
| `--rdv-conf PATH` | `rdv_conf` | Path to nextflow config |
| `--outdir-base PATH` | `outdir_base` | Base output directory |
| `--seq2neo PATH` | `seq2neo_root` | Data root directory |
| `--dry-run` | `dry_run` | Preview without executing |
| `--no-resume` | `resume` | Disable `-resume` flag |
| `--project ID` | — | Filter by project ID |
| `--patient ID` | — | Filter by patient ID |
| `--disease STR` | — | Filter by disease (substring) |
| `--set N` | — | Filter by partition set (1–4) |
| `--status STR` | — | Filter by status: `standard`, `extra`, or `all` |

---

## Outputs

### Per-sample nextflow output

Each sample writes to `$seq2neo_root/output/<project>_<patient>/`, e.g.:
```
output/PRJNA298376_4060/
├── pipeline_info/
│   ├── execution_trace_*.txt        # nextflow task-level timing and resource usage
│   ├── execution_report_*.html      # visual execution report
│   └── execution_timeline_*.html    # timeline view
├── <sample>/
│   ├── *.filtered.vcf.gz            # final filtered variant calls (completion marker)
│   ├── *.filtered.vcf.gz.tbi
│   └── ...                          # BAMs, intermediate VCFs, QC files
```

A sample is considered complete when ALL `completion_artifacts` globs match
at least one file under its output directory (configurable in `runner.yaml`).

### Run state

`$seq2neo_root/runs/run_state.json` — tracks status per sample:
```json
{
  "PRJNA298376::4060": {
    "status": "succeeded",
    "finished": "2024-01-15T10:23:45"
  }
}
```
Possible status values: `running`, `succeeded`, `failed`.

### Per-sample input CSVs

`$seq2neo_root/runs/csv/<project>_<patient>.csv` — generated before each run:
```
patient,status,sample,lane,fastq_1,fastq_2
PRJNA298376_4060,0,PRJNA298376_4060DN,LX,/path/DN_1.fastq.gz,/path/DN_2.fastq.gz
PRJNA298376_4060,1,PRJNA298376_4060DT,LX,/path/DT_1.fastq.gz,/path/DT_2.fastq.gz
PRJNA298376_4060,2,PRJNA298376_4060RT,LX,/path/RT_1.fastq.gz,/path/RT_2.fastq.gz
```
Status codes: DN=0 (normal DNA), DT=1 (tumor DNA), RT=2 (tumor RNA).

### Parse outputs

`$seq2neo_root/data/processed/`:
- `merged.json` — full sample database with modality paths, status, partition assignment
- `set1_samples.tsv` through `set4_samples.tsv` — TSV listing project, patient, disease, status per set

---

## Resume / continue

State is tracked in `$seq2neo_root/runs/run_state.json`. Re-running the same command is safe — finished samples are auto-skipped.

**How completion is detected (strict):**
- `$seq2neo_root/output/<project>_<patient>/` must exist
- All globs in `completion_artifacts` (runner.yaml) must match ≥1 file

**Resume after interruption:**
```bash
# Just re-run — already-finished samples are skipped automatically
bash $SEQ2NEO_SCRIPTS/run_set2.sh
```

**Force re-run one sample:**
```bash
DATA=/your/seq2neo_root

python3 -c "
import json
p = '$DATA/runs/run_state.json'
s = json.load(open(p))
s.pop('PRJNA298376::4060', None)
json.dump(s, open(p, 'w'), indent=2)
"
rm -rf $DATA/output/PRJNA298376_4060
PROJECT=PRJNA298376 PATIENT=4060 bash $SEQ2NEO_SCRIPTS/run_single.sh
```

**Reset all failed samples for retry:**
```bash
python3 -c "
import json
p = '$DATA/runs/run_state.json'
s = json.load(open(p))
for v in s.values():
    if v.get('status') == 'failed':
        v['retries'] = 0
json.dump(s, open(p, 'w'), indent=2)
"
bash $SEQ2NEO_SCRIPTS/run_set2.sh
```

---

## Manual nextflow invocation (single sample)

Each sample's CSV is written to `$seq2neo_root/runs/csv/<project>_<patient>.csv`.
To run manually:

```bash
REPO=/your/path/rnadnavar
DATA=/your/seq2neo_root

HTTPS_PROXY="http://10.233.17.241:3128" \
NXF_CONDA_CACHEDIR="/path/to/nf_conda_envs" \
NXF_CONDA_USEMAMBA=true \
micromamba run -n nextflow nextflow run \
    $REPO/main.nf \
    -c $REPO/examples/seq2neo/seq2neo.shared.config \
    --input  $DATA/runs/csv/PRJNA298376_4060.csv \
    --outdir $DATA/output/PRJNA298376_4060 \
    -offline -with-conda -resume
```

To change the nextflow work directory for this invocation:
```bash
NXF_WORK=/scratch/nextflow_work \
micromamba run -n nextflow nextflow run ...
```
