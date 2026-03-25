# seq2neo — Usage Guide

## Split-root layout

Scripts self-locate using `$BASH_SOURCE` — no hardcoded repo path anywhere.
The only two things to configure before first use are in `config/runner.yaml`:

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

## Step 1 — Parse manifests

Edit `config/runner.yaml` first (see above), then:

```bash
# From anywhere — script self-locates
bash /your/path/rnadnavar/examples/seq2neo/parse.sh

# Or if you're already in the examples/seq2neo directory:
bash parse.sh
```

Writes to `$seq2neo_root/data/processed/`: `merged.json` + `set{1-4}_samples.tsv`.
Prints a validation report: status counts, set counts, disease exclusivity, incomplete list.

---

## Step 2 — Config

Edit `config/runner.yaml` — only three fields need to be set:

```yaml
main_nf:      /your/path/rnadnavar/main.nf
rdv_conf:     /your/path/rnadnavar/examples/seq2neo/seq2neo.shared.config
seq2neo_root: /your/data/output/root
```

Everything else (conda paths, proxy, lane, retry policy, completion artifacts) has defaults you can tune.

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

Each sample's CSV is written to `$seq2neo_root/runs/csv/<project>_<patient>.csv`:

```
patient,status,sample,lane,fastq_1,fastq_2
PRJNA298376_4060,0,PRJNA298376_4060DN,LX,/path/DN_1.fastq.gz,/path/DN_2.fastq.gz
PRJNA298376_4060,1,PRJNA298376_4060DT,LX,/path/DT_1.fastq.gz,/path/DT_2.fastq.gz
PRJNA298376_4060,2,PRJNA298376_4060RT,LX,/path/RT_1.fastq.gz,/path/RT_2.fastq.gz
```

To run manually:
```bash
REPO=/your/path/rnadnavar
DATA=/your/seq2neo_root

HTTPS_PROXY="http://10.233.17.241:3128" \
NXF_CONDA_CACHEDIR="/t9k/mnt/joey/nf_conda_envs" \
NXF_CONDA_USEMAMBA=true \
micromamba run -n nextflow nextflow run \
    $REPO/main.nf \
    -c $REPO/examples/seq2neo/seq2neo.shared.config \
    --input  $DATA/runs/csv/PRJNA298376_4060.csv \
    --outdir $DATA/output/PRJNA298376_4060 \
    -offline -with-conda -resume
```
