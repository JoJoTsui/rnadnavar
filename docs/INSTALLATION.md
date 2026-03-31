# rnadnavar — Installation Guide

## Overview

This guide covers a full installation of the rnadnavar pipeline on a Linux server
(Ubuntu 22.04) without Docker, using conda/micromamba for environment management.

Three components need to be installed:

1. **Nextflow** — the workflow engine
2. **rnadnavar pipeline** — the pipeline code and conda environments
3. **DeepSomatic** — local binary installation (no Docker required)

The shared storage (reference databases, pipeline code, test data, pre-built
DeepSomatic binaries) is accessible to all teammates. However, **each user must
install their own runtime environment** — Nextflow, micromamba, and conda
environments are per-user and not shared.

---

## Shared infrastructure (shared storage — accessible to all)

| Resource | Path |
|----------|------|
| Pipeline repo | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/` |
| Reference databases | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/` |
| Test dataset (COO8801) | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/` |
| DeepSomatic binaries (rsync source) | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/DeepVariant-1.9.0/` |
| DeepSomatic models | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/models/` |

> **Note:** Conda environments and Nextflow runtime are **not** shared. Each user
> must install their own. See [NXF_CONDA_CACHEDIR](#nxf_conda_cachedir) below.

---

## Caveats

### Nextflow invocation

The pipeline scripts use `micromamba run -n nextflow nextflow run ...` to invoke
Nextflow. This wrapper is only needed when Nextflow is installed inside a
micromamba/conda environment. If Nextflow is already on your `PATH` (e.g.,
installed system-wide or via the official installer), call it directly:

```bash
nextflow run /path/to/main.nf -c /path/to/config ...
```

### micromamba vs conda

`micromamba` is a fast, standalone drop-in replacement for `conda`. All commands
in this guide use `micromamba`. If you only have `conda`, replace `micromamba`
with `conda` — the syntax is identical.

### Output directory

Default `OUTDIR` values in the example scripts are absolute paths under the shared
`rnadnavar_test` directory. For your own runs, always use an absolute path outside
the repo to avoid writing output into the codebase.

### Offline mode

The pipeline runs with `-offline` to prevent Nextflow from fetching remote
resources. The first run on a new machine requires internet access (or proxy) to
populate `NXF_CONDA_CACHEDIR` with conda environments. Subsequent runs work
offline.

### NXF_CONDA_CACHEDIR

Set this to a directory **in your own user space** where Nextflow will build and
cache conda environments. Do not point this at a shared location — conda
environments are per-user and must be built locally.

```bash
# Example: use a directory under your home or personal workspace
export NXF_CONDA_CACHEDIR="${HOME}/nf_conda_envs"
# Or a dedicated workspace path:
export NXF_CONDA_CACHEDIR="/t9k/mnt/hdd/work/${USER}/nf_conda_envs"
```

---

## Installation steps

Every user must complete all steps below. There is no shared environment to join.

### Requirements
- Python 3.10.12
- CPUs with SSE4 and AVX support
- `sudo` privilege
- Internet access (or proxy configured)
- Minimum 36 CPUs, 160 GB RAM recommended for the full pipeline

---

## Step 1 — Install micromamba

```bash
# Install micromamba (fast conda replacement)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# Or via conda-forge if conda is already available
conda install -c conda-forge micromamba

# Verify
micromamba --version
```

---

## Step 2 — Install Nextflow

Nextflow requires Java 11+.

### Option A — Official installer (recommended, no conda needed)

```bash
# Install Java if not present
sudo apt install -y default-jdk

# Install Nextflow via the official installer
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Verify
nextflow -version
```

### Option B — micromamba environment

```bash
micromamba create -n nextflow -c conda-forge -c bioconda nextflow>=24.10.5

# Verify
micromamba run -n nextflow nextflow -version
```

Set the conda environment cache to your **personal** directory:

```bash
# Add to ~/.bashrc or ~/.zshrc — use your own path, not a shared location
export NXF_CONDA_CACHEDIR="${HOME}/nf_conda_envs"
export NXF_CONDA_USEMAMBA=true
```

---

## Step 3 — Clone the pipeline

```bash
# Clone to the shared location
git clone https://github.com/nf-core/rnadnavar.git \
    /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar

# Or pull latest changes if already cloned
cd /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar
git pull
```

---

## Step 4 — Prepare reference databases

All reference files are already available at the shared `bio_db`. If setting up
from scratch, the following databases are required:

| Database | Shared path |
|----------|-------------|
| GRCh38 FASTA + FAI | `bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/` |
| BWA index | `bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/` |
| STAR index | `bio_db/star/` |
| HISAT2 index | `bio_db/hisat2/` |
| dbSNP | `bio_db/references/.../GATKBundle/dbsnp_146.hg38.vcf.gz` |
| Panel of Normals | `bio_db/references/.../GATKBundle/1000g_pon.hg38.vcf.gz` |
| gnomAD germline | `bio_db/references/.../GATKBundle/af-only-gnomad.hg38.vcf.gz` |
| VEP cache (v115) | `bio_db/vep/` |
| GTF (Gencode) | `bio_db/references/.../Genes.gencode/genes.gtf` |
| REDIportal (RNA editing) | `bio_db/rna_editing/REDIportal/` |
| COSMIC | `bio_db/COSMIC/` |
| gnomAD exomes | `bio_db/gnomAD/exomes/` |
| Intervals BED | `bio_db/intervals/ukb.pad50.broad.pad50.union.bed` |

---

## Step 5 — Install DeepSomatic locally (no Docker)

DeepSomatic r1.9 must be installed as local binaries because the pipeline runs
without Docker on this cluster.

### Requirements

- Ubuntu 22.04
- Python 3.10.12
- CPUs with SSE4 and AVX support
- `sudo` privilege

### Option A — Fast install via rsync (recommended)

Pre-configured binaries are available at the shared path:

```bash
sudo rsync -avP /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/DeepVariant-1.9.0/ /opt/deepvariant/bin/
```

Verify:

```bash
/opt/deepvariant/bin/run_deepsomatic --help 2>&1 | head -5
```

### Option B — Full manual installation

#### 5.1 Install system dependencies

```bash
sudo apt update
sudo apt install -y apt-utils build-essential python3-dev python3-pip python3-pip-whl \
    libcairo2-dev libgirepository1.0-dev pkg-config libdbus-1-dev parallel
```

#### 5.2 Download prebuilt DeepSomatic binaries

Binaries are available at the shared path:
```
/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/DeepVariant-1.9.0
```

Or download from Google Cloud:
```bash
gsutil -m rsync -r "gs://deepvariant/binaries/DeepVariant/1.9.0/DeepVariant-1.9.0" .
```

#### 5.3 Download DeepSomatic models

Models are available at the shared path:
```
/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/models
```

Or download from Google Cloud:
```bash
DEST="./models/deepsomatic/1.9.0/"
mkdir -p "${DEST}"
gsutil -m rsync -r "gs://deepvariant/models/DeepSomatic/1.9.0/" "${DEST}"
```

#### 5.4 Prepare system environment

```bash
# Run the prerequisite setup script
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/dv_tf/run-prereq.sh
```

#### 5.5 Prepare Python environment for DeepSomatic

**Option B1 — Use system Python:**

```bash
/usr/bin/python3 -m pip install -r requirements.txt --no-deps \
    -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple
```

**Option B2 — Use micromamba environment (recommended if using custom Python):**

```bash
micromamba create -n tf \
    -c conda-forge -c nvidia \
    tensorflow=2.13.1=cuda118py310h189a05f_1 \
    python=3.10.12 cudatoolkit cudnn uv \
    pygobject=3.42.1 pkg-config pkgconfig

PIP_USER=false micromamba run -n tf uv pip install -r requirements.txt --no-deps \
    -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple
```

> If you encounter tensorflow installation issues, check `~/.condarc` for
> any custom channels that may conflict:
> ```bash
> cat ~/.condarc
> ```

#### 5.6 Install DeepSomatic to /opt

```bash
mkdir -p /opt/deepvariant/bin/deepsomatic/

# Copy run_deepsomatic.py from shared location
cp /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/deepvariant/scripts/run_deepsomatic.py \
    /opt/deepvariant/bin/deepsomatic/

# Copy prebuilt binaries
rsync -avP /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/DeepVariant-1.9.0/ \
    /opt/deepvariant/bin/

# Patch binaries if using non-system Python
# Edit BIN_DIR and PYTHON variables in the script before running
bash patch/patch_dv_stub.sh

# Install models
mkdir -p /opt/models/deepsomatic
# Edit SRC_BASE to:
#   /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/models/deepsomatic/1.9.0/savedmodels
# Change rsync -avPn to rsync -avP in the script before running
bash rsync_install_deepsomatic_models.sh

# Generate CLI wrappers
# Edit PYTHON variable if using non-system Python
bash make_cli.sh
```

#### 5.7 Backup installed binaries

After a successful installation, back up to the shared location:

```bash
sudo rsync -avP /opt/deepvariant/ /t9k/mnt/joey/bio_utilities/deepvariant_opt/deepvariant/
sudo rsync -avP /opt/models/ /t9k/mnt/joey/bio_utilities/deepvariant_opt/models/
```

---

## Step 6 — Configure the pipeline

The shared config at
`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/input/test.rdv.shared.config`
is pre-configured for this cluster. For your own runs, copy and adapt
`examples/seq2neo/seq2neo.shared.config` or `examples/neoantigen/neoantigen.shared.config`.

Key environment variables to set in your shell or runner script:

```bash
# Use your own personal conda env cache — not a shared path
export NXF_CONDA_CACHEDIR="${HOME}/nf_conda_envs"
export NXF_CONDA_USEMAMBA=true
export HTTPS_PROXY="http://10.233.17.241:3128"
```

---

## Step 7 — Validate the installation

Run the shared test dataset (COO8801, small subset) to confirm everything works:

```bash
# Dry run first — preview the command
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/shared_test/run.sh \
    --dry-run

# Full run
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/shared_test/run.sh
```

Or manually (equivalent to what the script runs):

```bash
RDV_TEST_DIR="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test"

HTTPS_PROXY="http://10.233.17.241:3128" \
NXF_CONDA_CACHEDIR="${HOME}/nf_conda_envs" \
NXF_CONDA_USEMAMBA=true \
micromamba run -n nextflow nextflow run \
    /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/main.nf \
    -c "${RDV_TEST_DIR}/C008801/input/test.rdv.shared.config" \
    --input  "${RDV_TEST_DIR}/C008801/input/test.rdv.shared.csv" \
    --outdir "${RDV_TEST_DIR}/output/COO8801.shared" \
    -offline -with-conda -resume
```

> If Nextflow is on your `PATH` directly, omit `micromamba run -n nextflow`:
> ```bash
> nextflow run /t9k/mnt/.../main.nf -c ... --input ... --outdir ...
> ```

### Expected outputs

A successful run produces (based on the COO8801 test dataset with neoantigen workflow):

```
${RDV_TEST_DIR}/output/COO8801.shared/
├── preprocessing/
│   ├── mapped/
│   │   ├── COO8801DN/  COO8801DN.sorted.bam(.bai)
│   │   ├── COO8801DT/  COO8801DT.sorted.bam(.bai)
│   │   └── COO8801RT/  COO8801RT.bam(.bai)
│   ├── markduplicates/
│   │   ├── COO8801DN/  COO8801DN.md.cram(.crai)
│   │   ├── COO8801DT/  COO8801DT.md.cram(.crai)
│   │   └── COO8801RT/  COO8801RT.md.cram(.crai)
│   ├── recal_table/    COO8801DN/DT/RT.recal.table
│   ├── recalibrated/   COO8801DN/DT/RT.recal.cram(.crai)
│   └── splitncigarreads/
│       └── COO8801RT/  COO8801RT.sncr.cram(.crai)
├── variant_calling/
│   ├── mutect2/
│   │   ├── COO8801DT_vs_COO8801DN/  *.mutect2.filtered.vcf.gz, *.mutect2.vcf.gz
│   │   └── COO8801RT_vs_COO8801DN/  *.mutect2.filtered.vcf.gz, *.mutect2.vcf.gz
│   ├── strelka/
│   │   ├── COO8801DT_vs_COO8801DN/  *.strelka.variants.vcf.gz
│   │   └── COO8801RT_vs_COO8801DN/  *.strelka.variants.vcf.gz
│   └── deepsomatic/
│       ├── COO8801DT_vs_COO8801DN/  *.deepsomatic.vcf.gz, *.deepsomatic.g.vcf.gz
│       └── COO8801RT_vs_COO8801DN/  *.deepsomatic.vcf.gz, *.deepsomatic.g.vcf.gz
├── normalized/
│   ├── mutect2/    COO8801DT_vs_COO8801DN/ + COO8801RT_vs_COO8801DN/
│   │               *.mutect2.variants.dec.norm.vcf.gz(.csi/.tbi)
│   ├── strelka/    *.strelka.variants.dec.norm.vcf.gz(.csi/.tbi)
│   └── deepsomatic/ *.deepsomatic.variants.dec.norm.vcf.gz(.csi/.tbi)
├── consensus/
│   ├── COO8801DT_vs_COO8801DN/  *.consensus.vcf.gz(.tbi)
│   └── COO8801RT_vs_COO8801DN/  *.consensus.vcf.gz(.tbi)
├── filtered/
│   ├── COO8801DT_vs_COO8801DN/  *.filtered.vcf.gz, *.filtered.vcf.stripped.vcf.gz
│   └── COO8801RT_vs_COO8801DN/  *.filtered.vcf.gz, *.filtered.vcf.stripped.vcf.gz
├── rescue/
│   └── COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/
│       ├── *.rescued.vcf.gz(.tbi)
│       ├── *.rescue.rna_annotated.vcf.gz(.tbi)
│       ├── *.rescue.cosmic_gnomad_annotated.cosmic.vcf.gz(.tbi)
│       ├── *.rescue.cosmic_gnomad_annotated.gnomad.vcf.gz(.tbi)
│       ├── *.rescue.cosmic_gnomad_annotated.cosmic_gnomad_annotated.vcf.gz(.tbi)
│       ├── *.rescue.cosmic_gnomad_annotated.final.vcf.gz(.tbi)
│       ├── *.rescue.filtered.stripped.vep.vcf.gz(.tbi)
│       ├── *.filtered.vcf.gz(.tbi)
│       └── *.filtered.vcf.stripped.vcf.gz(.tbi)
├── neoantigen/                              (only when enable_neoantigen_workflow=true)
│   ├── COO8801DT_vs_COO8801DN/
│   │   └── *.mutect2.variants.dec.norm.vcf.gz(.tbi)   — neoantigen-ready VCF
│   └── COO8801RT-LX/
│       ├── quant.sf, quant.tsv              — Salmon quantification
│       ├── lib_format_counts.json, cmd_info.json
│       └── aux_info/  meta_info.json, ambig_info.tsv, fld.gz, bias files
├── vcf_realignment/
│   ├── preprocessing/
│   │   ├── mapped/COO8801RT/
│   │   ├── markduplicates/COO8801RT_realign/
│   │   └── splitncigarreads/COO8801RT_realign/
│   ├── readids/COO8801RT/  COO8801RT_IDs_all.txt
│   ├── variant_calling/    mutect2/strelka/deepsomatic for COO8801RT_realign_vs_COO8801DN
│   ├── normalized/         *.dec.norm.vcf.gz for realigned RNA
│   ├── consensus/COO8801RT_realign_vs_COO8801DN/  *.consensus.vcf.gz(.tbi)
│   ├── filtered/COO8801RT_realign_vs_COO8801DN/   *.filtered.vcf.gz, *.stripped.vcf.gz
│   ├── rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN/
│   │   └── (same structure as rescue/ above)
│   └── vcf2bed/COO8801RT/  COO8801RT.bed
└── pipeline_info/
    └── execution_trace_*.txt
```

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `NoChannelsConfiguredError` in conda | No channels in conda env | Use `environment.yml` with explicit channels (already fixed in modules) |
| `--validateMappings` specified more than once | Duplicate flag in `ext.args` and script | Remove hardcoded flag from module script |
| `Conda environment file does not exist` | Missing `environment.yml` | Create the file in the module directory |
| DeepSomatic `command not found` | Not installed or not in PATH | Run Option A rsync install or check `/opt/deepvariant/bin/` |
| `CUDA out of memory` | GPU memory insufficient | DeepSomatic runs on CPU by default; check TF config |
| Nextflow `offline` mode fails | Missing cached conda envs | Run once with internet access to populate `NXF_CONDA_CACHEDIR` |
| `micromamba: command not found` | micromamba not installed | Follow Step 1, or use `conda` as a drop-in replacement |
| Pipeline writes output into repo | Relative `OUTDIR` used from repo root | Always pass `--outdir` with an absolute path |
