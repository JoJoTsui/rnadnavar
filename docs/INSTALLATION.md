# rnadnavar — Installation Guide

## Overview

This guide covers a full installation of the rnadnavar pipeline on a Linux server
(Ubuntu 22.04) without Docker, using conda/micromamba for environment management.

Three components need to be installed:

1. **Nextflow** — the workflow engine
2. **rnadnavar pipeline** — the pipeline code and conda environments
3. **DeepSomatic** — local binary installation (no Docker required)

All shared resources (reference databases, conda envs, pre-installed binaries) are
already available on this cluster. If you are setting up a new server, follow the
full guide. If you are joining an existing cluster, jump to
[Joining an existing cluster](#joining-an-existing-cluster).

---

## Shared infrastructure (existing cluster)

| Resource | Path |
|----------|------|
| Pipeline repo | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/` |
| Reference databases | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/` |
| Test dataset (COO8801) | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/` |
| Conda environments | `/t9k/mnt/joey/nf_conda_envs/` |
| DeepSomatic binaries | `/opt/deepvariant/` |
| DeepSomatic models | `/opt/models/deepsomatic/` |
| Pre-configured binaries (rsync source) | `/t9k/mnt/joey/bio_utilities/deepvariant_opt/` |

---

## Joining an existing cluster

If the shared infrastructure above is already in place, you only need to:

1. Verify micromamba and the `nextflow` environment are accessible
2. Run the shared test to confirm your access

```bash
# 1. Check micromamba is available
micromamba --version

# 2. Check the nextflow environment exists
micromamba env list | grep nextflow

# 3. Check nextflow works
micromamba run -n nextflow nextflow -version

# 4. Check DeepSomatic is installed
/opt/deepvariant/bin/run_deepsomatic --help 2>&1 | head -5

# 5. Run the shared test (dry run first)
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/shared_test/run.sh --dry-run
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/shared_test/run.sh
```

---

## Full installation (new server)

### Requirements

- Ubuntu 22.04
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

# Or via conda-forge
conda install -c conda-forge micromamba

# Verify
micromamba --version
```

---

## Step 2 — Install Nextflow

Nextflow requires Java 11+.

```bash
# Install Java if not present
sudo apt install -y default-jdk

# Install Nextflow into a dedicated micromamba environment
micromamba create -n nextflow -c conda-forge -c bioconda nextflow>=24.10.5

# Verify
micromamba run -n nextflow nextflow -version
```

> **Note:** If nextflow is already on your `PATH` (e.g., installed system-wide),
> you can skip the micromamba environment and call `nextflow` directly.
> The `micromamba run -n nextflow nextflow run ...` wrapper is only needed when
> nextflow lives inside a conda environment.

Set the conda environment cache to the shared location so all teammates share
pre-built environments:

```bash
# Add to ~/.bashrc or ~/.zshrc
export NXF_CONDA_CACHEDIR="/t9k/mnt/joey/nf_conda_envs"
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

| Database | Path |
|----------|------|
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

Pre-configured binaries are already available at the shared path. This is the
fastest option:

```bash
sudo rsync -avP /t9k/mnt/joey/bio_utilities/deepvariant_opt/ /opt/
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

#### 5.4 Prepare Python environment for DeepSomatic

**Option B1 — Use system Python:**

```bash
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/dv_tf/run-prereq.sh
/usr/bin/python3 -m pip install -r requirements.txt --no-deps \
    -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple
```

**Option B2 — Use micromamba environment (recommended):**

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

#### 5.5 Install DeepSomatic to /opt

```bash
mkdir -p /opt/deepvariant/bin/deepsomatic/

# Copy run_deepsomatic.py
cp /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/deepvariant/scripts/run_deepsomatic.py \
    /opt/deepvariant/bin/deepsomatic/

# Copy prebuilt binaries
rsync -avP /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/DeepVariant-1.9.0/ \
    /opt/deepvariant/bin/

# Patch binaries if using non-system Python (edit BIN_DIR and PYTHON in the script)
bash patch/patch_dv_stub.sh

# Install models
mkdir -p /opt/models/deepsomatic
# Edit SRC_BASE to: /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/models/deepsomatic/1.9.0/savedmodels
# Change rsync -avPn to rsync -avP in the script before running
bash rsync_install_deepsomatic_models.sh

# Generate CLI wrappers (edit PYTHON variable if needed)
bash make_cli.sh
```

#### 5.6 Backup installed binaries

After a successful installation, back up to the shared location:

```bash
sudo rsync -avP /opt/deepvariant/ /t9k/mnt/joey/bio_utilities/deepvariant_opt/deepvariant/
sudo rsync -avP /opt/models/ /t9k/mnt/joey/bio_utilities/deepvariant_opt/models/
```

---

## Step 6 — Configure the pipeline

The shared config at `rnadnavar_test/C008801/input/test.rdv.shared.config` is
pre-configured for this cluster. For your own runs, copy and adapt
`examples/seq2neo/seq2neo.shared.config` or `examples/neoantigen/neoantigen.shared.config`.

Key environment variables to set in your shell or runner script:

```bash
export NXF_CONDA_CACHEDIR="/t9k/mnt/joey/nf_conda_envs"
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

Or manually:

```bash
HTTPS_PROXY="http://10.233.17.241:3128" \
NXF_CONDA_CACHEDIR="/t9k/mnt/joey/nf_conda_envs" \
NXF_CONDA_USEMAMBA=true \
micromamba run -n nextflow nextflow run \
    /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/main.nf \
    -c /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/input/test.rdv.shared.config \
    --input  /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/input/test.rdv.shared.csv \
    --outdir /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/output/COO8801.shared \
    -offline -with-conda -resume
```

### Expected outputs

A successful run produces:

```
output/COO8801.shared/
├── variant_calling/
│   ├── mutect2/COO8801DT_vs_COO8801DN/   *.mutect2.filtered.vcf.gz
│   ├── strelka/COO8801DT_vs_COO8801DN/   *.strelka.variants.vcf.gz
│   └── deepsomatic/COO8801DT_vs_COO8801DN/ *.deepsomatic.vcf.gz
├── consensus/COO8801DT_vs_COO8801DN/     *.consensus.vcf.gz
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
