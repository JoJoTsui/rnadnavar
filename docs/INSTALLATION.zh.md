# rnadnavar — 安装指南

## 概述

本指南介绍如何在 Linux 服务器（Ubuntu 22.04）上无需 Docker、使用 conda/micromamba
管理环境的方式完整安装 rnadnavar 流水线。

需要安装三个组件：

1. **Nextflow** — 工作流引擎
2. **rnadnavar 流水线** — 流水线代码和 conda 环境
3. **DeepSomatic** — 本地二进制安装（无需 Docker）

共享存储（参考数据库、流水线代码、测试数据、预构建 DeepSomatic 二进制文件）对所有
团队成员开放。但是，**每位用户必须自行安装运行环境** — Nextflow、micromamba 和
conda 环境是每用户独立的，不共享。

---

## 共享基础设施（共享存储——所有人可访问）

| 资源 | 路径 |
|------|------|
| 流水线代码库 | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/` |
| 参考数据库 | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/` |
| 测试数据集（COO8801） | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/` |
| DeepSomatic 二进制文件（rsync 来源） | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/DeepVariant-1.9.0/` |
| DeepSomatic 模型 | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/models/` |

> **注意：** Conda 环境和 Nextflow 运行时**不共享**。每位用户必须自行安装。
> 详见下方 [NXF_CONDA_CACHEDIR](#nxf_conda_cachedir) 说明。

---

## 注意事项

### Nextflow 调用方式

示例脚本使用 `micromamba run -n nextflow nextflow run ...` 调用 Nextflow。
此包装器仅在 Nextflow 安装于 micromamba/conda 环境中时才需要。
如果 Nextflow 已在 `PATH` 中（例如系统级安装或通过官方安装器安装），可直接调用：

```bash
nextflow run /path/to/main.nf -c /path/to/config ...
```

### micromamba 与 conda

`micromamba` 是 `conda` 的快速独立替代品。本指南所有命令均使用 `micromamba`。
如果您只有 `conda`，将 `micromamba` 替换为 `conda` 即可——语法完全相同。

### 输出目录

示例脚本中的默认 `OUTDIR` 为共享 `rnadnavar_test` 目录下的绝对路径。
对于自己的运行，请始终使用代码库外的绝对路径，避免将输出写入代码库。

### 离线模式

流水线使用 `-offline` 运行，防止 Nextflow 获取远程资源。
在新机器上首次运行需要网络访问（或已配置代理）以填充 `NXF_CONDA_CACHEDIR` 中的 conda 环境。
后续运行可离线进行。

### NXF_CONDA_CACHEDIR

将此变量设置为**您个人用户空间**中的目录，Nextflow 将在此构建和缓存 conda 环境。
请勿指向共享位置——conda 环境是每用户独立的，必须在本地构建。

```bash
# 示例：使用 home 目录或个人工作空间下的目录
export NXF_CONDA_CACHEDIR="${HOME}/nf_conda_envs"
# 或专用工作空间路径：
export NXF_CONDA_CACHEDIR="/t9k/mnt/hdd/work/${USER}/nf_conda_envs"
```

---

## 安装步骤

每位用户必须完成以下所有步骤，没有可加入的共享环境。

### 环境要求
- Python 3.10.12
- CPU 支持 SSE4 和 AVX 指令集
- `sudo` 权限
- 网络访问（或已配置代理）
- 建议最低配置：36 CPU，160 GB 内存

---

## 步骤 1 — 安装 micromamba

```bash
# 安装 micromamba（conda 的快速替代品）
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# 或在已有 conda 的情况下通过 conda-forge 安装
conda install -c conda-forge micromamba

# 验证
micromamba --version
```

---

## 步骤 2 — 安装 Nextflow

Nextflow 需要 Java 11+。

### 方案 A — 官方安装器（推荐，无需 conda）

```bash
# 如未安装 Java，先安装
sudo apt install -y default-jdk

# 通过官方安装器安装 Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# 验证
nextflow -version
```

### 方案 B — micromamba 环境

```bash
micromamba create -n nextflow -c conda-forge -c bioconda nextflow>=24.10.5

# 验证
micromamba run -n nextflow nextflow -version
```

将 conda 环境缓存设置为**个人**目录：

```bash
# 添加到 ~/.bashrc 或 ~/.zshrc — 使用您自己的路径，不要使用共享位置
export NXF_CONDA_CACHEDIR="${HOME}/nf_conda_envs"
export NXF_CONDA_USEMAMBA=true
```

---

## 步骤 3 — 克隆流水线代码库

```bash
# 克隆到共享位置
git clone https://github.com/nf-core/rnadnavar.git \
    /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar

# 或如果已克隆，拉取最新更改
cd /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar
git pull
```

---

## 步骤 4 — 准备参考数据库

所有参考文件已在共享 `bio_db` 中。如需从头搭建，需要以下数据库：

| 数据库 | 共享路径 |
|--------|---------|
| GRCh38 FASTA + FAI | `bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/` |
| BWA 索引 | `bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/` |
| STAR 索引 | `bio_db/star/` |
| HISAT2 索引 | `bio_db/hisat2/` |
| dbSNP | `bio_db/references/.../GATKBundle/dbsnp_146.hg38.vcf.gz` |
| Panel of Normals | `bio_db/references/.../GATKBundle/1000g_pon.hg38.vcf.gz` |
| gnomAD 胚系变异 | `bio_db/references/.../GATKBundle/af-only-gnomad.hg38.vcf.gz` |
| VEP 缓存（v115） | `bio_db/vep/` |
| GTF（Gencode） | `bio_db/references/.../Genes.gencode/genes.gtf` |
| REDIportal（RNA 编辑） | `bio_db/rna_editing/REDIportal/` |
| COSMIC | `bio_db/COSMIC/` |
| gnomAD 外显子组 | `bio_db/gnomAD/exomes/` |
| Intervals BED | `bio_db/intervals/ukb.pad50.broad.pad50.union.bed` |

---

## 步骤 5 — 本地安装 DeepSomatic（无需 Docker）

DeepSomatic r1.9 必须以本地二进制文件方式安装，因为本集群不使用 Docker。

### 环境要求

- Ubuntu 22.04
- Python 3.10.12
- CPU 支持 SSE4 和 AVX 指令集
- `sudo` 权限

### 方案 A — 通过 rsync 快速安装（推荐）

预编译的二进制文件已在共享路径中：

```bash
sudo rsync -avP /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/DeepVariant-1.9.0/ /opt/deepvariant/bin/
```

验证：

```bash
/opt/deepvariant/bin/run_deepsomatic --help 2>&1 | head -5
```

### 方案 B — 完整手动安装

#### 5.1 安装系统依赖

```bash
sudo apt update
sudo apt install -y apt-utils build-essential python3-dev python3-pip python3-pip-whl \
    libcairo2-dev libgirepository1.0-dev pkg-config libdbus-1-dev parallel
```

#### 5.2 下载预编译的 DeepSomatic 二进制文件

二进制文件已在共享路径中：
```
/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/DeepVariant-1.9.0
```

或从 Google Cloud 下载：
```bash
gsutil -m rsync -r "gs://deepvariant/binaries/DeepVariant/1.9.0/DeepVariant-1.9.0" .
```

#### 5.3 下载 DeepSomatic 模型

模型已在共享路径中：
```
/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/models
```

或从 Google Cloud 下载：
```bash
DEST="./models/deepsomatic/1.9.0/"
mkdir -p "${DEST}"
gsutil -m rsync -r "gs://deepvariant/models/DeepSomatic/1.9.0/" "${DEST}"
```

#### 5.4 准备系统环境

```bash
# 运行前置条件安装脚本
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/dv_tf/run-prereq.sh
```

#### 5.5 为 DeepSomatic 准备 Python 环境

**方案 B1 — 使用系统 Python：**

```bash
/usr/bin/python3 -m pip install -r requirements.txt --no-deps \
    -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple
```

**方案 B2 — 使用 micromamba 环境（推荐，适用于自定义 Python）：**

```bash
micromamba create -n tf \
    -c conda-forge -c nvidia \
    tensorflow=2.13.1=cuda118py310h189a05f_1 \
    python=3.10.12 cudatoolkit cudnn uv \
    pygobject=3.42.1 pkg-config pkgconfig

PIP_USER=false micromamba run -n tf uv pip install -r requirements.txt --no-deps \
    -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple
```

> 如果安装 tensorflow 时遇到问题，请检查 `~/.condarc` 中是否有可能冲突的自定义 channels：
> ```bash
> cat ~/.condarc
> ```

#### 5.6 将 DeepSomatic 安装到 /opt

```bash
mkdir -p /opt/deepvariant/bin/deepsomatic/

# 从共享位置复制 run_deepsomatic.py
cp /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/deepvariant/scripts/run_deepsomatic.py \
    /opt/deepvariant/bin/deepsomatic/

# 复制预编译二进制文件
rsync -avP /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/DeepVariant-1.9.0/ \
    /opt/deepvariant/bin/

# 如果使用非系统 Python，对二进制文件打补丁
# 运行前修改脚本中的 BIN_DIR 和 PYTHON 变量
bash patch/patch_dv_stub.sh

# 安装模型
mkdir -p /opt/models/deepsomatic
# 将 SRC_BASE 修改为：
#   /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/models/deepsomatic/1.9.0/savedmodels
# 将脚本中的 rsync -avPn 改为 rsync -avP 后再运行
bash rsync_install_deepsomatic_models.sh

# 生成 CLI 包装器
# 如使用非系统 Python，修改 PYTHON 变量
bash make_cli.sh
```

#### 5.7 备份已安装的二进制文件

安装成功后，备份到共享位置：

```bash
sudo rsync -avP /opt/deepvariant/ /t9k/mnt/joey/bio_utilities/deepvariant_opt/deepvariant/
sudo rsync -avP /opt/models/ /t9k/mnt/joey/bio_utilities/deepvariant_opt/models/
```

---

## 步骤 6 — 配置流水线

共享配置文件
`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/input/test.rdv.shared.config`
已针对本集群预配置。对于自己的运行，可复制并修改
`examples/seq2neo/seq2neo.shared.config` 或
`examples/neoantigen/neoantigen.shared.config`。

在 shell 或运行脚本中设置以下关键环境变量：

```bash
# 使用您个人的 conda 环境缓存目录，不要使用共享路径
export NXF_CONDA_CACHEDIR="${HOME}/nf_conda_envs"
export NXF_CONDA_USEMAMBA=true
export HTTPS_PROXY="http://10.233.17.241:3128"
```

---

## 步骤 7 — 验证安装

运行共享测试数据集（COO8801，小型子集）以确认一切正常：

```bash
# 先预演 — 预览命令
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/shared_test/run.sh \
    --dry-run

# 完整运行
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/shared_test/run.sh
```

或手动运行（与脚本等效）：

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

> 如果 Nextflow 已在 `PATH` 中，省略 `micromamba run -n nextflow`：
> ```bash
> nextflow run /t9k/mnt/.../main.nf -c ... --input ... --outdir ...
> ```

### 预期输出

成功运行后将产生（基于 COO8801 测试数据集，含 neoantigen 工作流）：

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
├── neoantigen/                              （仅在 enable_neoantigen_workflow=true 时）
│   ├── COO8801DT_vs_COO8801DN/
│   │   └── *.mutect2.variants.dec.norm.vcf.gz(.tbi)   — neoantigen 就绪 VCF
│   └── COO8801RT-LX/
│       ├── quant.sf, quant.tsv              — Salmon 定量结果
│       ├── lib_format_counts.json, cmd_info.json
│       └── aux_info/  meta_info.json, ambig_info.tsv, fld.gz, 偏差模型文件
├── vcf_realignment/
│   ├── preprocessing/
│   │   ├── mapped/COO8801RT/
│   │   ├── markduplicates/COO8801RT_realign/
│   │   └── splitncigarreads/COO8801RT_realign/
│   ├── readids/COO8801RT/  COO8801RT_IDs_all.txt
│   ├── variant_calling/    mutect2/strelka/deepsomatic（COO8801RT_realign_vs_COO8801DN）
│   ├── normalized/         *.dec.norm.vcf.gz（重比对 RNA）
│   ├── consensus/COO8801RT_realign_vs_COO8801DN/  *.consensus.vcf.gz(.tbi)
│   ├── filtered/COO8801RT_realign_vs_COO8801DN/   *.filtered.vcf.gz, *.stripped.vcf.gz
│   ├── rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN/
│   │   └── （结构同上方 rescue/ 目录）
│   └── vcf2bed/COO8801RT/  COO8801RT.bed
└── pipeline_info/
    └── execution_trace_*.txt
```

---

## 常见问题排查

| 现象 | 可能原因 | 解决方法 |
|------|---------|---------|
| `NoChannelsConfiguredError` | conda 环境中未配置 channels | 使用包含显式 channels 的 `environment.yml`（模块中已修复） |
| `--validateMappings` 重复指定 | `ext.args` 和脚本中均有该标志 | 从模块脚本中删除硬编码的标志 |
| `Conda environment file does not exist` | 缺少 `environment.yml` | 在模块目录中创建该文件 |
| DeepSomatic `command not found` | 未安装或不在 PATH 中 | 运行方案 A rsync 安装，或检查 `/opt/deepvariant/bin/` |
| `CUDA out of memory` | GPU 内存不足 | DeepSomatic 默认在 CPU 上运行；检查 TF 配置 |
| Nextflow `offline` 模式失败 | 缺少缓存的 conda 环境 | 先在有网络的情况下运行一次，填充 `NXF_CONDA_CACHEDIR` |
| `micromamba: command not found` | micromamba 未安装 | 按步骤 1 安装，或用 `conda` 替代 |
| 流水线输出写入代码库 | 使用了相对 `OUTDIR` | 始终通过 `--outdir` 传入绝对路径 |
