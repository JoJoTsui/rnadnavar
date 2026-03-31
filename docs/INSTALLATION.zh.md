# rnadnavar — 安装指南

## 概述

本指南介绍如何在 Linux 服务器（Ubuntu 22.04）上无需 Docker、使用 conda/micromamba
管理环境的方式完整安装 rnadnavar 流水线。

需要安装三个组件：

1. **Nextflow** — 工作流引擎
2. **rnadnavar 流水线** — 流水线代码和 conda 环境
3. **DeepSomatic** — 本地二进制安装（无需 Docker）

集群上已有所有共享资源（参考数据库、conda 环境、预安装二进制文件）。
如果您在新服务器上安装，请按完整指南操作。
如果您加入现有集群，请直接跳至[加入现有集群](#加入现有集群)。

---

## 共享基础设施（现有集群）

| 资源 | 路径 |
|------|------|
| 流水线代码库 | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/` |
| 参考数据库 | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/` |
| 测试数据集（COO8801） | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/` |
| Conda 环境 | `/t9k/mnt/joey/nf_conda_envs/` |
| DeepSomatic 二进制文件 | `/opt/deepvariant/` |
| DeepSomatic 模型 | `/opt/models/deepsomatic/` |
| 预配置二进制文件（rsync 来源） | `/t9k/mnt/joey/bio_utilities/deepvariant_opt/` |

---

## 加入现有集群

如果上述共享基础设施已就绪，您只需：

1. 确认 micromamba 和 `nextflow` 环境可访问
2. 运行共享测试以确认访问权限

```bash
# 1. 检查 micromamba 是否可用
micromamba --version

# 2. 检查 nextflow 环境是否存在
micromamba env list | grep nextflow

# 3. 检查 nextflow 是否正常工作
micromamba run -n nextflow nextflow -version

# 4. 检查 DeepSomatic 是否已安装
/opt/deepvariant/bin/run_deepsomatic --help 2>&1 | head -5

# 5. 运行共享测试（先预演）
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/shared_test/run.sh --dry-run
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/shared_test/run.sh
```

---

## 完整安装（新服务器）

### 环境要求

- Ubuntu 22.04
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

# 或通过 conda-forge 安装
conda install -c conda-forge micromamba

# 验证
micromamba --version
```

---

## 步骤 2 — 安装 Nextflow

Nextflow 需要 Java 11+。

```bash
# 如未安装 Java，先安装
sudo apt install -y default-jdk

# 在专用 micromamba 环境中安装 Nextflow
micromamba create -n nextflow -c conda-forge -c bioconda nextflow>=24.10.5

# 验证
micromamba run -n nextflow nextflow -version
```

> **说明：** 如果 nextflow 已在 `PATH` 中（例如系统级安装），可跳过 micromamba
> 环境，直接调用 `nextflow`。`micromamba run -n nextflow nextflow run ...`
> 包装器仅在 nextflow 安装于 conda 环境中时才需要。

将 conda 环境缓存设置为共享位置，使所有团队成员共享预构建的环境：

```bash
# 添加到 ~/.bashrc 或 ~/.zshrc
export NXF_CONDA_CACHEDIR="/t9k/mnt/joey/nf_conda_envs"
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

| 数据库 | 路径 |
|--------|------|
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

预配置的二进制文件已在共享路径中，这是最快的安装方式：

```bash
sudo rsync -avP /t9k/mnt/joey/bio_utilities/deepvariant_opt/ /opt/
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

#### 5.4 为 DeepSomatic 准备 Python 环境

**方案 B1 — 使用系统 Python：**

```bash
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/dv_tf/run-prereq.sh
/usr/bin/python3 -m pip install -r requirements.txt --no-deps \
    -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple
```

**方案 B2 — 使用 micromamba 环境（推荐）：**

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

#### 5.5 将 DeepSomatic 安装到 /opt

```bash
mkdir -p /opt/deepvariant/bin/deepsomatic/

# 复制 run_deepsomatic.py
cp /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/deepvariant/scripts/run_deepsomatic.py \
    /opt/deepvariant/bin/deepsomatic/

# 复制预编译二进制文件
rsync -avP /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/DeepVariant-1.9.0/ \
    /opt/deepvariant/bin/

# 如果使用非系统 Python，对二进制文件打补丁（修改脚本中的 BIN_DIR 和 PYTHON 变量）
bash patch/patch_dv_stub.sh

# 安装模型
mkdir -p /opt/models/deepsomatic
# 将 SRC_BASE 修改为：/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/models/deepsomatic/1.9.0/savedmodels
# 将脚本中的 rsync -avPn 改为 rsync -avP 后再运行
bash rsync_install_deepsomatic_models.sh

# 生成 CLI 包装器（根据需要修改 PYTHON 变量）
bash make_cli.sh
```

#### 5.6 备份已安装的二进制文件

安装成功后，备份到共享位置：

```bash
sudo rsync -avP /opt/deepvariant/ /t9k/mnt/joey/bio_utilities/deepvariant_opt/deepvariant/
sudo rsync -avP /opt/models/ /t9k/mnt/joey/bio_utilities/deepvariant_opt/models/
```

---

## 步骤 6 — 配置流水线

`rnadnavar_test/C008801/input/test.rdv.shared.config` 中的共享配置已针对本集群预配置。
对于自己的运行，可复制并修改 `examples/seq2neo/seq2neo.shared.config` 或
`examples/neoantigen/neoantigen.shared.config`。

在 shell 或运行脚本中设置以下关键环境变量：

```bash
export NXF_CONDA_CACHEDIR="/t9k/mnt/joey/nf_conda_envs"
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

或手动运行：

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

### 预期输出

成功运行后将产生：

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

## 常见问题排查

| 现象 | 可能原因 | 解决方法 |
|------|---------|---------|
| `NoChannelsConfiguredError` | conda 环境中未配置 channels | 使用包含显式 channels 的 `environment.yml`（模块中已修复） |
| `--validateMappings` 重复指定 | `ext.args` 和脚本中均有该标志 | 从模块脚本中删除硬编码的标志 |
| `Conda environment file does not exist` | 缺少 `environment.yml` | 在模块目录中创建该文件 |
| DeepSomatic `command not found` | 未安装或不在 PATH 中 | 运行方案 A rsync 安装，或检查 `/opt/deepvariant/bin/` |
| `CUDA out of memory` | GPU 内存不足 | DeepSomatic 默认在 CPU 上运行；检查 TF 配置 |
| Nextflow `offline` 模式失败 | 缺少缓存的 conda 环境 | 先在有网络的情况下运行一次，填充 `NXF_CONDA_CACHEDIR` |
