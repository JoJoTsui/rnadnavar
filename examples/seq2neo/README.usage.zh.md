# seq2neo — 使用指南

## 共享基础设施

`seq2neo.shared.config` 和 `config/runner.yaml` 中的所有路径均指向集群上的共享位置，团队所有成员均可访问：

| 资源 | 共享路径 |
|------|---------|
| 流水线代码库 | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/` |
| 参考数据库 | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/` |
| Conda 环境 | `/t9k/mnt/joey/nf_conda_envs/` |
| Nextflow 自定义配置 | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/configs/` |

`seq2neo_root`（数据/输出/状态目录）**不是**共享的——每位用户在 `config/runner.yaml` 中设置自己的路径。

---

## 分离根目录布局

脚本通过 `$BASH_SOURCE` 自动定位自身路径，代码中没有任何硬编码的仓库路径。
首次使用前，只需在 `config/runner.yaml` 中设置以下三个值：

| 键 | 说明 |
|----|------|
| `main_nf` | `rnadnavar/main.nf` 的绝对路径 — 共享路径：`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/main.nf` |
| `rdv_conf` | `seq2neo.shared.config` 的绝对路径 — 共享路径：`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/seq2neo.shared.config` |
| `seq2neo_root` | 个人数据/输出/状态根目录（非共享，每位用户自行设置） |

仓库可以放在任意路径，脚本始终能自动定位：
```bash
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
```

### 仓库目录结构（版本控制，运行时只读）

```
$REPO/examples/seq2neo/
├── PRJNA298376.txt                  # 原始样本清单
├── PRJNA298330.txt
├── PRJNA298310.txt
├── PRJNA298330.disease.tsv          # 每个患者的疾病注释
├── seq2neo.shared.config            # nextflow 流程配置
├── config/
│   └── runner.yaml                  # 执行配置（seq2neo_root → DATA）
├── scripts/
│   ├── lib/
│   │   ├── __init__.py
│   │   └── common.py                # 共享逻辑（疾病规则、分类、分区）
│   ├── parse_projects_to_json.py    # 步骤1：解析清单 → 合并 JSON
│   └── run_batch_from_json.py       # 步骤2：批量 nextflow 运行器
├── parse.sh                         # 便捷脚本：执行解析步骤
├── dry_run.sh                       # 便捷脚本：预览任意过滤条件
├── run_set1.sh                      # 结直肠癌
├── run_set2.sh                      # 结肠癌
├── run_set3.sh                      # 壶腹部/胆管/胆管癌/食管/黑色素瘤
├── run_set4.sh                      # 胃/肺/胰腺/直肠
└── run_single.sh                    # 单个患者
```

### 工作区目录结构（运行时生成，不纳入版本控制）

```
$DATA/
├── data/processed/
│   ├── merged.json                  # 统一样本数据库
│   └── set{1-4}_samples.tsv         # 分区成员列表
├── runs/
│   ├── csv/                         # 每个样本的 nextflow 输入 CSV
│   └── run_state.json               # 断点续跑状态
└── output/
    └── PRJNA298376_4060/            # 每个样本的 nextflow 输出
        └── ...
```

---

## 样本状态

| 状态       | 含义                                                      | 是否可用？ |
|------------|-----------------------------------------------------------|-----------|
| standard   | 恰好有 1 个 DN + 1 个 DT + 1 个 RT 配对                  | 是        |
| extra      | 三种模态均存在，但至少一种有 >1 个配对                    | 是 — 使用第一个配对 |
| incomplete | 缺少 DN、DT 或 RT 中的任意一种                            | 否        |

---

## 分区集合

| 集合  | 疾病                                                                          |
|-------|-------------------------------------------------------------------------------|
| set1  | 结直肠癌（严格匹配 — 结肠癌和直肠癌不包含在此）                               |
| set2  | 结肠癌                                                                        |
| set3  | 壶腹部癌 / 胆管癌 / 胆管细胞癌 / 食管癌 / 黑色素瘤                           |
| set4  | 胃癌 / 肺癌 / 胰腺癌 / 直肠癌                                                |

set2–4 中每种疾病只出现在一个集合中（疾病互斥）。
集合大小均衡是次要目标。

---

## 配置参考

### config/runner.yaml — 所有选项

```yaml
# ── Nextflow 流程路径 ─────────────────────────────────────────────────────
main_nf:  /path/to/rnadnavar/main.nf
#   rnadnavar 流程入口文件的绝对路径。
#   如需切换 nextflow 版本或仓库位置，修改此路径即可。

rdv_conf: /path/to/rnadnavar/examples/seq2neo/seq2neo.shared.config
#   本项目共享 nextflow 配置文件的绝对路径。
#   包含所有流程参数（基因组、工具、参考文件等）。
#   可在运行时通过 --rdv-conf 覆盖。

# ── Conda / Nextflow 环境 ─────────────────────────────────────────────────
nxf_conda_cachedir: /path/to/nf_conda_envs
#   nextflow 存储/复用 conda 环境的目录（对应 NXF_CONDA_CACHEDIR）。
#   设置为持久化的共享目录，避免每次运行重新创建环境。
#   这是 nextflow 的 conda 环境工作/缓存目录。

nxf_conda_usemamba: "true"
#   使用 mamba 代替 conda 加速环境创建（对应 NXF_CONDA_USEMAMBA）。

micromamba_env: nextflow
#   安装了 nextflow 的 micromamba 环境名称。
#   实际调用方式：micromamba run -n <micromamba_env> nextflow run ...

# ── 网络代理 ──────────────────────────────────────────────────────────────
https_proxy: "http://10.233.17.241:3128"
#   nextflow 和 conda 下载使用的 HTTP/HTTPS 代理。
#   留空（""）表示禁用。

# ── 数据根目录 ────────────────────────────────────────────────────────────
seq2neo_root: /path/to/seq2neo_data
#   所有运行时数据的根目录。以下所有相对路径均相对于此目录解析。

merged_json: data/processed/merged.json
#   输入：由 parse.sh（步骤1）生成的样本数据库。

csv_dir: runs/csv
#   每次运行前写入的每个样本 nextflow 输入 CSV 的目录。

outdir_base: output
#   每个样本 nextflow 输出的基础目录。
#   每个样本写入：<outdir_base>/<project>_<patient>/

state_file: runs/run_state.json
#   记录每个样本运行状态（succeeded/failed/running）的 JSON 文件。
#   用于断点续跑逻辑 — 除非要重置所有状态，否则不要删除。

# ── 执行选项 ──────────────────────────────────────────────────────────────
lane: LX
#   写入 nextflow 输入 CSV 的泳道标识符。

dry_run: false
#   若为 true，打印命令但不执行。等同于 --dry-run 标志。

resume: true
#   向 nextflow 传递 -resume，启用 nextflow 内置的任务缓存。

offline: true
#   向 nextflow 传递 -offline，阻止远程资源获取。

max_parallel: 1
#   并发运行的样本数量。1 = 顺序执行（安全默认值）。

# ── 重试策略 ──────────────────────────────────────────────────────────────
max_retries: 1
#   失败样本在被永久跳过前的最大重试次数。

# ── 完成检查 ──────────────────────────────────────────────────────────────
completion_artifacts:
  - "**/*.filtered.vcf.gz"
  - "**/pipeline_info/execution_trace*.txt"
#   所有 glob 必须在 outdir/<project>_<patient>/ 下匹配到至少1个文件，
#   样本才被认为成功完成。如果流程输出文件不同，请相应调整。
```

### seq2neo.shared.config — 流程参数

此文件是通过 `-c` 传递给 nextflow 的配置文件，设置所有流程级参数。主要部分：

```
资源限制        — 每个进程的 cpus、memory、time
参考基因组      — fasta、fasta_fai、genome（GRCh38）
比对工具        — bwa index、star_index、hisat2_index、aligner
注释数据库      — dbsnp、pon、germline_resource、vep_cache、gtf
RNA 编辑        — rediportal_vcf、min_rna_support
COSMIC / gnomAD — cosmic_database、gnomad_database、阈值参数
流程工具        — tools = deepsomatic,mutect2,strelka,vep,norm,...
起始步骤        — step = mapping（从原始 FASTQ 开始）
```

如需修改参考文件路径或启用/禁用工具，直接编辑此文件。
如需使用完全不同的配置文件，更新 `runner.yaml` 中的 `rdv_conf`。

---

## 如何修改 nextflow 路径

nextflow 可执行文件通过 micromamba 调用：

```
micromamba run -n <micromamba_env> nextflow run <main_nf> ...
```

有两个独立的路径需要配置：

**1. 流程脚本（`main_nf`）** — 运行哪个 `main.nf`：
```yaml
# config/runner.yaml
main_nf: /new/path/to/rnadnavar/main.nf
```
或在运行时通过命令行覆盖，无需修改配置文件：
```bash
bash run_set2.sh --main-nf /new/path/to/rnadnavar/main.nf
```

**2. nextflow 可执行文件** — 使用哪个 nextflow 安装：
```yaml
# config/runner.yaml
micromamba_env: nextflow   # 包含 nextflow 的 micromamba 环境
```
如果 nextflow 在不同的 conda/micromamba 环境中，修改 `micromamba_env`。
如果使用其他启动方式（例如 PATH 中的 `nextflow`），需要编辑
`scripts/run_batch_from_json.py` 中的 `build_command()` 函数。

---

## 如何修改 nextflow 工作目录（缓存目录）

Nextflow 使用两种缓存/工作目录：

**1. Conda 环境缓存（`nxf_conda_cachedir`）** — conda 环境的存储位置：
```yaml
# config/runner.yaml
nxf_conda_cachedir: /shared/path/nf_conda_envs
```
对应 `NXF_CONDA_CACHEDIR` 环境变量。设置为持久化的共享目录，
避免每次运行重新创建环境。

**2. Nextflow 工作目录** — 任务中间文件的缓存位置。
Nextflow 默认使用命令运行位置的 `./work` 目录。
如需修改，可在 `seq2neo.shared.config` 中添加 `-work-dir` 参数，
或通过环境变量 `NXF_WORK` 设置：
```bash
export NXF_WORK=/scratch/nextflow_work
bash run_set2.sh
```
也可以在运行前临时设置：
```bash
NXF_WORK=/scratch/nextflow_work bash run_set2.sh
```

---

## 步骤1 — 解析清单

先编辑 `config/runner.yaml`（见上文），然后：

```bash
# 从任意位置运行 — 脚本自动定位
bash /your/path/rnadnavar/examples/seq2neo/parse.sh

# 或者已在 examples/seq2neo 目录中：
bash parse.sh
```

从 `seq2neo_root`（在 `runner.yaml` 中配置）读取清单文件：
- `PRJNA298376.txt` — 树形格式清单
- `PRJNA298330.txt` — 管道格式清单
- `PRJNA298310.txt` — 管道格式清单（全部为黑色素瘤）
- `PRJNA298330.disease.tsv` — 每个患者的疾病注释

写入 `$seq2neo_root/data/processed/`：
- `merged.json` — 包含模态路径、状态、分区集合的统一样本数据库
- `set1_samples.tsv` 至 `set4_samples.tsv` — 分区成员列表

打印验证报告：状态统计、集合统计、疾病互斥性检查、不完整样本列表。

---

## 步骤2 — 配置

编辑 `config/runner.yaml` — 首次使用前只需设置三个字段：

```yaml
main_nf:      /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/main.nf
rdv_conf:     /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/seq2neo.shared.config
seq2neo_root: /your/data/output/root   # 个人路径，每位用户自行设置
```

其他所有选项均有合理的默认值。完整配置参考见上文。

---

## 步骤3 — 预运行（预览，不执行）

```bash
SEQ2NEO_SCRIPTS=/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo

# 所有可用样本
bash $SEQ2NEO_SCRIPTS/dry_run.sh

# 仅 set1（结直肠癌）
bash $SEQ2NEO_SCRIPTS/dry_run.sh --set 1

# set1，仅 standard 样本
bash $SEQ2NEO_SCRIPTS/dry_run.sh --set 1 --status standard

# 单个患者
bash $SEQ2NEO_SCRIPTS/dry_run.sh --project PRJNA298376 --patient 4060

# 按疾病子字符串过滤
bash $SEQ2NEO_SCRIPTS/dry_run.sh --disease "Pancreatic"
```

预运行会打印每个样本将要执行的完整 nextflow 命令，包括 CSV 路径、
输出目录和所有环境变量。不会写入任何文件，也不会启动 nextflow 进程。

---

## 步骤4 — 运行

```bash
SEQ2NEO_SCRIPTS=/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo

# 按集合运行
bash $SEQ2NEO_SCRIPTS/run_set1.sh   # 结直肠癌
bash $SEQ2NEO_SCRIPTS/run_set2.sh   # 结肠癌
bash $SEQ2NEO_SCRIPTS/run_set3.sh   # 壶腹部/胆管/胆管癌/食管/黑色素瘤
bash $SEQ2NEO_SCRIPTS/run_set4.sh   # 胃/肺/胰腺/直肠

# 运行单个患者
PROJECT=PRJNA298376 PATIENT=4060 bash $SEQ2NEO_SCRIPTS/run_single.sh

# 不修改 runner.yaml，通过命令行覆盖配置
bash $SEQ2NEO_SCRIPTS/run_set2.sh \
    --main-nf /other/path/main.nf \
    --outdir-base /scratch/output
```

### 命令行参数（所有运行脚本通用）

所有运行脚本都会将额外参数传递给 `run_batch_from_json.py`。
可用的覆盖选项：

| 参数 | 覆盖配置项 | 说明 |
|------|-----------|------|
| `--main-nf PATH` | `main_nf` | `main.nf` 的路径 |
| `--rdv-conf PATH` | `rdv_conf` | nextflow 配置文件路径 |
| `--outdir-base PATH` | `outdir_base` | 输出基础目录 |
| `--seq2neo PATH` | `seq2neo_root` | 数据根目录 |
| `--dry-run` | `dry_run` | 预览模式，不执行 |
| `--no-resume` | `resume` | 禁用 `-resume` 标志 |
| `--project ID` | — | 按项目 ID 过滤 |
| `--patient ID` | — | 按患者 ID 过滤 |
| `--disease STR` | — | 按疾病过滤（子字符串匹配） |
| `--set N` | — | 按分区集合过滤（1–4） |
| `--status STR` | — | 按状态过滤：`standard`、`extra` 或 `all` |

---

## 输出说明

### 每个样本的 nextflow 输出

每个样本写入 `$seq2neo_root/output/<project>_<patient>/`，例如：
```
output/PRJNA298376_4060/
├── pipeline_info/
│   ├── execution_trace_*.txt        # nextflow 任务级别的时间和资源使用情况
│   ├── execution_report_*.html      # 可视化执行报告
│   └── execution_timeline_*.html    # 时间线视图
├── <sample>/
│   ├── *.filtered.vcf.gz            # 最终过滤后的变异调用结果（完成标志文件）
│   ├── *.filtered.vcf.gz.tbi
│   └── ...                          # BAM 文件、中间 VCF、QC 文件
```

当 `completion_artifacts` 中的所有 glob 在输出目录下至少匹配到一个文件时，
样本被认为已成功完成（可在 `runner.yaml` 中配置）。

### 运行状态文件

`$seq2neo_root/runs/run_state.json` — 记录每个样本的状态：
```json
{
  "PRJNA298376::4060": {
    "status": "succeeded",
    "finished": "2024-01-15T10:23:45"
  }
}
```
可能的状态值：`running`（运行中）、`succeeded`（成功）、`failed`（失败）。

### 每个样本的输入 CSV

`$seq2neo_root/runs/csv/<project>_<patient>.csv` — 每次运行前自动生成：
```
patient,status,sample,lane,fastq_1,fastq_2
PRJNA298376_4060,0,PRJNA298376_4060DN,LX,/path/DN_1.fastq.gz,/path/DN_2.fastq.gz
PRJNA298376_4060,1,PRJNA298376_4060DT,LX,/path/DT_1.fastq.gz,/path/DT_2.fastq.gz
PRJNA298376_4060,2,PRJNA298376_4060RT,LX,/path/RT_1.fastq.gz,/path/RT_2.fastq.gz
```
状态码：DN=0（正常 DNA），DT=1（肿瘤 DNA），RT=2（肿瘤 RNA）。

### 解析输出

`$seq2neo_root/data/processed/`：
- `merged.json` — 包含模态路径、状态、分区分配的完整样本数据库
- `set1_samples.tsv` 至 `set4_samples.tsv` — 每个集合的 TSV 文件，包含项目、患者、疾病、状态

---

## 断点续跑

状态记录在 `$seq2neo_root/runs/run_state.json` 中。重复运行同一命令是安全的 — 已完成的样本会自动跳过。

**完成检测方式（严格）：**
- `$seq2neo_root/output/<project>_<patient>/` 目录必须存在
- `runner.yaml` 中 `completion_artifacts` 的所有 glob 必须匹配到至少1个文件

**中断后续跑：**
```bash
# 直接重新运行 — 已完成的样本自动跳过
bash $SEQ2NEO_SCRIPTS/run_set2.sh
```

**强制重新运行单个样本：**
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

**重置所有失败样本以重试：**
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

## 手动 nextflow 调用（单个样本）

每个样本的 CSV 写入 `$seq2neo_root/runs/csv/<project>_<patient>.csv`。
手动运行方式：

> **关于 micromamba 的说明：** `micromamba run -n nextflow` 包装器仅在 nextflow 安装于 micromamba/conda 环境中时才需要。如果 nextflow 已在 `PATH` 中，可直接运行：
> ```bash
> nextflow run /path/to/main.nf -c /path/to/config ...
> ```

```bash
REPO=/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar
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

如需为此次调用指定 nextflow 工作目录：
```bash
NXF_WORK=/scratch/nextflow_work \
micromamba run -n nextflow nextflow run ...
```
