# Neoantigen 工作流 — 使用指南

## 概述

Neoantigen 准备工作流是 rnadnavar 流水线的一个可选分支，功能包括：

1. **FORMAT 字段标准化**：将 Mutect2、Strelka2、SAGE、DeepSomatic 等 DNA 变异检测工具的 FORMAT 字段统一为 Mutect2 兼容格式，并在共识 VCF 中添加 `AD_BY_CALLER` 和 `AF_BY_CALLER` INFO 字段。
2. **Salmon RNA 定量**：使用 `salmon quant` v1.11.4 对 RNA 肿瘤样本（status=2）的 FASTQ 文件进行准比对模式定量，与 STAR 比对并行运行。
3. **发布 neoantigen 就绪 VCF**：为每个 DNA 肿瘤样本（status=1）发布 VCF 文件，可直接用于 pVACseq、seq2neo 等工具。

所有新代码路径均通过 `enable_neoantigen_workflow = true` 控制。现有 `examples/seq2neo/` 配置不受影响。

---

## 共享基础设施

`neoantigen.shared.config` 和 `run.sh` 中的所有路径均指向集群上的共享位置，团队所有成员均可访问：

| 资源 | 共享路径 |
|------|---------|
| 流水线代码库 | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/` |
| 参考数据库 | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/` |
| 测试数据集（COO8801） | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/` |
| Conda 环境 | `/t9k/mnt/joey/nf_conda_envs/` |
| Salmon 索引（Gencode v49） | `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/salmon/salmon_index_gencode_v49_salmon_v1.11` |

运行测试数据集无需修改任何路径——所有内容均从上述共享位置解析。

---

## 输出目录注意事项

`run.sh` 中默认的 `OUTDIR` 为 `output/COO8801.neoantigen`，这是一个**相对路径**，相对于运行脚本时的当前目录。如果在代码库目录内运行，输出将写入代码库中，不建议这样做。

**推荐做法**：始终通过 `--outdir` 指定代码库外的绝对路径：

```bash
bash examples/neoantigen/run.sh \
    --outdir /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/neoantigen_output/COO8801
```

或在运行前设置输出根目录变量：

```bash
OUTROOT="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/neoantigen_output"
bash examples/neoantigen/run.sh --outdir "${OUTROOT}/COO8801"
```

相对路径默认值适用于快速本地测试。对于需要保留的正式运行，请使用绝对路径。

---

## 快速开始

```bash
# 1. 预演运行 — 预览 nextflow 命令，不实际执行
bash examples/neoantigen/run.sh --dry-run

# 2. 使用共享测试数据集运行（COO8801，小型子集）
#    推荐使用 --outdir 将输出写到代码库外
bash examples/neoantigen/run.sh \
    --outdir /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/neoantigen_output/COO8801

# 3. 使用自己的样本运行
bash examples/neoantigen/run.sh \
    --input /path/to/my_sample.csv \
    --outdir /path/to/output
```

`run.sh` 通过 `$BASH_SOURCE` 从脚本自身位置解析 `MAIN_NF` 和 `RDV_CONF`，因此只要代码库位于共享路径，从任意工作目录运行均可正常工作。

---

## 目录中的文件

```
examples/neoantigen/
├── neoantigen.shared.config   # 完整流水线配置（参考基因组 + 工具 + neoantigen 参数）
├── run.sh                     # 单样本启动脚本（调试 / 教程）
├── README.usage.md            # 英文使用指南
└── README.usage.zh.md         # 本文件（中文使用指南）
```

`neoantigen.shared.config` 包含所有标准流水线参数（基因组、工具、参考文件）以及 neoantigen 专用参数。所有参考路径均指向共享 `bio_db`。通常只需修改 `salmon_index`（如需使用不同于共享索引的版本）。

---

## 每次运行可能需要修改的参数

大多数参数已针对共享集群预配置。以下是可能需要调整的参数：

| 参数 | 位置 | 何时修改 |
|------|------|---------|
| `salmon_index` | `neoantigen.shared.config` | 仅在需要不同 Gencode 版本或自定义转录组时修改。`bio_db/salmon/` 中的共享索引可直接使用。 |
| `neoantigen_input_source` | `neoantigen.shared.config` | 仅在下游工具需要多检测工具联合 VCF 时改为 `'consensus'`（见下方 [VCF 来源选择](#vcf-来源选择mutect2-vs-consensus)）。 |
| `tools` | `neoantigen.shared.config` | 移除不需要的工具（例如，调试时移除 `realignment` 可加快运行速度）。 |
| `resourceLimits` | `neoantigen.shared.config` | 根据计算节点调整 `cpus`/`memory`/`time`。 |
| `--input` | `run.sh` CLI 参数 | 使用自己的样本 CSV。默认为共享 COO8801 测试数据集。 |
| `--outdir` | `run.sh` CLI 参数 | **对于需要保留的运行，始终设置为代码库外的绝对路径。** 默认 `output/COO8801.neoantigen` 是相对路径，从代码库根目录运行时会写入代码库内。 |

---

## VCF 来源选择：`mutect2` vs `consensus`

### 为什么推荐 `mutect2` 作为默认值

对于 pVACseq、seq2neo 等 neoantigen 预测工具，Mutect2 过滤后的 VCF 是正确的输入：

| 属性 | Mutect2 VCF | 共识 VCF |
|------|-------------|---------|
| 每样本 FORMAT 字段 | ✓ GT/AD/AF/DP/GQ/F1R2/F2R1/SB | ✗ FORMAT 列为 `.`（空） |
| 变异数量 | ~91（过滤后的体细胞变异） | ~26,000+（所有检测工具的未过滤联合结果） |
| pVACseq/seq2neo 兼容性 | ✓ | ✗（无 FORMAT 数据） |
| 过滤状态 | 已过滤（PASS + 软过滤） | 所有检测工具的所有变异 |

共识 VCF 是为数据标注设计的多检测工具聚合结果——所有变异证据存储在 INFO 字段（`VAF_BY_CALLER`、`DP_BY_CALLER` 等）而非每样本 FORMAT 字段。这对流水线的主要用途是正确的，但与需要每样本基因型数据的 neoantigen 工具不兼容。

### 各检测工具 FORMAT 字段差异

| 字段 | Mutect2 | DeepSomatic | Strelka2 |
|------|---------|-------------|----------|
| `GT` | ✓ | ✓ | ✗ |
| `AD`（ref,alt） | ✓ `Number=R` | ✓ `Number=R` | ✗（使用 `AU`/`CU`/`GU`/`TU` 或 `TAR`/`TIR`） |
| `AF` | ✓（字段名 `AF`） | ✗（字段名 `VAF`） | ✗（需从碱基计数推算） |
| `DP` | ✓ | ✓ | ✓ |
| `GQ` | ✓ | ✓ | ✗ |
| `F1R2`/`F2R1`/`SB` | ✓ | ✗ | ✗ |

`FORMAT_HARMONIZER` 在共识步骤前将 Strelka2 和 DeepSomatic 的 FORMAT 字段标准化为 Mutect2 格式。但共识 VCF 本身按设计没有每样本 FORMAT 数据——它无法对齐到 Mutect2 FORMAT，因为它代表多个检测工具的联合结果，而非单个检测工具的基因型调用。

### 何时使用 `consensus`

仅当下游工具能够处理基于 INFO 字段的格式，且确实需要多检测工具联合结果时，才使用 `neoantigen_input_source = 'consensus'`。带 `--neoantigen` 的共识 VCF 会添加 `AD_BY_CALLER` 和 `AF_BY_CALLER` INFO 字段，编码每个检测工具的等位基因深度和频率。

---

## 参数说明

### Neoantigen 专用参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `enable_neoantigen_workflow` | boolean | `false` | 启用 neoantigen 准备分支。 |
| `neoantigen_input_source` | string | `'mutect2'` | VCF 来源：`mutect2`（推荐）或 `consensus`。 |
| `salmon_index` | string | 共享索引 | 预构建的 Salmon v1.11.x SSHash 索引目录的绝对路径。 |
| `salmon_libtype` | string | `'A'` | 传递给 `salmon quant --libType` 的文库类型。`'A'` = 自动检测。 |
| `salmon_gc_bias` | boolean | `false` | 向 `salmon quant` 传递 `--gcBias` 以启用 GC 偏差校正（每样本约增加 2 分钟）。 |

### 预检验证

流水线在启动任何工作流步骤前进行验证：

- `enable_neoantigen_workflow = true` 且未设置 `salmon_index` → 报错退出。
- `salmon_index` 路径不存在或不是目录 → 报错退出。
- `neoantigen_input_source` 不是 `consensus` 或 `mutect2` → 报错退出。

---

## Salmon 索引

共享路径中已有预构建的 Gencode v49 索引：

```
/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/salmon/salmon_index_gencode_v49_salmon_v1.11
```

该路径已在 `neoantigen.shared.config` 中设为默认值，无需任何操作（除非需要不同的转录组）。

自行构建索引：

```bash
bash /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/scripts/build_salmon_index.sh \
    --outdir /path/to/output_dir
```

### 注意事项

- **Salmon v1.11.x 索引不兼容性**：v1.11.x 使用新的基于 SSHash 的 k-mer 索引。v1.10.x 或更早版本构建的索引**不兼容**，必须重新构建。
- 对 Gencode 转录组建索引时必须使用 `--gencode` 标志（Gencode FASTA 头部使用管道符分隔格式）。流水线的 `QUANT_TSV_NORMALIZE` 步骤会去除这些后缀，生成供下游工具使用的 `quant.tsv`。
- Gencode v49：GRCh38.p14，Ensembl 115，2025 年 9 月发布。

---

## 预期输出

```
${outdir}/neoantigen/
├── <DNA肿瘤样本ID>/                         # status=1 样本（如 COO8801DT_vs_COO8801DN）
│   ├── *.filtered.vcf.gz                   # neoantigen 就绪 VCF（默认为 Mutect2 过滤结果）
│   └── *.filtered.vcf.gz.tbi
│
└── <RNA肿瘤样本ID>/                          # status=2 样本（如 COO8801RT-LX）
    ├── quant.sf                             # Salmon 原始输出（Gencode 管道符分隔名称）
    ├── quant.tsv                            # 标准化转录本 ID（纯 ENST ID）
    ├── lib_format_counts.json
    ├── cmd_info.json
    └── aux_info/
        ├── meta_info.json
        ├── ambig_info.tsv
        └── fld.gz
```

`quant.tsv` 是供 pVACseq、seq2neo 等下游 neoantigen 工具使用的文件。

---

## run.sh 参数

```
--input  CSV     输入样本表 CSV（默认：共享 COO8801 测试数据集）
--outdir DIR     输出目录（默认：output/COO8801.neoantigen，相对路径）
--main-nf PATH   main.nf 路径（默认：通过 $BASH_SOURCE 从代码库解析）
--conf   PATH    配置文件路径（默认：同目录下的 neoantigen.shared.config）
--dry-run        打印 nextflow 命令但不执行
-h, --help       显示帮助
```

---

## 手动 nextflow 调用

> **关于 micromamba 的说明：** `micromamba run -n nextflow` 包装器仅在 nextflow 安装于 micromamba/conda 环境中时才需要。如果 nextflow 已在 `PATH` 中，可直接运行：
> ```bash
> nextflow run /path/to/main.nf -c /path/to/config ...
> ```

```bash
HTTPS_PROXY="http://10.233.17.241:3128" \
NXF_CONDA_CACHEDIR="/t9k/mnt/joey/nf_conda_envs" \
NXF_CONDA_USEMAMBA=true \
micromamba run -n nextflow nextflow run \
    /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/main.nf \
    -c /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/neoantigen/neoantigen.shared.config \
    --input  /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/input/test.rdv.shared.csv \
    --outdir /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/neoantigen_output/COO8801 \
    -offline -with-conda -resume
```
