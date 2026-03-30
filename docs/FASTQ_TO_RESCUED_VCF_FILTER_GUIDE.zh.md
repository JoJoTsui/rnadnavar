# FASTQ 到 Rescued VCF：标签、工作流与 FILTER 解读

## 目的

本文档追踪 rnadnavar 流程从原始 DNA/RNA FASTQ 输入到最终 rescued VCF 输出的完整过程，涵盖样本如何标记、流程如何逐步处理，以及如何解读输出文件中的 FILTER 值。

范围：
- 输入标签与元数据传播
- 从比对到 rescue 及可选二次 rescue 的端到端工作流
- consensus 和 rescue 输出中 FILTER 值的赋值方式
- 基于真实运行验证：COO8801.shared（2026-03-18）

不涉及：算法变更或代码行为变更。

---

## 1. 流程总览

该流程获取同一患者的 DNA 和 RNA 测序数据，在各自模态中独立检测变异，再跨模态整合证据，最终产出高可信度的变异集合。

```
DNA FASTQs ──► BWA 比对 ──► GATK 预处理 ──► 变异检测  ──┐
                                                        ├──► 模态内
RNA FASTQs ──► STAR 比对 ──► GATK 预处理 ──► 变异检测 ──┘   Consensus
                                                                  │
                                          ┌───────────────────────┘
                                          ▼
                                    跨模态 Rescue
                                (DNA consensus + RNA consensus)
                                          │
                                          ▼
                                    注释与过滤
                                          │
                          ┌───────────────┴───────────────┐
                          ▼                               ▼
                  [最终 rescued VCF]          可选：RNA 重比对
                                                          │
                                                          ▼
                                               二次 Rescue + 过滤
                                                          │
                                                          ▼
                                          [最终重比对 rescued VCF]
```

---

## 2. 输入标签与样本标识

### 2.1 Samplesheet 标签

Samplesheet 中每行定义一个样本，`status` 字段控制流程的路由方式。

运行输入（COO8801.shared）：

| patient  | status | sample    | lane | 含义           |
|----------|--------|-----------|------|----------------|
| COO8801  | 0      | COO8801DN | LX   | DNA 正常样本   |
| COO8801  | 1      | COO8801DT | LX   | DNA 肿瘤样本   |
| COO8801  | 2      | COO8801RT | LX   | RNA 肿瘤样本   |

路由规则：
- `status <= 1` → DNA 路径（BWA 比对）
- `status == 2` → RNA 路径（STAR 比对）

### 2.2 元数据传播

通道元数据（patient、sample、lane、status、id）贯穿所有阶段。FASTQ 阶段的 `id` 命名规则为 `sample-lane`，例如 `COO8801DT-LX`。

实现：`subworkflows/local/samplesheet_to_channel/main.nf`

---

## 3. 运行配置（COO8801.shared）

本次运行使用的关键参数：

| 参数                             | 值                                                           |
|----------------------------------|--------------------------------------------------------------|
| rna                              | true                                                         |
| dna                              | true                                                         |
| tools                            | deepsomatic, mutect2, strelka, vep, norm, consensus, rescue, filtering, rna_filtering, realignment |
| realignment_mode                 | vcf                                                          |
| rescue_snv_thr                   | 2                                                            |
| rescue_indel_thr                 | 2                                                            |
| enable_rna_annotation            | true                                                         |
| enable_cosmic_gnomad_annotation  | true                                                         |

来源：`pipeline_info/params_2026-03-18_20-24-15.json`

---

## 4. 工作流阶段

### 阶段 1 — 比对

- DNA reads → BWA-MEM2 → 排序后的 BAM
- RNA reads → STAR（两轮比对）→ 排序后的 BAM

### 阶段 2 — GATK 预处理

DNA 和 RNA BAM 均经过：
1. 标记重复序列
2. 碱基质量分数重校正（BQSR）
3. 仅 RNA：SplitNCigarReads（处理剪接位点）

### 阶段 3 — 变异检测

三个检测工具在每对肿瘤/正常样本上独立运行：

| 检测工具    | DNA | RNA |
|-------------|-----|-----|
| DeepSomatic | ✓   | ✓   |
| Mutect2     | ✓   | ✓   |
| Strelka     | ✓   | ✓   |

### 阶段 4 — 标准化

各检测工具的原始 VCF 在进入 consensus 前，使用 bcftools/vt 进行标准化（左对齐、分解复合变异）。

### 阶段 5 — 模态内 Consensus

在同一模态（DNA 或 RNA）内，获得足够多检测工具支持的变异被保留。详见第 6 节。

### 阶段 6 — 跨模态 Rescue

在一个模态中检测到的变异，可借助另一模态的证据被"救回"。详见第 7 节。

### 阶段 7 — 注释与过滤

Rescued VCF 经 VEP、RNA 编辑数据库、COSMIC 和 gnomAD 注释后过滤，产出最终输出。

### 阶段 8 — 可选 RNA 重比对 + 二次 Rescue

若启用 `realignment`，候选变异附近的 RNA reads 将被重新比对以提升灵敏度，随后对重比对后的 RNA consensus 进行第二轮 rescue。

---

## 5. 输出文件

所有路径均相对于运行输出目录（`COO8801.shared/`）。

### 5.1 模态内 consensus

```
consensus/COO8801DT_vs_COO8801DN/COO8801DT_vs_COO8801DN.consensus.vcf.gz
consensus/COO8801RT_vs_COO8801DN/COO8801RT_vs_COO8801DN.consensus.vcf.gz
```

### 5.2 第一次 rescue

```
rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/
  ├── ...rescued.vcf.gz
  └── ...rescue.filtered.stripped.vep.vcf.gz
```

### 5.3 重比对 + 二次 rescue

```
vcf_realignment/consensus/COO8801RT_realign_vs_COO8801DN/
  └── ...consensus.vcf.gz

vcf_realignment/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN/
  └── ...rescued.vcf.gz
```

### 5.4 Rescue 后注释阶段

每个 rescue VCF 依次经过注释步骤，产生带有以下后缀的中间文件：

```
.rescued.vcf.gz
.rescue.rna_annotated.vcf.gz
.rescue.cosmic_gnomad_annotated.*.vcf.gz
.rescue.filtered.stripped.vep.vcf.gz        ← 最终输出
```

### 5.5 报告

```
reports/multiqc_report.html
pipeline_info/execution_trace_2026-03-18_20-23-45.txt
pipeline_info/pipeline_dag_2026-03-18_20-23-45.html
```

---

## 6. 模态内 Consensus 逻辑

**目标：** 仅保留在同一模态（DNA 或 RNA）内获得足够多独立检测工具支持的变异。

**工作方式：**

```
对于某一模态中的每个变异：
  count = 独立检测工具（不含 consensus 工具本身）的支持数量
  threshold = rescue_snv_thr（SNV）或 rescue_indel_thr（indel）

  若 count < threshold：
      → NoConsensus

  若 count >= threshold：
      统计各检测工具赋予的生物学类别
      若某一类别占明确多数 → 赋予该类别
      若最高类别并列       → Artifact
```

INFO 字段 `PASSES_CONSENSUS` 记录是否达到阈值，但最终 `FILTER` 值始终由统一分类函数赋值。

实现：`bin/vcf_utils/variant_classifier_unified.py` 中的 `UnifiedVariantClassifier.classify_consensus_variant`

---

## 7. 跨模态 Rescue 逻辑

**目标：** 利用另一模态的证据，找回被某一模态遗漏的变异。仅在 RNA 中检测到的变异可通过 DNA 证据得到确认，反之亦然。

**工作方式：**

```
给定一个携带 DNA_consensus 标签和/或 RNA_consensus 标签的变异：

情况 1 — 两个模态均有 consensus 标签：
  ├── 两侧标签相同                    → 使用该标签
  ├── 两侧均为 Artifact               → Artifact
  ├── 两侧均非 Artifact 但不一致：
  │     ├── 两侧均有足够检测工具支持  → Artifact（冲突无法解决）
  │     ├── 仅一侧有足够支持          → 使用该侧标签
  │     └── 两侧均无足够支持          → NoConsensus
  └── 一侧为 Artifact，另一侧非 Artifact：
        ├── 非 Artifact 侧满足支持要求 → 使用非 Artifact 标签
        └── 否则                       → Artifact

情况 2 — 仅一个模态有 consensus 标签：
  → 直接使用该标签

情况 3 — 无任何 consensus 标签：
  ├── 存在跨模态检测工具支持 → 基于该支持进行分类
  └── 否则                   → NoConsensus
```

实现：`bin/vcf_utils/variant_classifier_unified.py` 中的 `UnifiedVariantClassifier.classify_rescue_variant`

工作流连接：`subworkflows/local/vcf_consensus_workflow/main.nf` 和 `subworkflows/local/second_rescue/main.nf`

---

## 8. 输出 VCF 中的 FILTER 值

### 8.1 FILTER 的含义

在 consensus 和 rescue 输出文件中，`FILTER` 列包含统一的生物学分类，而非任何单个检测工具的原始过滤字符串。

| FILTER 值     | 含义                                               |
|---------------|----------------------------------------------------|
| Somatic       | 肿瘤特异性变异，正常样本中不存在                   |
| Germline      | 肿瘤和正常样本中均存在的变异（遗传性）             |
| Reference     | 被判定为参考等位基因或可能的测序假阳性             |
| Artifact      | 技术假阳性，或检测工具间存在无法解决的分歧         |
| NoConsensus   | 检测工具支持不足，无法做出判断                     |
| RNAedit       | 已知 RNA 编辑位点（在注释阶段赋值）                |

### 8.2 原始检测工具过滤信息保留在 INFO 中

各检测工具的原始信息不会丢失，而是存储在 INFO 字段中以供溯源：

| INFO 字段           | 内容                                         |
|---------------------|----------------------------------------------|
| FILTERS_ORIGINAL    | 各检测工具的原始 FILTER 字符串               |
| FILTERS_NORMALIZED  | 标准化后的过滤字符串                         |
| FILTERS_CATEGORY    | 各检测工具的生物学类别                       |
| UNIFIED_FILTER      | 最终统一分类                                 |
| UNIFIED_FILTER_DNA  | 仅基于 DNA 检测工具的统一分类               |
| UNIFIED_FILTER_RNA  | 仅基于 RNA 检测工具的统一分类               |
| PASSES_CONSENSUS    | 是否达到 consensus 阈值                      |

**使用原则：** 用 `FILTER` 获取最终判断；用 `INFO` 字段了解判断依据。

实现：`bin/vcf_utils/classification.py`、`bin/vcf_utils/variant_classifier_unified.py`、`bin/vcf_utils/io_utils.py`

---

## 9. 各检测工具 FILTER 到生物学类别的映射

### 9.1 DeepSomatic

| 检测工具 FILTER | 生物学类别 |
|-----------------|------------|
| PASS / （无）   | Somatic    |
| GERMLINE        | Germline   |
| RefCall         | Reference  |
| 其他            | Artifact   |

### 9.2 Mutect2

| 检测工具 FILTER                                    | 生物学类别 |
|----------------------------------------------------|------------|
| PASS / （无）                                      | Somatic    |
| germline、haplotype                                | Germline   |
| panel_of_normals、contamination、possible_numt     | Reference  |
| 其他                                               | Artifact   |

### 9.3 Strelka

| 条件                                   | 生物学类别 |
|----------------------------------------|------------|
| PASS / （无）                          | Somatic    |
| NT = het 或 hom，正常样本深度足够      | Germline   |
| NT = ref，正常样本深度足够             | Reference  |
| 其他情况                               | Artifact   |

---

## 10. 文件命名规则

| 输出类型              | 命名模式                                                   |
|-----------------------|------------------------------------------------------------|
| 模态内 consensus      | `{tumor}_vs_{normal}.consensus.vcf.gz`                     |
| 第一次 rescue         | `{dna_pair}_rescued_{rna_pair}.rescued.vcf.gz`             |
| 重比对 rescue         | `{dna_pair}_rescued_{rna_realign_pair}.rescued.vcf.gz`     |
| 最终过滤输出          | `...rescue.filtered.stripped.vep.vcf.gz`                   |

---

## 11. 验证清单

更新本文档时，请确认：

1. 所有列出的输出文件存在于真实运行目录中。
2. 所有 FILTER 规则与当前分类器和写入器代码一致。
3. 措辞明确区分统一输出 `FILTER` 与 `INFO` 中各检测工具的来源信息。
4. Rescue 描述涵盖第一次 rescue 和可选的重比对二次 rescue。
5. 样本标签示例反映基于 status 的路由。

---

## 12. 面向不同受众的使用说明

面向开发者：结合本文档与分类器源文件获取实现细节。

面向演示受众：`docs/presentation` 中的幻灯片提供简化流程概述。幻灯片有意将 FILTER 类别简化为 Somatic、Germline 和 Reference；完整分类体系保留在本文档中。

---

验证基准：COO8801.shared 运行，pipeline_info 条目日期为 2026-03-18。
