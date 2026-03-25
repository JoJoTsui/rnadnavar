# FASTQ 到 Rescued VCF：标签、工作流与 FILTER 解读

## 目的
本文档是追踪 rnadnavar 从原始 DNA/RNA FASTQ 输入到最终 rescued VCF 输出的唯一参考文档，FILTER 解读与当前流程实现保持一致。

范围：
- 输入标签与元数据传播
- 从比对到 rescue 及可选二次 rescue 的端到端工作流
- consensus 和 rescue 输出中 FILTER 值的赋值方式
- 使用一次真实运行输出进行验证：COO8801.shared

非范围：
- 不涉及算法变更
- 不涉及代码行为变更

## 1. 输入标签与样本标识

### 1.1 Samplesheet 标签
流程期望每行包含样本元数据，其中包含一个用于通道路由的模态状态字段。

观测到的运行输入：
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/params_2026-03-18_20-24-15.json](../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/params_2026-03-18_20-24-15.json)
  - input = /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/input/test.rdv.shared.csv

输入行：
- patient = COO8801, status = 0, sample = COO8801DN, lane = LX
- patient = COO8801, status = 1, sample = COO8801DT, lane = LX
- patient = COO8801, status = 2, sample = COO8801RT, lane = LX

工作流路由使用的解释：
- status 0：DNA 正常样本
- status 1：DNA 肿瘤样本
- status 2：RNA 肿瘤样本

### 1.2 元数据传播
在 [../subworkflows/local/samplesheet_to_channel/main.nf](../subworkflows/local/samplesheet_to_channel/main.nf) 中，通道元数据携带 patient、sample、lane、status 及派生 id。DNA 和 RNA 按 status 分流，RNA 具有 STAR 专用的 read-group 处理。

关键行为：
- DNA 路径：status <= 1
- RNA 路径：status == 2
- FASTQ 阶段的 id 命名规则：sample-lane（例如 COO8801DT-LX）

## 2. 从 FASTQ 到 Rescued VCF 的工作流路径

顶层编排位于 [../workflows/rnadnavar.nf](../workflows/rnadnavar.nf)。

### 2.1 阶段
1. 参考基因组/索引准备
2. 比对
   - DNA 通过 BWA 系列工具
   - RNA 通过 STAR
3. GATK 预处理
4. 变异检测（本次运行：DeepSomatic、Mutect2、Strelka）
5. 标准化
6. 模态内 consensus
7. 跨模态 rescue
8. 注释与过滤
9. 可选 RNA 重比对，然后进行二次 rescue
10. 报告（MultiQC、traces、DAG、timelines）

### 2.2 运行配置证据（COO8801.shared）
来自 [../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/params_2026-03-18_20-24-15.json](../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/params_2026-03-18_20-24-15.json)：
- rna = true
- dna = true
- tools = deepsomatic,mutect2,strelka,vep,norm,consensus,rescue,filtering,rna_filtering,realignment
- realignment_mode = vcf
- rescue_snv_thr = 2
- rescue_indel_thr = 2
- enable_rna_annotation = true
- enable_cosmic_gnomad_annotation = true

## 3. 具体输出文件（已验证）

以下所有示例均存在于附带的运行输出目录中。

### 3.1 Consensus 输出
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801DT_vs_COO8801DN/COO8801DT_vs_COO8801DN.consensus.vcf.gz](../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801DT_vs_COO8801DN/COO8801DT_vs_COO8801DN.consensus.vcf.gz)
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801RT_vs_COO8801DN/COO8801RT_vs_COO8801DN.consensus.vcf.gz](../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801RT_vs_COO8801DN/COO8801RT_vs_COO8801DN.consensus.vcf.gz)

### 3.2 第一次 rescue 输出
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN.rescued.vcf.gz](../sequencing/aim_exp/rdv_test/output/COO8801.shared/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN.rescued.vcf.gz)
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN.rescue.filtered.stripped.vep.vcf.gz](../sequencing/aim_exp/rdv_test/output/COO8801.shared/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN.rescue.filtered.stripped.vep.vcf.gz)

### 3.3 重比对 + 二次 rescue 输出
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/vcf_realignment/consensus/COO8801RT_realign_vs_COO8801DN/COO8801RT_realign_vs_COO8801DN.consensus.vcf.gz](../sequencing/aim_exp/rdv_test/output/COO8801.shared/vcf_realignment/consensus/COO8801RT_realign_vs_COO8801DN/COO8801RT_realign_vs_COO8801DN.consensus.vcf.gz)
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/vcf_realignment/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN.rescued.vcf.gz](../sequencing/aim_exp/rdv_test/output/COO8801.shared/vcf_realignment/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN.rescued.vcf.gz)

### 3.4 报告证据
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/reports/multiqc_report.html](../sequencing/aim_exp/rdv_test/output/COO8801.shared/reports/multiqc_report.html)
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/execution_trace_2026-03-18_20-23-45.txt](../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/execution_trace_2026-03-18_20-23-45.txt)
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/pipeline_dag_2026-03-18_20-23-45.html](../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/pipeline_dag_2026-03-18_20-23-45.html)

## 4. FILTER 分类：输出 VCF 中的含义

核心实现要点：
- 检测工具级别的映射位于 [../bin/vcf_utils/classification.py](../bin/vcf_utils/classification.py)
- Consensus 投票位于 [../bin/vcf_utils/variant_classifier_unified.py](../bin/vcf_utils/variant_classifier_unified.py)
- 最终输出 FILTER 赋值位于 [../bin/vcf_utils/io_utils.py](../bin/vcf_utils/io_utils.py)

### 4.1 输出 FILTER 是统一的生物学分类
对于 consensus/rescue 输出记录，FILTER 由统一分类赋值，而非直接复制各检测工具的原始 FILTER 字符串。

统一类别：
- Somatic（体细胞突变）
- Germline（胚系突变）
- Reference（参考等位基因）
- Artifact（技术假阳性）
- NoConsensus（无共识）
- RNAedit（注释阶段赋值）

### 4.2 原始检测工具过滤信息保留在 INFO 中
脚本将各检测工具的详细信息写入 INFO 字段，包括：
- FILTERS_ORIGINAL
- FILTERS_NORMALIZED
- FILTERS_CATEGORY
- UNIFIED_FILTER
- UNIFIED_FILTER_DNA
- UNIFIED_FILTER_RNA
- PASSES_CONSENSUS

实际解读：
- 使用 FILTER 字段作为该输出文件中的最终分类决策。
- 使用 INFO 字段追溯来源及各检测工具的详细说明。

## 5. Consensus 逻辑（模态内）

通过 [../bin/vcf_utils/variant_classifier_unified.py](../bin/vcf_utils/variant_classifier_unified.py) 中的 `UnifiedVariantClassifier.classify_consensus_variant` 实现。

规则摘要：
1. 使用各独立检测工具（排除 consensus 工具本身）。
2. 按变异类型检查阈值：
   - SNV 阈值
   - indel 阈值
3. 若检测工具数量低于阈值：NoConsensus
4. 若达到阈值：
   - 明确多数分类：赋予该类别
   - 最高类别并列：Artifact

注意：
- PASSES_CONSENSUS 在 INFO 中为信息性字段。
- FILTER 仍由统一分类函数赋值。

## 6. Rescue 逻辑（跨模态）

通过 [../bin/vcf_utils/variant_classifier_unified.py](../bin/vcf_utils/variant_classifier_unified.py) 中的 `UnifiedVariantClassifier.classify_rescue_variant` 实现，工作流连接位于：
- [../subworkflows/local/vcf_consensus_workflow/main.nf](../subworkflows/local/vcf_consensus_workflow/main.nf)
- [../subworkflows/local/second_rescue/main.nf](../subworkflows/local/second_rescue/main.nf)

规则摘要：
1. 若存在，解析 DNA_consensus 和 RNA_consensus 标签。
2. 若两个标签均存在：
   - 标签相同：使用该标签
   - 均为 Artifact：Artifact
   - 不一致且均非 Artifact：
     - 若两个模态均有足够的独立检测工具支持：Artifact
     - 否则选择有足够支持的模态
     - 否则 NoConsensus
   - 一个为 Artifact，另一个非 Artifact：
     - 若非 Artifact 模态满足支持要求，优先使用该模态
     - 否则 Artifact
3. 若仅存在一个模态的 consensus：使用该标签。
4. 若无 consensus 标签：
   - 需要跨模态支持
   - 否则 NoConsensus
   - 其他不一致模式映射为 Artifact

## 7. 常见检测工具 FILTER 输入到生物学类别的映射

### 7.1 DeepSomatic
- PASS 或未过滤：Somatic
- GERMLINE：Germline
- RefCall：Reference
- 其他标签：Artifact

### 7.2 Mutect2
- PASS 或未过滤：Somatic
- germline 或 haplotype：Germline
- panel_of_normals、contamination 或 possible_numt：Reference
- 其他标签：Artifact

### 7.3 Strelka
- PASS 或未过滤：Somatic
- NT 指示 het/hom 且正常样本深度足够：Germline
- NT 指示 ref 且正常样本深度足够：Reference
- 其他情况：Artifact

## 8. 真实输出中观察到的文件命名规则

COO8801.shared 中观察到的命名规范：
- Consensus：{tumor}_vs_{normal}.consensus.vcf.gz
- Rescue：{dna_pair}_rescued_{rna_pair}.rescued.vcf.gz
- 重比对 rescue：{dna_pair}_rescued_{rna_realign_pair}.rescued.vcf.gz
- rescue 后续阶段输出包含以下后缀：
  - .rescue.rna_annotated.vcf.gz
  - .rescue.cosmic_gnomad_annotated.*.vcf.gz
  - .rescue.filtered.stripped.vep.vcf.gz

## 9. 未来文档更新的验证清单

1. 确认所有声明的输出文件存在于真实运行目录中。
2. 确认所有声明的 FILTER 规则与当前分类器/写入器代码一致。
3. 确认措辞区分：
   - 统一输出 FILTER
   - INFO 中原始检测工具 FILTER 的来源信息
4. 确认 rescue 描述包含第一次 rescue 和可选的重比对二次 rescue。
5. 确认样本标签示例反映基于 status 的路由。

## 10. 面向不同受众的使用说明

面向开发者：
- 结合本文档与规则页面获取实现细节。

面向演示受众：
- 使用 docs/presentation 中的幻灯片获取简洁流程概述。
- 幻灯片有意将 FILTER 讨论简化为 Somatic、Germline 和 Reference；完整分类体系保留在本文档中。

---
验证基准：COO8801.shared 运行元数据及输出，与 pipeline_info 条目日期对齐，日期为 2026-03-18。
