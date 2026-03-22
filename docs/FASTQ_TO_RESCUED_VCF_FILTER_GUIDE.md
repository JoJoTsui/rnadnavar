# FASTQ to Rescued VCF: Labeling, Workflow, and FILTER Interpretation

## Purpose
This is the single source document for tracing rnadnavar from raw DNA/RNA FASTQ inputs to final rescued VCF outputs, with FILTER interpretation aligned to the current pipeline implementation.

Scope:
- Input labeling and metadata propagation
- End-to-end workflow from mapping to rescue and optional second rescue
- How FILTER values are assigned in consensus and rescue outputs
- Verification using one real run output: COO8801.shared

Non-scope:
- No algorithm changes
- No code behavior changes

## 1. Input Labeling and Sample Identity

### 1.1 Samplesheet labels
The pipeline expects per-row sample metadata, including a modality status field used throughout channel routing.

Observed run input:
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/params_2026-03-18_20-24-15.json](../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/params_2026-03-18_20-24-15.json)
  - input = /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/input/test.rdv.shared.csv

Input rows:
- patient = COO8801, status = 0, sample = COO8801DN, lane = LX
- patient = COO8801, status = 1, sample = COO8801DT, lane = LX
- patient = COO8801, status = 2, sample = COO8801RT, lane = LX

Interpretation used by workflow routing:
- status 0: DNA normal
- status 1: DNA tumor
- status 2: RNA tumor

### 1.2 Metadata propagation
In [../subworkflows/local/samplesheet_to_channel/main.nf](../subworkflows/local/samplesheet_to_channel/main.nf), channel metadata carries patient, sample, lane, status, and derived id. DNA and RNA diverge by status, and RNA has STAR-specific read-group handling.

Key behavior:
- DNA path: status <= 1
- RNA path: status == 2
- id convention at FASTQ stage: sample-lane (for example COO8801DT-LX)

## 2. Workflow Path from FASTQ to Rescued VCF

Top-level orchestration is in [../workflows/rnadnavar.nf](../workflows/rnadnavar.nf).

### 2.1 Stages
1. Reference/index preparation
2. Alignment
   - DNA via BWA family
   - RNA via STAR
3. GATK preprocessing
4. Variant calling (DeepSomatic, Mutect2, Strelka in this run)
5. Normalization
6. Within-modality consensus
7. Cross-modality rescue
8. Annotation and filtering
9. Optional RNA realignment, then second rescue
10. Reporting (MultiQC, traces, DAG, timelines)

### 2.2 Run configuration evidence (COO8801.shared)
From [../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/params_2026-03-18_20-24-15.json](../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/params_2026-03-18_20-24-15.json):
- rna = true
- dna = true
- tools = deepsomatic,mutect2,strelka,vep,norm,consensus,rescue,filtering,rna_filtering,realignment
- realignment_mode = vcf
- rescue_snv_thr = 2
- rescue_indel_thr = 2
- enable_rna_annotation = true
- enable_cosmic_gnomad_annotation = true

## 3. Concrete Output Artifacts (Verified)

All examples below exist in the attached run output tree.

### 3.1 Consensus outputs
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801DT_vs_COO8801DN/COO8801DT_vs_COO8801DN.consensus.vcf.gz](../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801DT_vs_COO8801DN/COO8801DT_vs_COO8801DN.consensus.vcf.gz)
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801RT_vs_COO8801DN/COO8801RT_vs_COO8801DN.consensus.vcf.gz](../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801RT_vs_COO8801DN/COO8801RT_vs_COO8801DN.consensus.vcf.gz)

### 3.2 First rescue outputs
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN.rescued.vcf.gz](../sequencing/aim_exp/rdv_test/output/COO8801.shared/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN.rescued.vcf.gz)
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN.rescue.filtered.stripped.vep.vcf.gz](../sequencing/aim_exp/rdv_test/output/COO8801.shared/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN.rescue.filtered.stripped.vep.vcf.gz)

### 3.3 Realignment + second rescue outputs
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/vcf_realignment/consensus/COO8801RT_realign_vs_COO8801DN/COO8801RT_realign_vs_COO8801DN.consensus.vcf.gz](../sequencing/aim_exp/rdv_test/output/COO8801.shared/vcf_realignment/consensus/COO8801RT_realign_vs_COO8801DN/COO8801RT_realign_vs_COO8801DN.consensus.vcf.gz)
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/vcf_realignment/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN.rescued.vcf.gz](../sequencing/aim_exp/rdv_test/output/COO8801.shared/vcf_realignment/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN.rescued.vcf.gz)

### 3.4 Reporting evidence
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/reports/multiqc_report.html](../sequencing/aim_exp/rdv_test/output/COO8801.shared/reports/multiqc_report.html)
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/execution_trace_2026-03-18_20-23-45.txt](../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/execution_trace_2026-03-18_20-23-45.txt)
- [../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/pipeline_dag_2026-03-18_20-23-45.html](../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/pipeline_dag_2026-03-18_20-23-45.html)

## 4. FILTER Classification: What It Means in Output VCF

Core implementation points:
- Caller-level mapping is in [../bin/vcf_utils/classification.py](../bin/vcf_utils/classification.py)
- Consensus voting is in [../bin/vcf_utils/variant_classifier_unified.py](../bin/vcf_utils/variant_classifier_unified.py)
- Final output FILTER assignment is in [../bin/vcf_utils/io_utils.py](../bin/vcf_utils/io_utils.py)

### 4.1 Output FILTER is unified biological classification
For consensus/rescue output records, FILTER is assigned from unified classification, not copied directly from per-caller raw FILTER strings.

Unified categories:
- Somatic
- Germline
- Reference
- Artifact
- NoConsensus
- RNAedit (annotation stage)

### 4.2 Original caller filters are preserved in INFO
The script writes caller-specific details into INFO fields such as:
- FILTERS_ORIGINAL
- FILTERS_NORMALIZED
- FILTERS_CATEGORY
- UNIFIED_FILTER
- UNIFIED_FILTER_DNA
- UNIFIED_FILTER_RNA
- PASSES_CONSENSUS

Practical interpretation:
- Use FILTER for final class decision in that output file.
- Use INFO fields for provenance and caller-by-caller explanation.

## 5. Consensus Logic (Within Modality)

Implemented through UnifiedVariantClassifier.classify_consensus_variant in [../bin/vcf_utils/variant_classifier_unified.py](../bin/vcf_utils/variant_classifier_unified.py).

Rule summary:
1. Use individual callers (exclude consensus callers).
2. Check threshold by variant type:
   - SNV threshold
   - indel threshold
3. If caller count is below threshold: NoConsensus
4. If threshold is met:
   - clear majority classification: assign that class
   - tie among top classes: Artifact

Notes:
- PASSES_CONSENSUS is informational in INFO.
- FILTER is still set by unified classification function.

## 6. Rescue Logic (Cross Modality)

Implemented through UnifiedVariantClassifier.classify_rescue_variant in [../bin/vcf_utils/variant_classifier_unified.py](../bin/vcf_utils/variant_classifier_unified.py), with workflow wiring in:
- [../subworkflows/local/vcf_consensus_workflow/main.nf](../subworkflows/local/vcf_consensus_workflow/main.nf)
- [../subworkflows/local/second_rescue/main.nf](../subworkflows/local/second_rescue/main.nf)

Rule summary:
1. Parse DNA_consensus and RNA_consensus labels if present.
2. If both labels exist:
   - same label: use it
   - both Artifact: Artifact
   - disagree and both non-Artifact:
     - if both modalities have enough supporting individual callers: Artifact
     - else pick modality with enough support
     - else NoConsensus
   - one Artifact and one non-Artifact:
     - prefer non-Artifact modality if it meets support requirement
     - otherwise Artifact
3. If only one modality consensus exists: use that label.
4. If no consensus labels:
   - require cross-modality support
   - otherwise NoConsensus
   - additional disagreement patterns map to Artifact

## 7. Mapping of Common Caller FILTER Inputs to Biological Classes

### 7.1 DeepSomatic
- PASS or unfiltered: Somatic
- GERMLINE: Germline
- RefCall: Reference
- other labels: Artifact

### 7.2 Mutect2
- PASS or unfiltered: Somatic
- germline or haplotype: Germline
- panel_of_normals or contamination or possible_numt: Reference
- other labels: Artifact

### 7.3 Strelka
- PASS or unfiltered: Somatic
- NT indicates het/hom with sufficient normal depth: Germline
- NT indicates ref with sufficient normal depth: Reference
- otherwise: Artifact

## 8. File Naming Patterns Seen in Real Outputs

Observed conventions in COO8801.shared:
- Consensus: {tumor}_vs_{normal}.consensus.vcf.gz
- Rescue: {dna_pair}_rescued_{rna_pair}.rescued.vcf.gz
- Realignment rescue: {dna_pair}_rescued_{rna_realign_pair}.rescued.vcf.gz
- Post-rescue staged outputs include suffixes:
  - .rescue.rna_annotated.vcf.gz
  - .rescue.cosmic_gnomad_annotated.*.vcf.gz
  - .rescue.filtered.stripped.vep.vcf.gz

## 9. Verification Checklist for Future Documentation Updates

1. Confirm all claimed outputs exist in a real run directory.
2. Confirm all claimed FILTER rules map to current classifier/writer code.
3. Confirm wording distinguishes:
   - unified output FILTER
   - original caller FILTER provenance in INFO
4. Confirm rescue descriptions include first rescue and optional realignment second rescue.
5. Confirm sample labeling examples reflect status-based routing.

## 10. Audience-Specific Usage

For developers:
- Use this document plus the rules pages for implementation details.

For presentation audiences:
- Use the slide deck at docs/presentation for concise flow.
- Slides intentionally simplify FILTER discussion to Somatic, Germline, and Reference; full taxonomy remains here.

---
Verification baseline: COO8801.shared run metadata and outputs, date-aligned with pipeline_info entries on 2026-03-18.
