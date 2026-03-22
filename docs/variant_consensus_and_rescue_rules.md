# Variant Consensus and Cross-Modality Rescue Rules

## Purpose
This page is a synchronized rules summary for consensus and rescue classification.
It is aligned to current code behavior and to observed output naming in COO8801.shared.

Primary implementation references:
- ../bin/vcf_utils/classification.py
- ../bin/vcf_utils/variant_classifier_unified.py
- ../bin/vcf_utils/io_utils.py
- ../subworkflows/local/vcf_consensus_workflow/main.nf
- ../subworkflows/local/second_rescue/main.nf

Primary narrative reference:
- ./FASTQ_TO_RESCUED_VCF_FILTER_GUIDE.md

## Canonical biological classes
- Somatic
- Germline
- Reference
- Artifact
- NoConsensus
- RNAedit (annotation stage)

## Consensus mode rules
Consensus mode is used for within-modality aggregation.

Input: caller VCFs for one modality.

Decision rules:
1. Compute caller-level biological classifications.
2. Count supporting individual callers by variant type.
3. If below threshold -> NoConsensus.
4. If threshold met -> majority class.
5. If tie at top -> Artifact.

Important:
- Tie is disagreement, not priority promotion.
- FILTER in output is unified class from classifier logic.
- PASSES_CONSENSUS is informational in INFO.

## Rescue mode rules
Rescue mode combines DNA and RNA consensus plus individual caller support.

Decision rules (high level):
1. Use DNA_consensus and RNA_consensus labels when available.
2. Agreement between modalities keeps agreed class.
3. Disagreement between non-Artifact classes can map to Artifact, or to a modality-specific class if support thresholds justify it.
4. Missing consensus labels relies on cross-modality support and disagreement checks.
5. Insufficient evidence yields NoConsensus.

Important:
- Rescue logic uses support-aware conditions, not a fixed DNA-only override.

## FILTER and INFO semantics
In consensus/rescue output VCF:
- FILTER = unified biological class for final interpretation.
- Original caller filters and normalized caller classes are retained in INFO:
  - FILTERS_ORIGINAL
  - FILTERS_NORMALIZED
  - FILTERS_CATEGORY
  - UNIFIED_FILTER
  - UNIFIED_FILTER_DNA
  - UNIFIED_FILTER_RNA
  - PASSES_CONSENSUS
  - RESCUED

## Verified output naming examples
From COO8801.shared:
- Consensus:
  - ../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801DT_vs_COO8801DN/COO8801DT_vs_COO8801DN.consensus.vcf.gz
- First rescue:
  - ../sequencing/aim_exp/rdv_test/output/COO8801.shared/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN.rescued.vcf.gz
- Second rescue (realigned RNA):
  - ../sequencing/aim_exp/rdv_test/output/COO8801.shared/vcf_realignment/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN.rescued.vcf.gz

## Run-verified configuration snapshot
From ../sequencing/aim_exp/rdv_test/output/COO8801.shared/pipeline_info/params_2026-03-18_20-24-15.json:
- rna = true
- dna = true
- tools include consensus,rescue,realignment,vep,filtering,rna_filtering
- rescue_snv_thr = 2
- rescue_indel_thr = 2
- realignment_mode = vcf
