# Rescue VCF Rules (Code-Accurate)

## Scope
This page documents cross-modality rescue behavior from:
- ../bin/run_rescue_vcf.py
- ../bin/vcf_utils/variant_classifier_unified.py
- ../bin/vcf_utils/io_utils.py
- ../subworkflows/local/vcf_consensus_workflow/main.nf
- ../subworkflows/local/second_rescue/main.nf

For full FASTQ to final output interpretation and run evidence, see:
- ./FASTQ_TO_RESCUED_VCF_FILTER_GUIDE.md

## Summary
Rescue mode merges DNA and RNA evidence and assigns unified biological FILTER classes with modality-aware logic.

Unified classes:
- Somatic
- Germline
- Reference
- Artifact
- NoConsensus

## Rescue inputs
Required in workflow wiring:
- DNA consensus VCF
- RNA consensus VCF

Recommended and used in this pipeline:
- Individual DNA caller VCFs
- Individual RNA caller VCFs

Individual caller support is important for modality-specific support counts and for parts of disagreement resolution.

## Modality tagging
Variant records are tagged with caller modality via a modality map and caller names such as:
- DNA_consensus
- RNA_consensus
- DNA_mutect2
- RNA_strelka

## Unified FILTER logic in rescue mode
Implemented in UnifiedVariantClassifier.classify_rescue_variant.

Decision outline:
1. Extract DNA and RNA consensus labels when present.
2. If both consensus labels exist:
   - same label -> use that label
   - both Artifact -> Artifact
   - both non-Artifact but different ->
     - if both modalities have enough individual caller support -> Artifact
     - else choose modality with sufficient support
     - else NoConsensus
   - one Artifact and one non-Artifact -> prefer non-Artifact label only when support threshold is met, else Artifact
3. If only one consensus label exists -> use that label.
4. If no consensus labels exist -> apply cross-modality support checks and disagreement rules; may return Artifact or NoConsensus.

Key implementation detail:
- Non-Artifact disagreement is not blindly resolved to DNA. It can return Artifact or modality-selected result depending on support counts.

## Output semantics
In ../bin/vcf_utils/io_utils.py:
- FILTER is assigned from unified rescue classification.
- INFO keeps provenance:
  - FILTERS_ORIGINAL
  - FILTERS_NORMALIZED
  - FILTERS_CATEGORY
  - UNIFIED_FILTER
  - UNIFIED_FILTER_DNA
  - UNIFIED_FILTER_RNA
  - PASSES_CONSENSUS
  - RESCUED

## First rescue and second rescue in workflow
- First rescue (consensus stage): ../subworkflows/local/vcf_consensus_workflow/main.nf
- Second rescue with realigned RNA: ../subworkflows/local/second_rescue/main.nf

## Verified naming examples from COO8801.shared

First rescue:
- ../sequencing/aim_exp/rdv_test/output/COO8801.shared/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN.rescued.vcf.gz

Second rescue after realignment:
- ../sequencing/aim_exp/rdv_test/output/COO8801.shared/vcf_realignment/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN.rescued.vcf.gz

Common downstream rescue files:
- .rescue.rna_annotated.vcf.gz
- .rescue.cosmic_gnomad_annotated.*.vcf.gz
- .rescue.filtered.stripped.vep.vcf.gz