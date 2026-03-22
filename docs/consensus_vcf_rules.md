# Consensus VCF Rules (Code-Accurate)

## Scope
This page documents current within-modality consensus behavior implemented by:
- ../bin/run_consensus_vcf.py
- ../bin/vcf_utils/aggregation.py
- ../bin/vcf_utils/classification.py
- ../bin/vcf_utils/variant_classifier_unified.py
- ../bin/vcf_utils/io_utils.py

For end-to-end FASTQ to rescue interpretation with real output evidence, see:
- ./FASTQ_TO_RESCUED_VCF_FILTER_GUIDE.md

## Summary
Consensus mode aggregates per-caller VCF records from one modality and writes a union VCF where each output record gets a unified biological FILTER class.

Unified classes:
- Somatic
- Germline
- Reference
- Artifact
- NoConsensus

## Input and aggregation
1. Collect VCF files from input directory.
2. Parse variants and classify each caller record to biological classes.
3. Aggregate by normalized variant key chrom:pos:ref:alt.
4. Mark passes_consensus from support counts against thresholds.

Threshold inputs:
- snv_thr for SNVs
- indel_thr for indels

## Caller-level biological mapping
Implemented in ../bin/vcf_utils/classification.py.

- DeepSomatic:
  - PASS or unfiltered -> Somatic
  - GERMLINE -> Germline
  - RefCall -> Reference
  - others -> Artifact
- Mutect2:
  - PASS or unfiltered -> Somatic
  - germline or haplotype -> Germline
  - panel_of_normals or contamination or possible_numt -> Reference
  - others -> Artifact
- Strelka:
  - PASS or unfiltered -> Somatic
  - NT het/hom with sufficient normal depth -> Germline
  - NT ref with sufficient normal depth -> Reference
  - others -> Artifact

## Unified FILTER logic in consensus mode
Implemented in UnifiedVariantClassifier.classify_consensus_variant.

Rules:
1. Build list of individual caller classifications (exclude any _consensus callers).
2. If caller count < threshold for variant type -> NoConsensus.
3. Else compute majority class among caller classifications.
4. If clear majority -> return that class.
5. If tie for top class -> Artifact.

Important:
- Tie does not use priority to force Somatic/Germline/Reference.
- Tie is treated as disagreement and mapped to Artifact.

## Output semantics
Implemented in ../bin/vcf_utils/io_utils.py.

- FILTER field is set to unified biological class.
- Original per-caller filter strings are preserved in INFO fields:
  - FILTERS_ORIGINAL
  - FILTERS_NORMALIZED
  - FILTERS_CATEGORY
- PASSES_CONSENSUS is informational in INFO and does not independently override FILTER.

## Verified naming examples from COO8801.shared
- ../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801DT_vs_COO8801DN/COO8801DT_vs_COO8801DN.consensus.vcf.gz
- ../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801RT_vs_COO8801DN/COO8801RT_vs_COO8801DN.consensus.vcf.gz