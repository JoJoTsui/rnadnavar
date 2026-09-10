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
The [evidence contract](CONSENSUS_RESCUE_EVIDENCE_CONTRACT.md) defines current missing-value, eligible-vote and paired-v2 provenance rules. With a positive alternate-read floor, missing/invalid AD cannot vote; a measured zero is not a missing value.

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

## Caller-support INFO field
Written by `write_union_vcf` in consensus mode only; rescue output is unchanged.
Every consensus record carries `ENS_SUPPORT`, the descriptive caller-detection
support fraction `k/n` (for example, `2/3`). It is not a calibrated probability
or confidence interval for the final biological label.

Semantics:
- k = number of supporting callers (non-Artifact record clearing the min
  alt-read floor — the same support set behind N_SUPPORT_CALLERS).
- n = the explicit expected caller panel for the sample and modality, not the
  VCF files that happened to be discovered. A caller absent at a site counts
  as a non-support vote.
- Missing, duplicate, and unexpected caller VCFs fail the consensus invocation.
- `ENS_CONF_LO` and `ENS_CONF_HI` were removed by ADR-0003 because the fixed,
  correlated three-caller panel does not satisfy the binomial interpretation.

## Verified naming examples from COO8801.shared
- ../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801DT_vs_COO8801DN/COO8801DT_vs_COO8801DN.consensus.vcf.gz
- ../sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801RT_vs_COO8801DN/COO8801RT_vs_COO8801DN.consensus.vcf.gz

## Opt-in caller-aware baseline policy

`--preserve-baseline-callers` is an experiment-only option. A named caller's Somatic record is retained only when its record is in the eligible support set (non-Artifact and satisfying available tumor AD floor). The legacy threshold policy remains the default. The output rationale identifies `rule:preserve_baseline`; this does not make a blanket quality claim about every baseline call.

### Native-evidence SNV policy (default-on)

The consensus process supports `--native-evidence-snv` (also exposed as
`params.native_evidence_snv`, default `true`). Pass the explicit opt-out value
when reproducing the legacy threshold policy. When enabled, SNVs may be
classified Somatic by native caller evidence rather than the ordinary caller
vote alone:

- a qualified DeepSomatic Somatic record is retained; or
- a candidate locus with a DeepSomatic record that is not accepted as a
  DeepSomatic Somatic call is admitted only when that record has QUAL > 0 and
  Mutect2 TLOD is at least 12 and GERMQ is at least 60, unless Mutect2 has one
  of the explicit contamination/germline artifact combinations; a locus with
  no DeepSomatic record is not admitted by this rule.

This policy is intentionally SNV-only. Indels are not restricted to DeepSomatic:
indels continue to use the configured ordinary consensus threshold, but they
are not promoted by this native-evidence rule. The policy is default-on for new runs; the explicit opt-out preserves legacy
behavior for controlled rollback and comparison. Caller and alignment caches
remain reusable because the policy is consumed after caller VCF production.
