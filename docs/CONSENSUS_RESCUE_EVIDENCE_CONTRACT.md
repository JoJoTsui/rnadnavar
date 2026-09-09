# Consensus and rescue evidence contract

Updated 2026-09-09 after the WES_LL runtime and logic review. These rules apply to both FASTQ-derived and caller-ready BAM-derived VCFs. They change only downstream evidence handling and decisions; mapping, callers, references, regions and resource configuration are unchanged.

## Eligible votes

A caller record is observed whenever present. It is eligible only when not Artifact, not an undecomposed multi-ALT record, and supported by valid allele-specific alternate depth meeting `min_alt_support`. With a positive floor, unavailable AD does not vote. Setting the floor to zero explicitly disables the depth requirement, but not the Artifact/multi-ALT exclusions. Sampleless consensus labels do not create independent caller votes; rescue evaluates those labels separately.

Depth/count measurements must be nonnegative integers. VAF must be finite and in [0,1]. Missing values, nonfinite numbers and negative integer missing sentinels are unavailable, never inferred as zero. Invalid fields are retained in detailed evidence under `INVALID_FIELDS`. Valid zero is distinct from unavailable. VAF may be derived from valid AD when reported AF is unavailable, with `VAF_SOURCE=derived`; Strelka tier-one SNV/indel counts follow the same rule. Aggregate statistics exclude unavailable/invalid measurements.

## Modality counts and RNA editing

Rescue emits `N_DNA_CALLERS_OBSERVED`, `N_DNA_CALLERS_ELIGIBLE`, and `N_DNA_CALLERS_SOMATIC`, and the corresponding RNA fields. Somatic counts are eligible callers classified Somatic. Consensus records are excluded. Legacy `N_*_CALLERS_SUPPORT`/`DNA_SUPPORT`/`RNA_SUPPORT` retain their observed-count meanings for compatibility.

RNA-editing annotation prefers the explicit eligible counts: an Artifact or insufficient-AD DNA record cannot alone prevent an RNA-only editing decision. Older VCFs lacking explicit counts use the legacy fallback; that fallback cannot reconstruct eligibility and should not be interpreted as the new contract. All supplied individual caller evidence should accompany rescue when these counts are required.

## Classification and verification

Within-modality thresholds, biological class vocabulary, and tie-to-Artifact behavior remain. `ENS_SUPPORT=k/n` describes eligible support in the configured caller panel; it is not a calibrated confidence interval. Baseline preservation remains opt-in and is not ordinary threshold consensus.

Rescue first computes its biological outcome, including Artifact veto and disagreement resolution. Optional DNA verification can then withhold only an otherwise-Somatic outcome that lacks independent DNA support. A passed DNA Somatic consensus or eligible DNA-only evidence satisfying the configured consensus rules establishes that independent support; one DNA vote in a cross-modality promotion does not automatically suffice. Inconclusive/rejected verification leaves such RNA-dependent outcomes NoConsensus with a rationale. Established Artifact/Germline/Reference labels are preserved. Non-Somatic final outcomes clear rescued/promotion flags. Omitting verification preserves the existing unverified workflow behavior.

## Lossless provenance

`EVIDENCE_SCHEMA=paired-v2` identifies `CALLER_EVIDENCE`, a Number=1 String containing percent-encoded JSON. Decode once with URL decoding, then parse JSON. The payload has `version` and `entries`; each entry identifies modality, caller, sample role, sample identity, normalized allele key and alignment round, plus measurements, availability and source provenance.

Exact observations coalesce deterministically with their source lists. Conflicting measurements remain distinct observations, never extra independent caller votes. Phased GT and comma-separated AD survive encoding. Unavailable normal evidence is explicit. Legacy evidence enriches missing detail; ambiguously concatenated legacy strings remain raw provenance instead of fabricated measurements. Legacy rounded VAF comparisons tolerate only their four-decimal serialization precision. Do not treat old summary INFO fields as the complete provenance source.

Unknown sample/modality/round is explicitly unknown. Rescue binds only unknown identity to supplied context: DNA first round and RNA `--alignment-round`; already known origins are not rewritten. Standalone consensus without explicit round context may retain unknown round until contextual rescue. Realignment is a correlated reassessment, not an independent biological sample or caller vote.

## Reporting

The old `rescue_rate = rescued / cross_modality` is removed: promotions could make it exceed 100%. `rescued_union_fraction = rescued / total_variants` has an explicit denominator and is an output composition statistic, not precision, recall, or rescue success.

Rescue reports DNA Somatic baseline, final Somatic count, retained DNA Somatic sites, newly Somatic sites versus DNA, and lost DNA Somatic sites. These compare labels at variant keys, not truth status. Report benchmark precision/recall/F1 and paired FP/FN transitions separately using a fixed truth set, regions and variant-type definition. No performance improvement over DeepSomatic is claimed by these correctness changes alone.

## Cache boundary and validation

No Nextflow alignment/caller module or configuration change is required. Preserve original caller files, indexes, work directory and run history. Only downstream tasks should be invalidated; actual cache hits must be checked on resume. New Python helper imports alone are not evidence that a previously successful task will be invalidated.

Validation covers numeric edge cases, no-AD votes, verification versus veto, count semantics, source conflicts and round trips, and replay using existing WES_LL caller VCFs. The original seq2neo helper tests were updated to their existing `_AUDITED` exports and `audit_fai` variable without modifying those production tasks.
