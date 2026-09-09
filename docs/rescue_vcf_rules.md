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
The [evidence contract](CONSENSUS_RESCUE_EVIDENCE_CONTRACT.md) specifies observed/eligible/Somatic counts, lossless paired-v2 provenance, verification precedence and replacement rescue metrics. Verification is applied after biological classification and cannot erase established negative labels.

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

A "NoConsensus" consensus label means the modality reached no consensus for
the site and is treated as absent (it never blocks the rules below).

Decision outline:
1. Extract DNA and RNA consensus labels when present.
2. If both consensus labels exist:
   - same label -> use that label
   - both Artifact -> Artifact
   - both non-Artifact but different ->
     - if both modalities have enough individual caller support -> Artifact
     - else choose modality with sufficient support
     - else NoConsensus
   - one Artifact and one non-Artifact -> Artifact veto (see below), then
     prefer the non-Artifact label only when its support threshold is met,
     else Artifact
3. If only one consensus label exists -> use that label.
4. If no consensus labels exist -> apply cross-modality support checks,
   promotion, and disagreement rules (see below).

Key implementation detail:
- Non-Artifact disagreement is not blindly resolved to DNA. It can return Artifact or modality-selected result depending on support counts.

### Cross-modality promotion (audit M1)
When neither modality has a consensus label (the exact case rescue exists
for — a variant that failed within-modality consensus), individual callers
can rescue the site:

- If at least `rescue_promotion_min_dna_callers` DNA callers AND at least
  `rescue_promotion_min_rna_callers` RNA callers internally agree on Somatic
  (each modality's individual callers must be consistent), the variant is
  classified Somatic and tagged `RESCUE_PROMOTED=YES` / `RESCUED=YES` in INFO.
- Only Somatic agreement promotes; Germline/Reference agreement still returns
  NoConsensus.
- Disagreement between modalities' individual-caller classes still -> Artifact.
- Gated by `rescue_promotion_enabled` (default: ON).

### Artifact veto (audit M2)
When one modality's consensus label is Artifact and the other's is not, the
`rescue_veto_direction` config decides the outcome:

- `dna` (default): a DNA Artifact label outranks RNA non-Artifact evidence —
  the result stays Artifact; RNA can no longer flip a DNA-flagged artifact to
  Somatic. RNA Artifact vs DNA Somatic keeps the legacy behavior (the DNA
  non-Artifact label still wins when `cross_modality_min_callers_for_artifact`
  DNA callers support it, else Artifact).
- `rna`: mirrors the veto (RNA Artifact outranks DNA non-Artifact).
- `none`: legacy behavior — the non-Artifact side wins when it has at least
  `cross_modality_min_callers_for_artifact` callers, RNA checked first.

### Truthful rescue flags (audit M3)
`RESCUED`, `CROSS_MODALITY`, `PASSES_CONSENSUS_DNA`, and
`PASSES_CONSENSUS_RNA` are computed from records that PASSED as Somatic
(the record's FILTER/classification is Somatic), never from mere presence in
a union file — the union includes NoConsensus/Artifact/Germline records:

- `PASSES_CONSENSUS_DNA` / `PASSES_CONSENSUS_RNA`: YES iff the modality's
  consensus record is present AND its label is Somatic.
- `CROSS_MODALITY`: YES iff each modality contributed at least one record
  (consensus or individual caller) classified Somatic.
- `RESCUED`: YES iff the site passed as Somatic in both consensus sets, or
  was rescued by cross-modality promotion.
- `RESCUE_PROMOTED`: YES iff the site failed within-modality consensus and
  was promoted by agreeing DNA+RNA individual callers (see above).

## Configuration
Defaults live in `bin/vcf_utils/classification_config.py` (`DEFAULT_THRESHOLDS`)
and are exposed as `run_rescue_vcf.py` CLI flags:

| Config knob | CLI flag | Default |
|---|---|---|
| `rescue_promotion_enabled` | `--disable_rescue_promotion` | True |
| `rescue_promotion_min_dna_callers` | `--rescue_min_dna_callers` | 1 |
| `rescue_promotion_min_rna_callers` | `--rescue_min_rna_callers` | 1 |
| `rescue_veto_direction` | `--rescue_veto {dna,rna,none}` | `dna` |

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
  - PASSES_CONSENSUS_DNA / PASSES_CONSENSUS_RNA (passed-as-Somatic, see above)
  - RESCUED (passed-as-Somatic in both modalities, or promoted)
  - RESCUE_PROMOTED
  - CROSS_MODALITY (both modalities contributed a Somatic record)

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

## Evidence contract through annotation

Observed caller presence, eligible caller votes and Somatic agreement are separate. Downstream annotation must use allele-specific tumor AD when available and cannot promote an observed caller that failed the configured alternate-read floor. `DNA_VERIFICATION=rejected` vetoes annotation-based Somatic promotion; `confirmed`, `rejected` and `inconclusive` outcomes remain in the final rationale. Missing AD or normal evidence is unavailable, not zero.

RNA-only consensus remains a nomination for an opt-in DNA-verification policy. Realignment is a correlated reassessment of the same RNA reads and does not add an independent vote. The second rescue config records `ALIGNMENT_ROUND=realignment`; first-round defaults remain unchanged.
