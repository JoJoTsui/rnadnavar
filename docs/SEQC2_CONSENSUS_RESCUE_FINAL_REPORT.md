# SEQC2 Consensus and Rescue Final Report

## Scope

This report freezes the current consensus/rescue interpretation for the SEQC2
WES-LL hybrid benchmark. Caller VCFs, alignment outputs, and workflow caches
remain read-only for this report. The validated comparisons used the same truth
VCF, high-confidence BED, UKB target BED, reference genome, and `som.py -N`
benchmark contract.

## Decision diagrams

### Consensus (within one modality)

```text
caller VCF panel
      |
      v
normalize CHROM:POS:REF:ALT + classify each caller
      |
      v
eligible vote?  -- no --> Artifact/NoConsensus evidence
      | yes
      v
alt-read floor + non-Artifact support count
      |
      +--> SNV and native policy OFF: support >= snv_thr?
      |
      +--> SNV and native policy ON: qualified DeepSomatic
      |                         OR TLOD>=12 & GERMQ>=60
      |                         with Mutect2 artifact vetoes
      |
      +--> indel: support >= indel_thr (native override never applies)
      |
      v
majority biological class
      |
      +--> unique winner: Somatic/Germline/Reference
      +--> tie: Artifact
      +--> insufficient support: NoConsensus

## Consensus rules

### Default consensus

- Aggregate the explicit caller panel by normalized `CHROM:POS:REF:ALT`.
- Count only eligible non-Artifact caller records that satisfy the configured
  tumor alternate-read floor (`min_alt_support`, default 3 when AD exists).
- Require `snv_thr` callers for SNVs and `indel_thr` callers for indels (both
  default 2).
- Use the unified biological FILTER vocabulary: `Somatic`, `Germline`,
  `Reference`, `Artifact`, and `NoConsensus`.
- A top-class tie is `Artifact`; it is never forced to Somatic.
- Preserve per-caller FILTER, support, genotype, depth, VAF, and rationale INFO.

### Native-evidence SNV policy

The implemented policy is opt-in via `--native-evidence-snv` or
`params.native_evidence_snv=true`; the default remains disabled to preserve
existing workflow behavior and caller-cache reuse.

When enabled for SNVs:

1. Retain a qualified DeepSomatic Somatic record.
2. Admit a non-DeepSomatic candidate only when DeepSomatic QUAL is positive and
   Mutect2 TLOD is at least 12 and GERMQ is at least 60.
3. Reject the candidate when Mutect2 carries either explicit contamination /
   germline / haplotype / panel-of-normals or contamination / orientation /
   weak-evidence artifact combination.
4. Never apply this override to indels.

The policy emits `CLASSIFICATION_RATIONALE=rule:native_evidence_snv` and leaves
all original caller evidence intact.

## Rescue rules

### Rescue (cross-modality)

```text
DNA consensus + RNA consensus + individual caller evidence
                         |
                         v
             ignore NoConsensus as positive evidence
                         |
                         v
     both modality labels present?
       | yes                                | no
       v                                    v
 same label -> keep label             one label -> apply it
 conflicting labels -> support/veto     no labels -> promotion test
       |                                    |
       v                                    v
 DNA Artifact veto (default)             DNA Somatic callers >= threshold
       |                                    AND RNA Somatic callers >= threshold
       v                                    |
 Artifact / supported non-Artifact        v
                                         RESCUE_PROMOTED=YES

validated opt-in rescue gate for rescue-only SNVs:
N_DNA_CALLERS_SUPPORT >= 1 AND N_RNA_CALLERS_SOMATIC >= 2
(indels are not promoted by this gate)

### Current workflow rescue contract

- Inputs are DNA consensus, RNA consensus, and individual caller VCFs.
- `NoConsensus` is absence of modality consensus, not positive evidence.
- DNA artifact veto is the default when DNA and RNA disagree.
- Cross-modality promotion requires internally consistent Somatic evidence from
  the configured DNA and RNA caller thresholds.
- `RESCUED`, `RESCUE_PROMOTED`, `CROSS_MODALITY`, and
  `PASSES_CONSENSUS_*` are based on Somatic-labelled evidence only.
- Realignment is a correlated reassessment of RNA reads, not an independent
  caller vote.

### Validated gated-rescue policy

For rescue-only SNVs, the validated WES-LL gate requires:

- `N_DNA_CALLERS_SUPPORT >= 1`; and
- `N_RNA_CALLERS_SOMATIC >= 2`.

Existing native-consensus records are always retained. Indels are not rescued by
this gate.

This gate is validated as a candidate policy; the current production rescue
classifier still uses its documented configurable promotion contract. Enabling
the gate as the default requires a code-path change and cross-cohort validation.

## WES-LL benchmark result

| Output | SNP TP/FP/FN | Indel TP/FP/FN | Overall TP/FP/FN | P / R / F1 |
| --- | ---: | ---: | ---: | ---: |
| DNA Mutect2 | 906 / 147 / 1299 | 32 / 4 / 63 | 938 / 164 / 1362 | 0.8512 / 0.4078 / 0.5514 |
| DNA Strelka2 | 944 / 1524 / 1261 | 38 / 41 / 57 | 982 / 1565 / 1318 | 0.3856 / 0.4270 / 0.4052 |
| DNA consensus | 938 / 96 / 1267 | 35 / 1 / 60 | 973 / 97 / 1327 | 0.9093 / 0.4230 / 0.5774 |
| DNA DeepSomatic | 1007 / 34 / 1198 | 41 / 4 / 54 | 1048 / 38 / 1252 | 0.9650 / 0.4557 / 0.6190 |
| ClairS | 237 / 0 / 1968 | 0 / 0 / 95 | 237 / 0 / 2063 | 1.0000 / 0.1030 / 0.1868 |
| Native consensus + gated rescue | 1021 / 36 / 1184 | 41 / 4 / 54 | 1062 / 40 / 1238 | 0.9637 / 0.4617 / 0.6243 |

The validated policy improves TP, recall, and F1 over DNA DeepSomatic, ties it
on indels, and adds two overall FPs. Precision is therefore slightly lower and
must remain an explicit acceptance criterion in future rescue tuning.

## Reproducibility and cache safety

- The native consensus flag is opt-in and changes only the consensus process
  when enabled.
- Caller alignment and variant-calling processes are not modified by this
  policy.
- The WES-LL gated-rescue comparison was generated from completed VCFs; no full
  workflow rerun is required.
- WES-IL and WGS-IL currently lack completed rescue VCF artifacts, so their
  rescue validation is intentionally deferred.

Implementation references: [consensus rules](consensus_vcf_rules.md), [rescue
rules](rescue_vcf_rules.md), and [optimization follow-up](review/SEQC2_OPTIMIZATION_FOLLOWUP.md).
