# SEQC2 optimization follow-up

Review date: 2026-09-08. Inspected branch: `seqc2-consolidated`, HEAD `152d304`.
This document records accepted direction and investigation findings, not an implemented redesign.

## Accepted direction

The user accepted the recommendations in the new interview round Q1–Q3:

- Evaluate DNA consensus, first rescue, and realignment rescue separately against DNA DeepSomatic: higher overall F1 with precision at least as high, on the same benchmark domain. Report SNP and INDEL separately, with uncertainty; do not conceal an INDEL regression behind SNP performance.
- Use WES-LL for discovery, reserve genomic blocks before tuning, and validate frozen rules on WES-IL/WGS-IL. These are same-cell-line robustness checks, not independent biological validation. Preserve chr1 as the previously agreed holdout; prior exploratory inspection means it is not a wholly unseen dataset.
- Establish artifact provenance before regenerating a current-code baseline. Attribute errors before optimizing rules. Preserve original inputs/outputs, biological FILTER vocabulary, and the seq2neo DN/DT/RT FASTQ triplet workflow.

The latest Q1 acceptance sets the overall acceptance gate for this follow-up. The earlier review's stricter simultaneous precision-and-recall improvement per variant type remains a separately reported superiority target, not a result already demonstrated.

Prior accepted design choices remain applicable: DeepSomatic-preserving integration as an experiment; RNA nominations need DNA tumor/normal verification; inconclusive evidence is not Somatic; realignment is not an independent vote; freeze numerical rules after development analysis. See [previous review](SEQC2_COMPREHENSIVE_REVIEW.md) and [ADR 0005](../adr/0005-verify-rna-nominations-before-somatic-labeling.md).

## Evidence gathered

The comprehensive comparison is under `examples/seqc2/hybrid/comparison/comprehensive_realign/WES_LL_T_1_vs_WES_LL_N_1/`.

| Callset | Records TP | FP | FN | Precision | Recall | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| DNA DeepSomatic | 1048 | 38 | 1252 | .9650 | .4557 | .6190 |
| DNA consensus | 985 | 107 | 1315 | .9020 | .4283 | .5808 |
| First rescue | 1002 | 539 | 1298 | .6502 | .4357 | .5217 |
| Realignment rescue | 1002 | 508 | 1298 | .6636 | .4357 | .5260 |

These are historical results. Net count differences are not paired allele transitions.

### Provenance and metric limitations

- The final execution trace (`hybrid/output/seqc2.wes.ll.hybrid.realign.full/pipeline_info/execution_trace_2026-09-08_09-14-59.txt`) reports DNA consensus `2d/a02415` cached from September 7 08:32 and initial RNA consensus `0a/1fc393` cached from 08:30. Both predate the 10:37 commit `076bb42` connecting eligible support to classification. Realigned consensus completed September 8. Commit times alone do not identify the exact staged script contents; source snapshots/checksums must establish that binding. The comparison cannot be assumed to evaluate a uniform current-code baseline.
- The comprehensive directory retains metrics and PASS query copies, but no per-site TP/FP/FN VCFs or feature tables. Native precision/recall intervals are omitted by the aggregate CSV.
- Records are not necessarily SNP plus INDEL: DNA Mutect2 FP is 164 versus 147 SNP plus 4 INDEL. Account for other/complex records explicitly rather than forcing equality.
- Installed som.py uses normalization plus allele intersection. Representation-sensitive indel/complex disagreements need a separate haplotype-aware check.
- Its reported FP region size is 2,875,001,522 bases, inconsistent with the HC/target intersection. This affects FP-per-region reporting, not evidence that the region-restricted precision/recall counts are wrong.
- Index-existence cache checks do not bind benchmark derivatives to source content; version/checksum provenance is incomplete. The comprehensive orchestrator is a temporary script rather than a tracked reusable entry point.

### MAPQ evidence

A bounded inspection at `chr1:1000000-1100000` in executed task `hybrid/work/b3/04de60d5111151e63a3fa3dd40ac46` found:

| Stage | MAPQ 255 | MAPQ 60 |
| --- | ---: | ---: |
| MarkDuplicates CRAM | 26917 | 0 |
| SplitNCigarReads CRAM | 0 | 29410 |

The command does not disable the default MAPQ transformation. Splitting/filtering changes record counts, so this is not a read-for-read conversion count. It establishes conversion in this region and path only. Inspect final caller-ready alignments and bypass paths before choosing an intervention. GATK documents the default [255-to-60 transformation](https://gatk.broadinstitute.org/hc/en-us/articles/30332077175963-SplitNCigarReads).

### Current method concerns requiring attribution

- Equal categorical voting can reject DeepSomatic singleton Somatic calls and admit agreement from less precise callers. Historical paired evidence exists, but must be refreshed for this output and current code.
- Strelka normal-reference status can become a Reference class although a reference normal is compatible with a tumor somatic event. Detection support and somatic evidence must remain distinct.
- Missing alternate-count evidence can bypass the support floor.
- RNA-only consensus can pass through when DNA is absent or NoConsensus without verifying DNA callability and normal evidence.
- DNA Artifact vetoes and internal-unanimity promotion may suppress recoverable variants; relaxing them without evidence risks label precision.
- Rescue filtering records QC in `RaVeX_FILTER` while preserving biological FILTER. A Somatic-only benchmark and a QC-selected benchmark evaluate different products.
- Prior standalone verification/backbone/read-recovery helpers are not integrated into the workflow. Inspection found deletion counting, candidate TSV parsing, evidence provenance, and VCF header compatibility gaps; do not treat those helpers as validated production implementations.

## Next evidence phase

1. Freeze hashes, executed selectors, tool identities, task sources, and benchmark domains; regenerate consensus/rescue from frozen callers in an isolated location under one code version.
2. Retain normalized scored alleles and pair transitions to original evidence. Report representation disagreements and other/complex records separately.
3. Separate truth absent from candidate inputs, caller-rejected truth, integration-rejected truth, and post-filter losses. For FP, distinguish inherited DNA errors, consensus additions, RNA-only pass-through, explicit promotion, and annotation/QC transitions.
4. Join tumor/normal depth, alternate support, VAF, mapping/base quality, strand/read-position, editing, repeat, and splice context. Report missing measurements as unknown.
5. Compare frozen single-factor ablations against the refreshed baseline, then evaluate held-out regions. Do not tune using held-out results.

No workflow, classifier, benchmark implementation, or original artifact has been changed by this follow-up. Product selector and representation-adjudication decisions are the next interview frontier.
