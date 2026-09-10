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


## Native-evidence consensus experiment (2026-09-10)

An isolated, consensus-only policy was evaluated from cached caller VCFs; no
alignment or variant-calling process was rerun. The policy anchors labels on
DeepSomatic PASS calls, removes two explicit Mutect2 artifact combinations
(`contamination;germline;haplotype;panel_of_normals` and
`contamination;orientation;weak_evidence`), and admits SNV candidates absent
from the DeepSomatic set only when raw DeepSomatic QUAL is greater than zero
and Mutect2 TLOD is at least 12 with GERMQ at least 60. The rescue is
SNV-only because a WGS-IL indel candidate added an FP without a TP.

The som.py results, using the same truth, reference, HC regions, and UKB
target regions, were:

| Dataset | DeepSomatic TP/FP/F1 | Policy TP/FP/F1 |
| --- | ---: | ---: |
| WES-LL | 1048 / 38 / 0.6190 | 1051 / 36 / 0.6206 |
| WES-IL | 1365 / 21 / 0.7406 | 1371 / 21 / 0.7427 |
| WGS-IL | 2168 / 19 / 0.9663 | 2169 / 19 / 0.9666 |

The improvement is SNP-driven; indel performance is intentionally unchanged.
The thresholds were discovered on WES-LL and then held fixed for WES-IL and
WGS-IL, so this is evidence for an opt-in rule, not yet a default-policy
change. The next phase must compare indel-specific consensus and existing
first/realignment rescue VCFs from completed runs without touching workflow
outputs or caller caches.


## Cached rescue and indel analysis (2026-09-10)

This comparison used only completed WES-LL caller/rescue VCFs and existing som.py
scratch artifacts. No mapping, variant calling, or full workflow process was
rerun, and no caller cache or workflow output was modified.

| VCF | SNP TP/FP/FN | Indel TP/FP/FN | Record TP/FP/FN |
| --- | ---: | ---: | ---: |
| DNA consensus | 938 / 96 / 1267 | 35 / 1 / 60 | 973 / 97 / 1327 |
| DNA DeepSomatic | 1007 / 34 / 1198 | 41 / 4 / 54 | 1048 / 38 / 1252 |
| first rescue | 975 / 485 / 1230 | 36 / 30 / 59 | 1011 / 515 / 1289 |
| realignment rescue | 975 / 290 / 1230 | 36 / 20 / 59 | 1011 / 314 / 1289 |

The rescue branch does not currently outperform the DNA baseline: it adds 38
TPs but 418 FPs (first rescue) or 217 FPs (realignment rescue). Realignment
changes the error profile rather than recall: relative to first rescue it keeps
all 1011 TPs, removes 201 FPs, and introduces no net TP. The remaining
realignment-rescue set still has 290 SNP FPs and 20 indel FPs, far above the
DeepSomatic baseline.

Indels are the clearest negative result. DNA consensus is 35/1/60, already
more precise than DeepSomatic (41/4/54); first rescue becomes 36/30/59 and
realignment rescue 36/20/59. Thus indel rescue should remain disabled or use a
separate, substantially stricter policy. A broad high-confidence indel rescue
was tested separately and added an FP without a TP.

The FP attributes identify the rescue failure mode: among first-rescue FPs,
402/515 have zero DNA caller support and 416/515 have zero DNA callers
classified Somatic; among realignment-rescue FPs these counts are 199/314 and
213/314. Requiring at least one DNA Somatic caller still leaves 99--101 FPs,
while requiring three DNA Somatic callers reduces FPs to 9 but loses 85--121
TPs. This makes an unconditional RNA-only promotion unsuitable. The next
opt-in rescue experiment should preserve the DNA consensus and admit only
SNVs with explicit DNA evidence plus independent RNA corroboration/quality;
indels should require a separate evidence floor and be evaluated on their own.


## Indel candidate audit (2026-09-10)

The cached caller-level audit was extended to the held-out WES-IL and WGS-IL
inputs. The native-evidence policy was benchmarked without rerunning any
workflow stage:

| Dataset | DeepSomatic indels TP/FP/FN | Native policy indels TP/FP/FN |
| --- | ---: | ---: |
| WES-IL | 53 / 3 / 42 | 53 / 3 / 42 |
| WGS-IL | 86 / 12 / 9 | 85 / 11 / 10 |

Raw caller unions are not a viable indel strategy. On WES-IL, Mutect2 and
Strelka expose thousands of raw indel candidates but contribute only one
additional truth indel beyond DeepSomatic while adding thousands of FPs. Even
the filtered consensus adds only one TP and one FP. WGS-IL shows the same
pattern at smaller scale: one added TP is paired with one FP. The current
policy therefore leaves indels anchored to the strongest baseline caller; any
future indel rescue must use an independently validated, indel-specific
evidence model rather than the SNP rescue thresholds.
