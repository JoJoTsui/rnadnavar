# Updated consensus/rescue validation — interview record

Date: 2026-09-14. Starting commit: `76c9f47e`.
Status: decisions Q1–Q9 accepted; awaiting final confirmation to execute. This is not a
completed validation report or implementation authorization.

## Accepted boundaries

- Q1: audit actual defaults and integration gaps, but do not activate or change
  global defaults pending independent validation. Names of enabled flags do
  not prove equivalence to the standalone experimental policy.
- Q2: existing completed WES-LL/WGS-IL outputs only; targeted alignment reads
  and isolated consensus/rescue/annotation runs may use fresh output locations.
  No mapping, calling, cohort run, HG008 interference, original-output edits,
  or cache changes. Verify all ingress modes with configuration/focused tests.
- Q3: reproduce frozen candidate SNP/indel/overall metrics in both regions;
  integration must not lose TP or add FP. Check alleles and the complete
  training-label contract, not only metrics. Report DeepSomatic separately;
  ties are acceptable. Historical indel copying remains provisional. Audit
  missing/conflicting evidence, multiallelic and duplicate records, sample
  identity, normal contamination, and RNA evidence dependence.

- Q4: freeze experimental alleles, thresholds, candidate universe, and
  provisional indels. Generalized implementations are separate candidates,
  not replacements for the reference experiment.
- Q5: retain the frozen baseline for parity; independently flag common-AF,
  matched-normal, or other biological contradictions as potential training-label
  defects. Removal is a separately evaluated alternative, not a silent edit.
- Q6: inspect rescued TPs, remaining FPs, excluded sites, and matched controls
  before considering new thresholds. Missing/low coverage and ambiguity stay
  explicit. RNA realignment is not independent biological confirmation.
- Q7: do not extrapolate DNA-calibrated native rules automatically to RNA or
  RNA realignment. Validate those branches and both rescue rounds separately;
  consistency across input formats does not require identical modality policies.

- Q8: validation may create reproductions and regression tests and correct
  validation scripts. Production fixes, policy semantics, and default changes
  are follow-on work requiring implementation approval.
- Q9: wrong sample roles, invented evidence, unexplained allele/FILTER changes,
  annotation loss, and incorrect final-rescue lineage block approval regardless
  of F1. Insufficient read coverage is inconclusive, not passed. Independent
  validation and the provisional indel policy remain explicitly unresolved.
Architectural boundary: [ADR-0008](../adr/0008-validate-frozen-policy-before-default-parity-claims.md).

## Read-only defaults audit

Current defaults are **not equivalent** to the successful experimental rules.

| Layer | Evidence | Consequence |
| --- | --- | --- |
| Root defaults | `nextflow.config:109` enables native SNV and rescue promotion, DNA/RNA floors 1/2 | Older flags are on, not the newly tested gate |
| Forwarding | `conf/modules/consensus/vcf_consensus_workflow.config:12`; consensus/rescue local modules | Generic flags propagate across input modes |
| Native classification | `bin/vcf_utils/classification.py:483` and `:518`; `aggregation.py:1126` | Native override is additive; ordinary majority remains a fallback |
| Historical replay | `examples/seqc2/scripts/replay_historical_native_gate.py:46` | M2 combination veto precedes baseline retention; no ordinary-majority fallback |
| Indels | `bin/vcf_utils/aggregation.py:1094` versus historical replay | Production threshold consensus differs from historical DS-derived retention |
| Updated rescue gate | `examples/seqc2/scripts/test_rescue_nomination_gate.py:36` and `:47` | Nomination/population/editing predicates are experimental only |
| Baseline retention | `bin/vcf_utils/variant_classifier_unified.py:309` versus experimental `native | kept` | Production may reclassify DNA Somatic during modality disagreement |
| Vote semantics | `variant_classifier_unified.py:384`–`:392` | Comment says DNA Somatic; condition counts eligible lists, which can include Reference |
| Annotation order | `subworkflows/local/vcf_rescue_post_processing/main.nf:52` and `:79` | Biological annotations arrive after rescue; final annotated experimental gate is not raw rescue classification |

The generic workflow options apply to hybrid and FASTQ-triplet consensus;
DNA-only mode has no cross-modality rescue without RNA. First and second rescue
share the rescue workflow. Standalone consensus Python defaults also differ
from Nextflow defaults, so direct CLI and workflow checks must be separate.

These are code-inspection findings, not evidence of completed execution parity.
They do not authorize editing production policy before the interview concludes.

## Execution handoff

The decision frontier is empty. Await final user confirmation of the shared
understanding before starting validation execution. Then verify frozen inputs,
reproduce the benchmark matrix, inspect targeted read evidence, exercise
consensus/rescue and annotation contracts, and report per-check pass/fail/
inconclusive outcomes with reproducible commands and follow-on defect records.
Do not claim overall validation complete while required evidence is missing.

No new workflow execution or production code changes were performed in this
interview round. The factual audit used a read-only sub-agent as directed by
the grilling skill; the main agent owns the decisions and documentation.
