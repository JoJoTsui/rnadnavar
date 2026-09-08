# SEQC2 review: consensus, RNA tumor, rescue, and realignment

Date: 2026-09-07. Review branch: `review/seqc2-consensus-rt-rescue`; worktree: `/tmp/rnadnavar-seqc2-review`; code baseline: `931c476a081162e6398130a5fd28b1fd8aced11c`.

## Conclusion

There are two measured problems, with different remedies. DNA consensus sometimes purchases precision with lost DeepSomatic true positives, but in WES-LL its equal categorical vote also adds many errors. First-round hybrid rescue then admits a largely false-positive RNA-only set. Its main failure is **RNA-only consensus pass-through**, not the explicit 1-DNA + 1-RNA promotion rule. Separately, a confirmed support/FILTER inconsistency undermines alt-read threshold enforcement, and process argument wiring prevents ordinary parameter ablation.

The user clarified that downstream training currently consumes only the **second-round rescue output using realigned RNA**, based on their observation that it contains fewer false positives. The first-round SEQC2 findings below diagnose intermediate behavior; they are not measurements of the training artifact's precision. No completed SEQC2 second-round rescue benchmark exists in the inspected evidence. The reported lower false-positive burden in existing realignment results motivates retaining that endpoint, while its magnitude and recall trade-off require a matched comparison.

This is a completed evidence review and a proposed experiment plan, not an implemented classifier redesign. Original inputs, results, and local checkout changes were preserved. The glossary was extended with benchmark-domain and paired-transition terms. The user accepted per-type simultaneous precision/recall improvement as the superiority criterion and DNA verification before RNA-only nominations become Somatic labels. The latter boundary is recorded in [ADR 0005](../adr/0005-verify-rna-nominations-before-somatic-labeling.md). All seven review-level policy decisions are accepted; numerical rules and implementation remain future work. The second-round rescue artifact is the training endpoint.

## 1. Current results and their actual domain

The main comparisons use SEQC2 v1.2.1 high-confidence truth, high-confidence regions, and the padded UKB/Broad target union. All below use 2,300 truth records (2,205 SNVs; 95 indels). **WGS here means WGS input evaluated inside those target regions.** The hybrid rescue row is first-round only and is not the user's downstream training endpoint. Percentages are from native som.py metrics, not from the older files named `corrected`.

| Input / final query | TP | FP | FN | Precision | Recall |
| --- | ---: | ---: | ---: | ---: | ---: |
| WES-LL DeepSomatic DNA | 1,048 | 38 | 1,252 | 96.50% | 45.57% |
| WES-LL DNA consensus | 985 | 107 | 1,315 | 90.20% | 42.83% |
| WES-IL DeepSomatic DNA | 1,365 | 21 | 935 | 98.48% | 59.35% |
| WES-IL DNA consensus | 1,269 | 17 | 1,031 | 98.68% | 55.17% |
| WGS-IL DeepSomatic DNA | 2,168 | 19 | 132 | 99.13% | 94.26% |
| WGS-IL DNA consensus | 2,131 | 5 | 169 | 99.77% | 92.65% |
| WES-LL hybrid RT DeepSomatic | 305 | 679 | 1,995 | 31.00% | 13.26% |
| WES-LL hybrid RT consensus | 315 | 611 | 1,985 | 34.02% | 13.70% |
| WES-LL hybrid filtered first rescue | 1,002 | 715 | 1,298 | 58.36% | 43.57% |

The successful hybrid baseline is `examples/seqc2/hybrid/output/seqc2.wes.ll.hybrid.pooling_fix`, traced on September 6. Its DNA metrics match the earlier WES-LL baseline. Cached DNA/RT benchmark queries match their production consensus records, and the saved rescue query exactly matches the transformed **filtered** first-round rescue output, not raw rescue. [Content verification](seqc2-evidence/query_cache_verification.json).

Per-type results show why a pooled score is insufficient:

| Input | DeepSomatic SNV P / R | Consensus SNV P / R | DeepSomatic indel P / R | Consensus indel P / R |
| --- | --- | --- | --- | --- |
| WES-LL DNA | 96.73 / 45.67% | 89.96 / 43.08% | 91.11 / 43.16% | 97.22 / 36.84% |
| WES-IL DNA | 98.65 / 59.50% | 98.94 / 55.28% | 94.64 / 55.79% | 92.59 / 52.63% |
| WGS-IL DNA | 99.62 / 94.47% | 99.95 / 93.33% | 88.54 / 89.47% | 94.81 / 76.84% |

First rescue adds **zero indel TP and 27 indel FP** to DNA consensus (35 TP remain); SNVs gain 17 TP and 581 FP. Thus a gain in combined recall says nothing favorable about indel rescue in this run.

For context, WES-LL Mutect2 has 938 TP/164 FP and Strelka 982 TP/1,565 FP. RT Mutect2 has 215 TP/595 FP and Strelka 375 TP/5,287 FP. The panel members have very different error rates here; treating their categorical votes as interchangeable discards this information. [All main per-type metrics](seqc2-evidence/main_metrics.csv); [native metric paths and executed scoring commands](seqc2-evidence/benchmark_inventory.json).

## 2. What actually changes the errors

The following paired attribution uses exact CHROM/POS/REF/ALT sets and point membership in HC∩target. It is a diagnostic alongside normalized som.py: a few indel boundary/representation differences change totals. “Nonmatch” therefore means no exact truth key, not independent biological validation as false. The rescue increments agree with the official aggregate increments.

| Change | Added truth matches | Added nonmatches | Lost truth matches | Removed nonmatches |
| --- | ---: | ---: | ---: | ---: |
| DeepSomatic DNA → WES-LL DNA consensus | 8 | 80 | 71 | 9 |
| DNA consensus → filtered first rescue | 17 | 608 | 0 | 0 |
| RNA-only consensus rule within rescue additions | 17 | 606 | — | — |
| Explicit cross-modality promotion within rescue additions | 0 | 2 | — | — |

All 80 consensus-added nonmatches have `Somaticx2+Referencex1` majority rationale. The panel can therefore reject useful DeepSomatic calls while admitting two-caller agreements that DeepSomatic did not pass. The same mechanism need not fail identically in WGS or another WES replicate; the measured WGS precision gain confirms that sample and variant context matter.

Of the 625 rescue additions, 623 use `rna_consensus_only`. Their incremental exact-match yield is 17/623 = **2.73%**; all additions yield 17/625 = **2.72%**. Tightening 1+1 promotion would address only two of these added nonmatches. Even three RNA Somatic votes are not sufficient protection: the RNA-only addition stratum with no DNA votes and three RNA Somatic votes contains two truth matches and 90 nonmatches.

The default DNA Artifact veto cannot protect sites with absent DNA consensus or DNA NoConsensus. Those are treated as no substantive DNA label, allowing the RNA label through. This distinction is crucial: “no DNA call” is not proof of adequate DNA coverage, normal-reference evidence, or a true RNA-specific event. [Reproducible paired audit](seqc2-evidence/exact_allele_audit.json); [audit script](seqc2-evidence/reproduce_exact_allele_audit.py).

## 3. Consensus and rescue defects versus policy trade-offs

**Confirmed correctness defect: the alt-read floor is disconnected from final FILTER.** Aggregation calculates eligible support using the alt-read floor, but the classifier counts all normalized labels. Real full VCFs contain Somatic/PASSES_CONSENSUS=NO records: WES-LL 27/2,511 Somatic; WES-IL 9/3,239; WGS 0/6,033; RT 43/1,293. Some have ENS_SUPPORT=0/3. Synthetic direct calls reproduce both consensus and promotion accepting two-alt-read evidence with zero eligible support.

This is not a complete explanation of the accuracy loss. Inside the diagnostic domain, 12 of 22 inconsistent DNA records and eight of 35 inconsistent RT records match truth. Enforcing the declared floor may improve consistency while sacrificing some sensitivity. Measure that trade-off rather than claiming a free precision-and-recall gain.

**Confirmed control defect: `ext.args` is assigned but omitted from both module commands.** Normal Nextflow overrides for `--min_alt_support`, promotion thresholds, promotion disablement, and veto direction do not reach the scripts. Consensus config also sets `ext.thr` while its module reads `ext.snv_thr`/`ext.indel_thr`. The effective defaults are 2/2 observed-record thresholds, alt floor 3, promotion enabled at 1/1, DNA veto. Raising the threshold to three currently means three observed caller records followed by majority class; it does not mean three Somatic votes.

**Policy limitations:** singleton DeepSomatic calls fail the two-record threshold; categorical ties become Artifact; two weaker calls can outvote the stronger caller. A Strelka filtered call with normal NT=ref can be normalized to Reference, although reference normal is compatible with a tumor somatic event. Reference/Germline labels can still count as non-Artifact detection support. These categories are not interchangeable measures of caller rejection or somatic confidence.

**Rescue restrictions also lose opportunities:** the DNA veto can lock in an Artifact that originated from a tie; promotion requires internal class unanimity, so extra contradictory records block a 1+1 opportunity. These are intentional policies to ablate, not reasons to disable protections wholesale.

**Filtering semantics:** `filter_rescue_vcf.py` preserves biological FILTER and records, attaching QC reasons in `RaVeX_FILTER`. Selecting only FILTER=Somatic does not apply those INFO exclusions. COSMIC/gnomAD and editing annotations can change labels separately. Score each stage and selector explicitly. [Full code references, reproductions, effective controls, and tests](seqc2-logic-evidence.md).

## 4. RT tools and versions

The completed RT baseline used STAR **2.7.11b**, GATK **4.6.1.0**, Strelka **2.9.10**, reported DeepSomatic **1.9.0**, Nextflow **25.10.0**, vt **2015.11.10**, and bcftools norm **1.22**. Other annotation/stats tasks report bcftools **1.21**. DeepSomatic's reported version is written as a constant by the task: verify executable/image identity before relying on it as an independent version check.

Material executed settings:

- RT DeepSomatic uses `--model_type=WES` against DNA normal. The documented model interface does not list an RNA-seq model. This is an unvalidated model/input distribution assumption; the observed RT precision of 31% is direct evidence that its DNA performance does not transfer here. [Official model interface](https://github.com/google/deepsomatic).
- RT Strelka uses `--exome`, with no Manta candidate indels. Manta candidates are recommended upstream for somatic indel performance, but re-enabling them needs its own controlled experiment; do not invent a somatic `--rna` fix. [Strelka guidance](https://github.com/Illumina/strelka).
- STAR uses two-pass alignment, sjdbOverhang 75, up to 20 multimaps, match/score fractions 0.33, and splice overhang minimum 1. RNA then undergoes MarkDuplicates, SplitNCigarReads and BQSR. Mutect2 excludes soft-clipped bases, uses PoN/gnomAD plus orientation/contamination filtering, and FilterMutectCalls sets `--max-events-in-region 5`.
- **REDIportal annotation was disabled** in the successful hybrid. `rna_filtering` in the tools list does not activate that annotation. Editing is a plausible contributor, but its fraction among errors has not been measured here.

The three RT callers share the same alignment and RNA biology; their agreement cannot be treated as independent confirmation. Do not upgrade every tool or change every threshold together. First isolate RNA admission, editing annotation, mapping/context effects, and model assumptions. [Executed command paths, version manifest, and primary-source detail](seqc2-rt-tools-evidence.md).

## 5. Realignment: execution blocker and limits

The latest inspected full realignment attempt fails with a Nextflow `MissingMethodException` on an ArrayList passed to a map closure (September 7 log). The source's remainder-emitting join followed by a fixed-arity closure is consistent with this failure; the exact unmatched tuple still needs a focused reproduction. This is distinct from the earlier STAR exit 104, whose later retry succeeded. No measured HISAT2 benefit can be inferred.

The current design also limits possible recovery:

- Candidates are all first-pass RT-consensus records, including negative classes; DNA-only sites absent from RT are not systematically nominated. Candidate BEDs have no padding.
- Read names are selected at candidates, then recovered by name, so mates elsewhere in the supplied alignment can be retained. However, the supplied RT alignment is already processed by SplitNCigarReads/BQSR. Recovery cannot restore hard-clipped bases; singleton/other FASTQ outputs are not forwarded to the paired remapping route.
- HISAT2 is declared as 2.2.1, with `--no-mixed --no-discordant`; this is source configuration, not a successful measured second-pass version. Read-subset processing and strict pairing may remove artifacts and true support alike.
- Second rescue replaces first-pass RNA evidence with second-pass evidence; it does not inherently demand independent first/second-pass agreement. The same RNA-only admission policy would remain unless deliberately changed.

Fixing execution is necessary to evaluate realignment, but does not solve RNA-only label admission. Annotation is enabled in the realignment configuration, so compare both rounds under matched annotation policy to isolate alignment effects. [Detailed second-pass trace and source constraints](seqc2-rt-tools-evidence.md).

The local training manifest audit supports the user's endpoint clarification: all 66 original rescue entries use realigned RNA, and the rerun manifest retains that lineage. Its 63 eligible training entries comprise 56 rerun rescue VCFs and seven QC-cleaned derivatives; three samples are excluded. A rerun output directory named simply `rescue` does not mean first-round RNA was used. See [training-label lineage evidence and helper-policy caveats](seqc2-training-lineage.md).

## 6. Benchmark interpretation limits

The low absolute WES recall includes capture/coverage and target-domain effects. Historical `no_target`/`corrected` scores use 41,072 truth records, while `capture_bed` replaces HC restriction and counts many out-of-truth-domain calls as FP. Some historical alternate consensus comparisons contain zero query records. They are not comparable evidence of the current workflow's performance.

DeepSomatic r1.9 HCC1395 examples evaluate **chr1 held out during training**. Whole-target HCC1395 is useful for regression but not a clean held-out superiority claim. Our lightweight chr1 diagnostic already changes the ranking of raw TP counts: DNA consensus has 127 truth matches versus DeepSomatic's 126, while precision remains lower (138 versus 130 calls). Do not tune on this holdout or use these point-key diagnostics as a substitute for normalized official scoring. [Version-specific WGS example](https://raw.githubusercontent.com/google/deepsomatic/r1.9/docs/deepsomatic-case-study-wgs.md).

The evaluator is installed hap.py/som.py **0.3.15**, performing normalization and allele intersection. Complex-event sensitivity should also be checked with haplotype-aware comparison. The script's index-existence cache can become stale in future runs, although the inspected hybrid queries were content-verified. [Methodology audit](seqc2-benchmark-methodology.md).

## 7. Prioritized experiment plan

| Order | Question / bounded experiment | Evidence required before accepting a change |
| --- | --- | --- |
| 1 | Emit actual process arguments; connect eligible support to class decisions | Counterexample regression, missing-AD policy, exact emitted command, per-rule TP/FP changes |
| 2 | Compare baseline first rescue with RNA-only additions excluded; promotion unchanged | Paired transition table showing which of the 17 gained TP and 606 RNA-only nonmatches are removed |
| 3 | Apply editing annotation alone to fixed baseline artifacts | REDIportal overlap and label transitions; no simultaneous realignment/configuration changes |
| 4 | Compare current DNA vote with a DeepSomatic backbone plus separately qualified additions | Retention of DeepSomatic TP, incremental precision of additions, selected FP removals, per-type results |
| 5 | Verify RNA-nominated sites using original DNA tumor/normal evidence | Coverage, alt reads, strand/read-position/mapping context; distinguish unavailable evidence from convincing absence |
| 6 | Resolve realignment tuple failure and validate recovered reads on a bounded case | Correct candidate identity, pair/singleton/read-length accounting, successful exact second-pass artifacts |
| 7 | Compare realigned versus initial RT with identical annotation and domain | Added/lost TP and FP, candidate-conditional and end-to-end recall, splice/indel strata |
| 8 | Tune tool options or replace models one factor at a time | Fixed versions/checkpoints, held-out evaluation, stable improvement across relevant input/variant strata |

A DeepSomatic backbone is an experiment, not a guarantee of dominance. If all baseline calls are retained, added calls must have precision above baseline to improve precision while increasing recall. Otherwise some baseline FP must be selectively removed without losing the recall gains. Caller voting alone does not supply that discrimination. The present RNA-only incremental yield is far below any such requirement.

## 8. Interview frontier and documentation status

### Round 1 — accepted

The user answered “all use recommendations,” accepting both recommendations presented at the end of the review.

- **Q1 — Success criterion:** precision and recall must both improve over DeepSomatic separately for SNVs and indels on the same evaluation domain to claim superiority. Precision-first alternatives may be reported separately for training labels; F1 does not replace this criterion.
- **Q2 — RNA-only admission:** RNA-only evidence nominates sites for DNA tumor/normal read-level verification. Unresolved nominations remain outside the Somatic label set. This is a desired policy, not the behavior of the reviewed code; see [ADR 0005](../adr/0005-verify-rna-nominations-before-somatic-labeling.md).

### Round 2 — accepted

The user again answered “all use recommendations,” accepting Q3–Q5.

- **Q3 — DNA integration:** use DeepSomatic DNA PASS as the starting callset for the redesign experiment. Add or remove calls only through explicit, auditable evidence rules; another caller's absence or categorical majority alone does not overturn DeepSomatic. Keep current consensus as the comparison control. This experimental policy does not guarantee improvement.
- **Q4 — DNA verification:** direct evidence from original caller-ready DNA tumor/normal alignments can verify an RNA nomination even without a DNA caller PASS. Insufficient coverage or ambiguous evidence remains inconclusive and outside Somatic. Numerical gates remain to be specified.
- **Q5 — Evaluation scope:** assess WES-LL, WES-IL and WGS-IL separately with matched per-pair domains and SNV/indel precision/recall gates; gains in one pair do not compensate for regressions elsewhere. Preserve historical domains for regression, add declared capture/coverage strata, keep chr1 out of tuning, and require independent validation for generalization claims. Report uncertainty for small indel strata.

### Design tree

```text
Accepted Q1: precision AND recall improvement, separately SNV/indel
  ├─ Accepted Q3: DeepSomatic DNA starting set + evidence-based changes
  └─ Accepted Q5: separate input-pair gates; held-out validation
       └─ Accepted Q7: threshold development and freezing strategy
Accepted Q2: RNA nomination requires DNA verification
  └─ Accepted Q4: direct DNA reads can verify; inconclusive is not Somatic
       ├─ Accepted Q6: realignment evidence and read-source policy
       └─ Accepted Q7: verification threshold development
```

### Round 3 — accepted

**Q6 — Realignment role (accepted):** A successful second RNA alignment is not an additional independent vote for Somatic. Use realignment to reassess RNA nominations and identify alignment-sensitive evidence; DNA verification remains required. For the revised route, recover selected read pairs from original RT reads, preserving mate and library identity, rather than relying on sequence recovered from processed CRAM. Preserve the existing candidate policy initially so candidate expansion can be measured separately, and compare initial/realigned RNA under identical annotation settings. These are accepted design changes; the existing realignment route has not been fixed or rerun.

**Q7 — Evidence thresholds (accepted):** Do not prescribe universal numeric cutoffs before the development analysis. Use a bounded, declared development analysis to select verification and override rules from tumor/normal depth, alt support, mapping/base quality, and artifact context. Freeze rules before evaluating chr1 or an independent validation sample. Explicitly define missing-evidence handling and the relationship between eligible caller support and final FILTER; demonstrate both with regression examples. Report precision/recall trade-offs and reject a superiority claim if the accepted gates fail. Current alt-floor 3 remains a baseline to test, not an accepted optimum or an automatic override of every DeepSomatic PASS call.

The user accepted Q6 and Q7 and clarified that only the realignment/second-round rescue output is used for downstream training because they have observed fewer false positives there. This establishes the final training endpoint; it does not turn the first-round SEQC2 metrics into an evaluation of that endpoint.

### Consolidated review direction

All seven review-level policy recommendations are accepted. Use DeepSomatic DNA as the starting set for evidence-based integration, RNA as a nomination source requiring DNA verification, and realignment to reassess RNA evidence. Apply the accepted training-label policy to the **final second-round rescue artifact**. First-round outputs remain diagnostic controls; they are not substitutes for a missing or failed second-round training artifact. The historical endpoint derived from realigned RNA remains the baseline when evaluating the revised route; downstream re-consensus and QC derivatives retain that lineage.

The primary superiority comparison is final second-round rescue versus DeepSomatic DNA, separately for SNVs/indels and each input pair in the same declared domain. Compare first-round versus second-round rescue under matched annotation as a mechanistic ablation, recording gained/lost TP and FP. Report the user's lower-FP observation as such until a matched quantitative result is available; also distinguish fewer FP in absolute number from higher precision or improved recall.

Numerical threshold selection, detailed implementation specifications, and workflow execution are subsequent work; this review does not invent their results. No classifier or realignment implementation has been changed.

## Validation and evidence boundaries

Existing scoped classifier/rescue/FILTER contract tests: **50 passed** (15 NumPy deprecation warnings). Synthetic direct-function probes reproduced the support/FILTER defect despite those passing tests. Metrics were checked against native stats/commands; hybrid query contents were matched to production derivatives; exact allele attribution is reproducible from the saved script. Workflow execution and full benchmark reruns were not performed. BAM-level normal contamination, editing overlap, splice-associated FP fractions, and optimized parameter effects remain unmeasured.
