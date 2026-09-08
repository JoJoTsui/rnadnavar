# SEQC2 consensus and rescue evidence optimization

**Status:** draft — test boundaries pending confirmation
**Tracker:** local files

## Problem Statement

The current WES-LL integration outputs underperform DNA DeepSomatic. Consensus loses useful baseline calls while RNA-only rescue introduces many nonmatching variants. Ambiguous INFO support names and incomplete paired-sample evidence undermine explanation and safe optimization. Historical metrics and current source must be distinguished.

## Solution

Make per-caller sample evidence auditable from caller input through DNA consensus and both rescue rounds, then evaluate an explicit DeepSomatic-preserving strategy with DNA-verified RNA nominations. Retain existing labels and default FASTQ behavior. Supply reproducible benchmark outcomes and evaluate frozen improvements before adoption.

## User Stories

1. As a pipeline user, I want DNA consensus and rescue to improve on matched DNA DeepSomatic, so that integration adds measurable value.
2. As a model developer, I want biological FILTER labels to remain stable, so that existing label consumers continue working.
3. As a analyst, I want Somatic-only selection without obsolete RaVeX flags, so that the benchmark measures the intended product.
4. As a analyst, I want SNP, INDEL and records precision, recall and F1, so that pooled scores cannot conceal regressions.
5. As a analyst, I want retained per-allele TP/FP/FN and uncertainty, so that I can trace performance changes.
6. As a analyst, I want explicit region-overlap and normalization semantics, so that boundary and representation differences are not biological errors.
7. As a analyst, I want paired gains and losses between callsets, so that equal TP totals do not hide replaced true variants.
8. As a analyst, I want gnomAD from annotated rescue joined to caller evidence, so that annotation-stage differences are interpreted correctly.
9. As a analyst, I want missing database evidence distinguished from zero AF, so that unannotated sites are not assumed rare.
10. As a pipeline user, I want consistent tumor and normal sample resolution, so that classification and metrics describe the same samples.
11. As a analyst, I want full allele-specific AD and depth retained by caller, so that support can be checked against original observations.
12. As a analyst, I want caller AF distinguished from count-derived fractions, so that different estimators are not incorrectly treated as contradictions.
13. As a pipeline user, I want explicit unavailable normal evidence for tumor-only callers, so that missing data cannot falsely confirm a nomination.
14. As a pipeline user, I want safe multiallelic, missing-value and duplicate handling, so that alleles and evidence are not silently overwritten or swapped.
15. As a analyst, I want observed, eligible and Somatic support defined separately, so that similarly named INFO fields do not mislead downstream analysis.
16. As a model developer, I want compatible versioned INFO aliases, so that existing consumers survive the evidence expansion.
17. As a pipeline user, I want evidence preserved through consensus-only and full rescue, so that reruns do not lose source observations.
18. As a pipeline user, I want audited DeepSomatic retention and overrides, so that weaker votes do not automatically discard strong baseline calls.
19. As a pipeline user, I want RNA nominations verified against DNA tumor and normal, so that RNA-only artifacts do not automatically become Somatic labels.
20. As a pipeline user, I want confirmed, rejected and inconclusive verification outcomes, so that insufficient evidence is not mistaken for absence or confirmation.
21. As a pipeline user, I want correct insertion and deletion verification, so that indel performance is not sacrificed for SNP improvements.
22. As a pipeline user, I want realignment treated as reassessment of the same reads, so that correlated observations are not counted as independent support.
23. As a pipeline operator, I want hybrid-scoped upstream STAR MAPQ60, so that the upstream convention is explicit without changing original seq2neo defaults.
24. As a pipeline operator, I want original FASTQ routing and hybrid BAM checks preserved, so that both input modes remain usable.
25. As a pipeline operator, I want isolated runs with input checksums, so that original data and historical outputs remain unchanged.
26. As a researcher, I want frozen rules and held-out evaluation, so that reported improvement is not threshold overfitting.
27. As a researcher, I want truth/database conflicts retained for investigation, so that population AF does not silently rewrite the truth.
28. As a maintainer, I want already-fixed defects separated from remaining problems, so that implementation avoids redundant fixes.
29. As a researcher, I want failed and inconclusive experiments reported honestly, so that adoption follows evidence rather than a promised score.

## Implementation Decisions

- Review baseline is commit 152d304. Existing support-floor fixes are validated, not reimplemented; historical labels are reported separately.
- Preserve eight-column VCF and biological FILTER vocabulary. Select FILTER=Somatic for workflow benchmarks, native caller PASS for controls; ignore RaVeX_FILTER. PASS conversion is benchmark-only, without -P.
- Use a shared validated sample/allele evidence model keyed by caller, sample role, modality and alignment round. Preserve full available raw tumor/normal AD, DP, GT, AF, quality and caller diagnostics; keep derived fractions distinct. Never sum depth across callers.
- Distinguish absent, unavailable, invalid and zero evidence. Conflicting roles, truncated sample arrays, duplicate conflicting keys and unsupported allele mappings fail clearly. Do not synthesize normal evidence from tumor-only DeepSomatic.
- Expand canonical INFO alongside compatible legacy aliases. Define observed, eligible-detection and affirmative-Somatic counts separately; do not silently rename or reinterpret them. Retain available evidence through full and consensus-only rescue.
- Join gnomAD and other database annotations at the annotated rescue stage with source version/reference/allele provenance. Missing annotation is unknown; no unconditional truth override or locus-specific production exception.
- Keep legacy strategy as default control. Experimental DNA integration starts with DeepSomatic PASS; additions/removals require explicit evidence rules, not majority or caller absence alone.
- Follow the accepted RNA-nomination ADR: verify original DNA tumor/normal evidence; inconclusive nominations remain outside Somatic. Realignment is correlated reassessment, not another vote. Preserve mate/library identity and candidate provenance.
- Verification covers SNPs and indels, local allele context, read-filter provenance and input checksums. Existing standalone helpers need validation before reuse.
- Expose an optional STAR unique-MAPQ control; hybrid uses 60, unset retains the original command. Retain SplitNCigarReads and test intermediate effects.
- Preserve normalized benchmark outcomes and source-allele links, correct HC/target overlap, other/complex accounting, uncertainty and paired transitions. Haplotype-aware adjudication remains distinct from exact diagnostics.
- Develop a bounded rule grid outside chr1, then freeze code, parameters and source versions. Evaluate each product/input pair separately. Require higher records F1 without lower precision; expose and reject hidden per-type regressions. Report stricter simultaneous SNP/INDEL precision-recall superiority separately. Replicate SEQC2 data demonstrate robustness, not independent biology.

## Testing Decisions

- Primary existing boundary: caller/consensus VCF inputs through the consensus and rescue command interfaces to emitted VCFs. Assert observable labels, evidence, headers, errors and round trips, not internal helper calls.
- Benchmark command boundary: small truth/query/region fixtures through retained outcomes and final metric tables. Include HC-overlap deletion, equivalent allele representations, missing annotations and source selection.
- Existing bounded Nextflow workflow boundary: original DN/DT/RT FASTQ triplet and hybrid caller-ready DN/DT plus RT FASTQ inputs through expected artifacts. Verify BAM checks, option propagation, both rescue rounds, empty outcomes and original defaults.
- Reuse existing consensus module-control, policy/provenance, read-pair recovery, label-QC, hybrid runner/validation and seq2neo regression test patterns. Narrow parser tests supplement these only where malformed evidence is difficult to exercise at the higher boundary.
- Each ticket carries its own behavior tests. Held-out experiments are evaluation work, not a substitute for implementation tests.

## Out of Scope

Guaranteed superiority; rewriting truth; truth membership in production decisions; unrestricted tuning; new model training; global aligner upgrades; deleting legacy INFO; changing original inputs or outputs; default strategy or shared STAR changes before reviewed adoption; external issue publication.

## Further Notes

The audit found DeepSomatic-to-consensus changes of eight gained and 71 lost truth matches with 80 added exact nonmatches. RNA-only pass-through accounts for 397 of 401 realignment-rescue added nonmatches. Exact diagnostics are not haplotype-aware FP adjudication. The gnomAD/truth conflict at chr14:106324560 requires locus/source evidence, not a blanket assumption that truth is wrong. Numeric verification rules and actual coverage causes remain bounded experiment deliverables.

