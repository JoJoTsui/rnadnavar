# Consensus and rescue runtime and logic review

## Approved optimization follow-up

The user accepted the remaining recommendations. The implemented rules are now maintained in [Consensus and rescue evidence contract](../CONSENSUS_RESCUE_EVIDENCE_CONTRACT.md), with linked updates to the consensus, rescue and RNA-editing guides. Missing AD with a positive floor cannot vote; invalid numeric measurements are unavailable; modality counts distinguish observed/eligible/Somatic; optional verification preserves negative labels; paired-v2 preserves typed per-caller provenance; reporting uses rescued/union and DNA-relative Somatic transitions.

The investigation and first-stage results below are historical evidence, not outstanding approval requests. The five stale compatibility tests have been updated to existing production helper names, version labels and the FAI selection expression; production modules were not changed. Core suite after integrated optimizations: 197 passed. The affected helper test file passes all 9 tests; the preceding broad compatibility run passed the other tests (108 passed, 16 skipped before the last three stale version assertions were corrected).

Optimized validation outputs are isolated under `/tmp/seqc2-evidence-v2.4gYtjL`. Five replay paths on 1,000-record-per-caller WES_LL subsets passed: DNA/RNA/realigned-RNA consensus and both rescue rounds. Each rescue VCF was decoded as paired-v2, retained round identity, satisfied Somatic <= eligible <= observed counts, and carried classification rationale. First rescue contained 5,008 union records; realignment rescue 4,904. All nine source VCF SHA-256 values remained unchanged. These are correctness smoke tests, not a truth benchmark or a complete RNA rerun.

The complete optimized DNA replay also passed: 188,775 variant keys preserved and 2,484 Somatic records unchanged versus the preceding bounded-fix replay. Non-Somatic transitions were Reference->NoConsensus 9, Artifact->NoConsensus 50, Germline->NoConsensus 10, and Germline->Artifact 1. Final counts: Reference 3,645; NoConsensus 171,199; Artifact 3,611; Germline 7,836; Somatic 2,484. Every record retained paired-v2 evidence and classification rationale, with no sample column. These transitions are label changes, not measured FP/FN changes.

Separate findings not changed by the approved Q4-Q8 implementation: standalone reference/header compatibility guards, dormant helper API harmonization, and baseline-retention policy. Their limitations remain documented below; current WES_LL caller inputs use the same configured reference. No improvement in precision/recall/F1 is asserted without a subsequent fixed-domain benchmark.

Scope: current seqc2-consolidated code and WES_LL hybrid failure. User approved bounded implementation defects and isolated caller-VCF replay, with biological-rule changes reviewed separately. Variant-calling cache preservation is required. This review does not claim improved benchmark performance.

Reviewed baseline: `75ba1b5`, followed by the uncommitted bounded fixes described below.

## Confirmed runtime cause and bounded fixes

- `bin/vcf_utils/io_utils.py`, `write_union_vcf`: a missing normal genotype became the string `.` before numeric VAF formatting. The failed DNA task's DeepSomatic VCF has only WES_LL_T_1; its first record is chr1:16957:G:T. Missing normal evidence is legitimate. Keep it as None until serialization. This fixes the shared consensus/rescue writer without changing labels.
- `bin/vcf_utils/aggregation.py`, `extract_genotype_info`: both Strelka SNV and indel count-derived VAF branches omitted VAF_SOURCE. They now report derived, without changing the VAF calculation.
- `bin/run_rescue_vcf.py`: reject repeated modality/caller identities and missing individual caller files before opening output. The same caller in DNA and RNA remains valid. Consensus-only mode continues to ignore individual caller arguments.
- `bin/run_rescue_vcf.py`: accumulate statistics after writing finalizes classification/promotion flags. Previously a promoted record could have RESCUED=YES while the reported rescued count omitted it.

Existing rescue writer tests already exposed the float-formatting failure. Added regressions cover tumor-only normal evidence through both writer modes, Strelka SNV/indel provenance, duplicate paths and caller identities, missing files, and promotion log/output reconciliation.

## Comprehensive logic findings

| Area | Evidence and consequence | Disposition |
| --- | --- | --- |
| Caller independence | Rescue accepted repeated caller identities; list-based votes could satisfy minimum DNA=2 with one repeated caller. | CLI validation fixed; direct internal callers still require valid unique collections. |
| Support semantics | io_utils modality support counts use observed callers, while N_SUPPORT_CALLERS uses eligible voters. annotate_rna_editing consumes modality counts; a rejected DNA record can prevent RNA-only editing treatment. | Policy decision: distinguish observed, eligible, and Somatic support before changing downstream interpretation. |
| Missing alternate depth | aggregation._counts_toward_support explicitly allows absent AD through a nonzero floor. A Somatic record with genotype={} can vote at floor=3. | Existing legacy policy, not silently changed. Recommend explicit unavailable-evidence handling and sensitivity evaluation. |
| Verification and veto | classification.compute_unified_classification_rescue can return NoConsensus for unconfirmed verification before evaluating DNA Artifact. Reproduced DNA Artifact + RNA Somatic + inconclusive => NoConsensus. | Optional-manifest path. Recommend verification gate only for otherwise-Somatic RNA-dependent promotions, preserving established negative labels; approval needed. |
| Source evidence | aggregate_variants concatenates differing source strings with `||`; write_union_vcf only fills absent destination fields. Modality/round identity may be ambiguous and partial reconstructed evidence can mask richer source evidence. | Needs explicit typed evidence merge contract and DNA/RNA/realignment round-trip fixtures. |
| Numeric missing values | Extraction accepts nonfinite AF and negative missing sentinels; aggregation often tests only None. Tumor missing fields use empty entries, normal fields use dots. | Follow-up sanitization at extraction; never substitute zero. Check malformed/missing arrays and finite-value invariants. |
| INFO cardinality | AD strings contain commas and readers can return tuples even for conceptually one serialized field. | Do not blanket stringify tuples or redefine headers without reader round-trip tests. |
| Reference/header compatibility | Streaming chromosome traversal comes from first/DNA template. Other-source-only contigs can be omitted. | Current SEQC2 shares reference; add explicit compatibility validation for standalone inputs. |
| Rescue metrics | Finalized rescued count is fixed, but cross_modality and rescue-rate denominator reflect distinct existing semantics. | Define denominators explicitly before changing metric schema. |
| Baseline retention | Opt-in preserve-baseline-callers returns Somatic before ordinary consensus thresholds. | Retention experiment, not standard >=2 consensus; keep opt-in and report separately. |
| RNA-only labels | Without verification, an RNA-only Somatic consensus can pass rescue. | Decide callable DNA requirements explicitly; absence of a record is not proof of reference. |
| Correlated evidence | ENS_SUPPORT is descriptive; realignment reassesses the same RNA reads. | Do not interpret fixed caller votes or repeated alignments as independent calibrated confidence. |
| Dormant helpers | tagging.compute_unified_filter resolves ties toward Somatic; active consensus resolves ties to Artifact. reclassify_with_annotation has separate COSMIC/common-AF precedence. | No active production caller found in audit; align or deprecate separately, without attributing them to this failure. |

## Cache and compatibility boundary

Only downstream Python consensus/rescue implementation, tests, and documentation are changed. Alignment/caller modules, command arguments, input manifests, reference/BED settings, conda environments, resource configuration, work directories, and Nextflow history are not edited by this fix. The same helpers serve BAM/hybrid and FASTQ modes; no input-mode routing is changed.

This is a source-scope guarantee, not proof of a future cache hit. On eventual resume, use the same successful session/work directory and inspect cached statuses before allowing upstream execution. Imported Python helpers may not independently invalidate successful Nextflow tasks; downstream invalidation must be targeted if needed. Do not force rerun upstream callers or modify their cached outputs to invalidate consensus.

## Validation

Core suite after bounded fixes: 139 passed (54 warnings, including NumPy conversion deprecations).

SEQC2/seq2neo compatibility suites: 106 passed, 16 skipped, 5 failed. All five failures are in unchanged `tests/seq2neo/test_realignment_helpers.py`: three import old process names instead of the current `_AUDITED` exports; two extract the old `splice_fai` assignment while production uses the multiline `audit_fai` selection. The implicated production modules, test file, and reference subworkflow are identical to HEAD. These are pre-existing test/code drift, but they prevent claiming a fully green compatibility suite. No full FASTQ workflow was run.

Real-data replay outputs and logs are isolated under `/tmp/seqc2-consensus-review.lYI8FV`. Consensus replay succeeded for DNA (188,775 union records), RNA (836,095), and prior-run realigned RNA (273,614). These are union record counts, not Somatic counts or TP counts. The realigned inputs come from the prior completed `.realign.full` run, so they exercise downstream compatibility rather than claim completion of the current retry's realignment branch. The full hybrid workflow has not been resumed as part of this audit.

Additional validation:

- All five real-data replay paths completed successfully. SHA-256 checks of all nine source caller VCFs before/after rescue replay were identical. Temporary outputs only were indexed; original cached artifacts were read-only.
- First rescue succeeded with 998,969 records: 3,082 Somatic, 8,642 Germline, 18,513 Reference, 11,977 Artifact, and 956,755 NoConsensus. All records had CLASSIFICATION_RATIONALE and no sample column. Its log reconciles with 573 RESCUED records, including 4 promotions. Its legacy rescue rate is likewise invalid at 100.7%.

- Reversed caller order on a real WES_LL subset produced identical FILTER labels at 2,481 variant keys. This checks label invariance, not INFO byte equality.
- In-memory numeric reproductions confirmed that aggregate_genotypes retains NaN/Inf VAF means and the negative missing-depth sentinel DP=-2147483648.
- Realignment rescue succeeded with 440,008 records: 3,019 Somatic, 8,514 Germline, 10,332 Reference, 7,747 Artifact, and 410,396 NoConsensus. All records had CLASSIFICATION_RATIONALE and no sample column. The emitted RESCUED count reconciles with the log at 628, including 6 promotions.
- The legacy rescue-rate formula reports 101.0% (628/622), because its denominator excludes promoted sites. The corrected raw count is usable; this percentage is not a valid effectiveness metric and must not be used for benchmarking. Define its replacement explicitly.

## Next decision frontier

1. Introduce explicit observed/eligible/Somatic modality counts and select which RNA-editing decisions consume which count.
2. Preserve negative DNA classifications while applying verification only to otherwise-Somatic RNA-dependent results.
3. Define lossless source evidence keyed by modality, caller, sample role, allele, and alignment round; keep unavailable values distinct from zero.
4. Define missing-AD voting and compatible-contig validation separately from runtime serialization fixes.
5. Replace ambiguous rescue-rate reporting with explicitly named fractions and declared denominators; keep truth-based precision/recall/F1 in the separate benchmark domain.

These decisions may change biological labels and require agreement before implementation. Use fixed truth/region definitions and paired FP/FN transitions to evaluate them, rather than optimize thresholds directly against the WES_LL truth set.
