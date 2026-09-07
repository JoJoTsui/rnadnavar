# SEQC2 work while repaired realignment is pending

Status: waiting-period grilling complete. The user accepted both recommendations: integration first and lightweight execution only. Existing implementation authorization and accepted Q1–Q7 decisions remain in force.

## Established boundaries

The endpoint is final rescue derived from DNA and realigned RNA, including explicitly selected lineage-preserving derivatives. First-round rescue remains a diagnostic control. The user's observation of fewer false positives after realignment has not yet been quantified by the matched completed SEQC2 second-round comparison.

Existing tickets are under `.scratch/seqc2-realignment-label-quality/drafts/issues/`. Their draft headers and unchecked criteria have not been reconciled with implementation commits. A helper or passing unit test is not evidence of complete workflow integration.

## Work executable before final results

| Order | Existing tickets | Reviewable deliverable |
| --- | --- | --- |
| 1 | 01–12 | Evidence-linked status for each acceptance criterion: implemented, validated, remaining, or awaiting real-run evidence. |
| 2 | 02–03 | Executed consensus/rescue controls, zero/missing/conflicting settings, eligible-vote and FILTER/rationale regression cases. |
| 3 | 04–05 | Bounded ingress and second-pass tests: missing/duplicate pairings, empty candidate sets, failed calling, mate/library preservation and original-read recovery integration. |
| 4 | 06–08 | DNA-verification and final-admission integration, including ambiguous evidence and downstream annotation attempts to promote unverified nominations. |
| 5 | 01, 09 | Artifact identity, actual tool/model/database/reference/argument provenance, and explicit lineage/QC selection failures. |
| 6 | 11–12 | Fail-closed policy selection and held-out evaluation fixtures, declared data partitions and report schemas; no empirical policy selection yet. |

Independent analysis can reproduce historical first-round metrics, explain paired error transitions and audit installed arguments/versions. It must preserve inputs and identify the diagnostic stage and benchmark domain.

## Concrete review findings

Inspection of `examples/seqc2/scripts/select_development_policy.py` shows that an empty baseline imposes no slice gates, allowing a candidate to qualify without comparison evidence. Required slices, finite metrics, duplicate rows, data partition identity and immutable policy provenance need validation before this helper satisfies ticket 11.

Inspection of `examples/seqc2/scripts/evaluate_heldout_policy.py` shows that omitted precision/recall gate flags default to true. It does not independently recompute acceptance against a bound baseline or establish that the metrics came from the frozen policy. Ticket 12 therefore remains incomplete even though helper tests pass.

The read-only integration audit found no Nextflow references to `recover_rna_read_pairs.py`, `verify_rna_nominations.py`, or `build_deepsomatic_labels.py` under modules, subworkflows or workflows. Tickets 05–08 still need production-path wiring and observable integration tests.

The initial inspection of `bin/recover_rna_read_pairs.py` identified two pre-integration defects: character-based `rstrip` can corrupt read identifiers ending in 1 or 2, and singleton records are appended to paired FASTQ outputs. Preserve identifiers exactly except for an explicit mate suffix, and keep singleton disposition separate from synchronized pairs.

## Work that needs completed evidence

- Matched historical first/second rescue comparison (ticket 10).
- Measured effect of original-read recovery on final labels.
- Numerical verification/override selection using declared development data (ticket 11).
- Frozen-policy held-out performance acceptance (ticket 12).
- Any claim of improved final-label precision/recall over DeepSomatic.

Do not automatically adopt new training labels or run downstream training; those remain outside the accepted scope.

## Accepted waiting-period decisions

Q1: Prioritize integration and acceptance tests, with lightweight historical diagnostics alongside them.

Q2: Restrict execution to unit tests, tiny synthetic integration fixtures and bounded read-only VCF analysis while the full run is active. Defer additional calling and realignment experiments.

The decision frontier is empty for this waiting-period scope. Scientific policy thresholds and performance conclusions retain their existing real-data dependencies.

No new ADR is needed for reversible work ordering. Existing glossary definitions for training-label artifact, RNA-nominated candidate, benchmark domain and paired error transition apply unchanged.

## Implementation evidence after acceptance

Ticket 05 recovery helper: the numeric identifier and mixed-singleton output defects above are repaired. CLI regression fixtures verify exact mate-suffix removal, synchronized paired outputs, separate singleton outputs, duplicate selected-name rejection, and malformed/truncated FASTQ rejection before outputs are created. Tests: `tests/seqc2/test_read_pair_recovery.py`.

This does not complete ticket 05. Library-scoped candidate extraction, propagation of original FASTQs through both ingress modes, read-length accounting and the indexed second-round demonstration remain outstanding. The CLI requires its caller to supply IDs and FASTQs from the same library; a supplied library string alone cannot prove that relationship.

Ticket 02 numeric controls: module defaults now apply only to absent values, preserving explicit zero floors/minima and allowing invalid zero SNV/indel thresholds to reach CLI validation. The opt-in local fixture `tests/seqc2/test_consensus_module_controls.py` executes actual consensus and rescue processes on tiny synthetic VCFs with containers disabled, checks commands and indexed outputs, and distinguishes default versus zero-floor support and labels. It does not validate root-parameter aliases/conflicts or claim ticket 02 complete.

Validation for recovery commit `c08228a`: documented Python suite excluding the two vcf_stats directories completed with 233 passed, 2 skipped and 27 warnings. Two-axis review found no hard documented-standard violations; spec review confirmed the outstanding ticket 05 integration, library provenance and read-length requirements above.


Policy-selection readiness: the development selector now fails closed on empty or mismatched baselines, undeclared/missing required slices, duplicate slices, nonfinite or out-of-range metrics, malformed partitions and invalid minimum deltas. The held-out evaluator validates frozen slice binding and recomputes per-slice gates from the frozen baseline, ignoring supplied gate booleans. Public CLI fixtures cover qualified, no-policy, missing-policy, malformed and regression cases. This establishes engineering behavior only; real development and held-out evidence remain deferred.
