# Hybrid realignment pairing repair — 2026-09-07

## Root cause

The full run's first-pass RNA consensus exists and contains 835,637 records:
1,293 Somatic, 2,640 Germline, 30,927 Reference, 17,702 Artifact, and
783,075 NoConsensus. Missing candidate records did not cause the failure.

Active caller wrappers and consensus grouping emit pair metadata (`id`,
`patient`, `status`) without the alignment's `sample` field. Preparation joined
that VCF on `(WES_LL, null)` against the RNA CRAM on `(WES_LL, WES_LL_RT_1)`.
The unmatched `join(remainder: true)` result had six fields; the following
eight-argument closure failed before its intended validation could run.

## Repair

- Resolve candidate sample identity at the realignment boundary from explicit
  `sample`/`tumor_id`, or the existing unambiguous `tumor_vs_normal` pair ID.
  Do not rewrite metadata across all callers. Ambiguous identifiers fail.
- Keep patient/sample matching and duplicate rejection. Use `failOnMismatch`
  rather than passing partial tuples into the eight-argument closure.
- Expose raw consensus separately and use it for candidates, preserving the
  accepted all-RT-consensus-record policy. Existing filtered output remains
  available to its other consumers. No rescue-only or DNA-only candidates are added.
- Reuse the configured candidate BED evidence under
  `vcf_realignment/vcf2bed/<sample>/`. Empty BEDs do not enter read
  extraction or second-pass alignment. The completion validator requires a
  successful conversion plus empty BED evidence for zero candidates; absence
  of VCF2BED is a failure, not a zero-candidate outcome.

The candidate route changes apply only when VCF realignment is enabled.
Original seq2neo configuration and initial alignment/caller modules are unchanged.

## Verification

The new Nextflow-channel regression reproduced the exact historical six-field
closure error before the fix. Afterward it passes for hybrid and seq2neo RNA
sample names and rejects missing partners, duplicates, and cross-patient pairing.
Real VCF2BED executions with bcftools verify retention of Somatic and Artifact
records and publication/suppression of a valid empty candidate set.

Passing checks cover `tests/seq2neo/test_realignment_pairing.py`,
`tests/seq2neo/test_alignment_pooling.py`, both hybrid realignment validation/runner
test files under `tests/seqc2`, and `tests/bam_dictionary`.
Tests ran with Python 3.12 where PyYAML is available; the project's Python 3.10
environment skipped the runner file because PyYAML is absent.
A stale pre-existing runner-test filename assertion was corrected.

The Nextflow 25.10.0 preview inheriting `tests/config/seq2neo_regression.config`
completed successfully in `/tmp/hybrid-routing-preview.Vnahj5`.
This is graph/configuration evidence plus focused execution tests, not a new
full seq2neo cohort execution.

The isolated hybrid launcher was resumed with nohup/setsid using the existing
hybrid launch directory, cache, full-run output directory, and append-only log.
Existing source inputs and historical non-realignment outputs were not edited.
Final biological output validation remains dependent on full-run completion.

## Live resume evidence

The September 7 10:15 resume reused all three STAR alignments, including
`82/1a77b1` for SRR9134727. RNA consensus completed at 10:18:53. The repaired
pair reached `PREPARE_REALIGNMENT_VCF:VCF2BED` (`c5/2cfaec`), which completed
at 10:18:59 and published 835,637 BED rows. Candidate-derived interval
preparation completed next. First-pass RNA filtering finished at 10:21:31,
followed by first-pass rescue. This proves the original pairing failure is
cleared on the full inputs, but does not establish final realignment completion.
