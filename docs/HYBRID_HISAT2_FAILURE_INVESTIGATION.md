# Hybrid HISAT2 failure investigation

2026-09-07. Diagnosis/design only: no production implementation or relaunch. Q1–Q8 accepted; final scope confirmation precedes implementation.

## Confirmed cause

The full hybrid run stopped at 11:26:42 after 97 cached and 19 completed tasks. STAR was cached; read extraction and BAM-to-FASTQ completed. HISAT2 never launched (work directory null). Evidence: `examples/seqc2/hybrid/.nextflow.log` and `runs/realign_full_state.json` (mutable run artifacts).

Two independent defects in `subworkflows/local/prepare_reference_and_intervals/main.nf:37–38` violate the consumer declarations in `modules/nf-core/hisat2/align/main.nf:12–13`:

- Index: `.map{[meta,file]}.collect()` flattens eight tuples into sixteen elements. Only metadata and the first index component bind to the consumer.
- Splice sites: `.collect()` emits `[file]`, but the consumer expects `[meta,path]`. Its missing path causes `Path value cannot be null`.

The outer workflow bypasses already tuple-shaped supplied resources from `prepare_genome/main.nf:147–156`. This is not a HISAT2 runtime crash or evidence of another FASTQ mismatch.

## Causal reproduction

Temporary harness `/tmp/hisat2-tuple-repro.lWFJIR/main.nf` uses real resource paths, production producer expressions and actual consumer declarations, with a staging-only body asserting eight components. Nextflow version: 25.10.0.

```bash
cd /tmp/hisat2-tuple-repro.lWFJIR
micromamba run -n nextflow nextflow run main.nf -offline -work-dir work-test
```

| Expression variant | Observed result |
| --- | --- |
| Production | Exact tuple warning and null-path error; reproduced twice |
| `--fix_index true` only | Null-path failure persists for splice input |
| `--fix_splice true` only | One index component binds; eight-component assertion fails |
| Both flags | Exit 0, all eight components and splice sites staged |

Recommendation: a canonical reusable tuple for each supplied/generated reference. `collect(flat:false)` alone would give a list of tuples, not the required single tuple.

## Resource audit

- Eight nonempty `.ht2` components exist with one basename.
- Native `hisat2-inspect-s -s` reports 3,366 contigs: names, lengths and order exactly match the configured FASTA index. Sequence-content identity has not been verified.
- Splice sites: 376,071 rows; 1,246 name contigs absent from the FASTA/index, e.g. `KN196472.1`. The other 374,825 passed field, coordinate-bound and strand checks. Actual HISAT2 handling of the off-reference rows has not been tested.
- Cached aligner version is 2.2.1. The Python `hisat2-inspect` wrapper fails importing removed module `imp`; the native inspection binary works. This diagnostic issue is separate from the run failure.

Proposed policy: create an audited run-local reference-compatible splice subset, preserving source/checksum and reporting excluded rows. Reject malformed/shared-contig coordinates. Regeneration from a verified matching annotation is an alternative. Never edit supplied resources.

## Broader audit

| Finding | Evidence and status | Recommended action |
| --- | --- | --- |
| Completion validator rejects correct scopes | Actual validator reproduced rejecting an in-memory COMPLETED trace with actual workflow scopes. It demands nonexistent suffixes such as `VCF_CONSENSUS_REALIGN` (`validate_realign_run.py:198–202`) | Match qualified second-round scopes and require COMPLETED/CACHED |
| Final artifact path inconsistent | Static chain: RT gains `_realign`, callers add normal pair, rescue combines DNA/RNA pair IDs, VEP publishes with full rescue ID. Validator and runner expect a plain RT/normal ID | One shared explicit artifact identity; do not use arbitrary broad globs |
| Legacy candidate semantics changed | `workflows/rnadnavar.nf:247` globally uses raw consensus; previous filtered VCF excludes noncanonical contigs (`filter_vcf.py:267–271`) | Restore legacy selection and explicitly scope hybrid raw policy |
| Independent library provenance lost | HISAT2 creates one RG without original LB (`hisat2/align/main.nf:35`), then MarkDuplicates reruns (`bam_gatk_preprocessing/main.nf:115`). Duplicate impact not measured | Preserve library-specific reads/alignment/RGs, then pool |
| Multi-patient intervals risk unsorted tabix input | `prepare_intervals/main.nf:139–145` concatenates files without genomic record sorting/merging. Not reproduced in single-patient run | Two-patient overlapping-contig regression; define interval scope |
| Generated reference reuse risk | Generated process outputs are queue channels; multiple RT consumers not tested | Supplied/generated two-sample scheduling tests |
| Weak preflight | Eight-file count does not enforce numbers, basename, assembly or splice agreement | Validate exact resource contract with provenance |
| Large index discovery unsupported | Module searches only `*.1.ht2`, not `.ht2l` | Support or explicitly reject before launch |
| Nonempty BED with zero usable reads | No explicit outcome identified | Distinguish zero candidate loci from zero extractable read pairs |
| Incomplete versions | RNA realignment omits preprocessing versions aggregation | Cover second-pass provenance in tests |

This is a bounded audit, not proof that all future errors have been excluded. Downstream execution, pooled-library duplicate impact and multi-patient behavior remain unverified.

## Safety and decisions

Inspected routing pairs realigned RT with first-pass DN, and second rescue combines first-pass DNA with second-pass RNA. Those downstream paths have not yet completed. No original inputs/historical outputs were modified; no full seq2neo triplet regression was rerun. Behavior-identical seq2neo support cannot yet be claimed.

Accepted decisions Q1–Q4 (user: "all use recommendations"):

1. Repair resource tuples plus completion/artifact contracts and regression-test before rerun.
2. Preserve legacy candidate policy; make hybrid raw-candidate selection explicit.
3. Use an audited run-local splice subset without modifying originals.
4. Require library preservation and multi-sample tests before declaring implementation complete.

Accepted decisions Q5–Q8 (user: "all use recommendations"):

5. Nonempty candidate BED with zero usable paired reads: fail with an explicit diagnostic outcome; do not report zero candidates, promote first-pass output to final, or declare successful realignment.
6. Accept complete single-basename `.ht2` or `.ht2l` index sets, supplied as a directory or file pattern, plus generated indexes. Reject mixed formats, incomplete sets, and multiple basenames before alignment; prove reuse with two RT samples.
7. Preserve current interval membership/scope while sorting and merging equivalent BED coverage before indexing. Do not introduce patient-specific calling-region semantics in this repair; test overlapping multi-patient inputs.
8. Preserve original library identities through extraction/alignment and regroup by biological sample. Fail on ambiguous/missing library identity where provenance is required, rather than inventing independent libraries or silently collapsing known ones. Preserve the legacy single-library path.

All eight policy decisions are settled. The grilling session now awaits confirmation of the consolidated implementation scope, not another policy round. No production relaunch, commit, or push is implied by these design approvals.

Implementation order:

1. Lock down both observed resource-contract failures with regression tests; normalize supplied/generated reference tuples and validate complete index sets and splice resources.
2. Preserve legacy candidate behavior, independent library identities, interval coverage, and explicit zero-read outcomes; test single-library and pooled-library paths and two-patient reference reuse.
3. Align successful process-scope checks, final artifact identity, and second-pass provenance; test valid completion as well as missing/failed/stale artifacts.
4. Run focused hybrid and original seq2neo FASTQ-triplet regressions and report evidence and any remaining limits before a production resume.

Original inputs and historical outputs remain read-only. Do not change first-pass STAR or DNA processing unless a separately demonstrated defect requires a new scope decision.

## References

- [Nextflow collect semantics](https://docs.seqera.io/nextflow/reference/operator#collect): nested lists flatten by default.
- [HISAT2 manual](https://daehwankimlab.github.io/hisat2/manual/): index inspection and splice-site format.
- [Hybrid ingress ADR](adr/0004-converge-hybrid-inputs-at-caller-ready-boundary.md): library boundaries and legacy behavior.
