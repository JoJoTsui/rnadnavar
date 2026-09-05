# Hybrid RT/DN routing failure — 2026-09-06

The WES_LL hybrid run completed DNA calling but silently lost the RNA branch
after STAR. A successful Nextflow exit was not evidence of a complete hybrid run.

## Evidence and cause

The trace at
`examples/seqc2/hybrid/output/seqc2.wes.ll.hybrid.modality/pipeline_info/execution_trace_2026-09-05_21-06-33.txt`
contains 62 completed tasks. All three RNA libraries completed STAR and BAM
sorting/indexing, but none reached MarkDuplicates, SplitNCigarReads, BQSR,
RT/DN calling, RNA consensus, or rescue. All three DT/DN callers completed.

In `BAM_ALIGN`, the expected alignment count was calculated per
`(patient, sample, status)` — three for WES_LL_RT_1. The later `groupKey(meta, 3)`
still contained row-specific `lane` and `library`. Each of the three distinct
keys received one BAM, while waiting for three. `groupTuple()` discarded those
incomplete groups on channel closure. The RNA channel was therefore empty before
preprocessing; changing tumor/normal pairing could not restore it.

A real Nextflow operator regression reproduced zero RT pools before the fix.
The same test passed for the original single-lane DN/DT/RT FASTQ triplet.

Two related issues were also verified:

- Matching grouping counts by sample name alone allowed same-named samples from
  different patients to mix. Matching now includes patient and status.
- STAR used `LB:WES_LL_RT_1` for each accession despite distinct manifest libraries.
  Staged raw reads now use the declared library for LB, preserving sample SM.
  Legacy reads retain sample-based LB.

The production-config seq2neo regression additionally exposed an existing
preflight defect: nf-schema emits `[]` for absent optional metadata. Converting
that to a string before checking emptiness made legacy manifests look staged.
Preflight now checks emptiness first; partially staged manifests still fail.

## Correction and boundaries

Alignment grouping excludes row-level lane/library fields and retains a sorted
list of contributing libraries in sample metadata. BAM contents retain their
read groups. Incomplete nonempty pools now fail explicitly rather than vanishing.
The DN/DT caller-ready dictionary audit and caller-pair construction are unchanged.

Neither reference/BED selection nor the original seq2neo production config is
changed. Original input files and previous SEQC2 output directories remain
read-only. Because LB changes the STAR command, old RNA alignments must not be
treated as equivalent cache entries; a verification run must regenerate them.

## Verification

Run the fast regressions with:

```bash
.venv/bin/python -m pytest tests/seq2neo/test_alignment_pooling.py -q
```

These execute production Nextflow channel operators and the actual samplesheet
subworkflow, testing three-library pooling, both caller pairs, legacy triplets,
multiple lanes, patient isolation, split chunks, incomplete pools, library tags,
and absent/partial stage metadata. They do not run biological variant callers.

The bounded integration run uses `tests/run_seq2neo_regression.sh`, which inherits
`examples/seq2neo/seq2neo.shared.config`. Its new output directory is
`.artifacts/seq2neo-hybrid-pooling-regression`, with log
`.artifacts/seq2neo-hybrid-pooling-regression.log`.
End-to-end acceptance requires both modality caller panels, two consensus VCFs,
rescue, and the enabled legacy realignment/annotation stages; launch or exit
status alone is insufficient. Integration results are recorded below when available.

Verification so far: the documented Python suite passed (220 tests); the
additional hybrid completion-gate test passed in the SEQC2 runner's Python
environment. The original pooling reproduction failed with zero RT pools before
the fix and passes afterward. The bounded seq2neo run completed successfully:
167 tasks, both modality caller panels, two first-round consensus VCFs, filtered
rescue, RNA realignment and second rescue, and the enabled annotation stages.
The acceptance validator reported PASS. No external-alignment dictionary audit
or normalization ran on its FASTQ-derived BAMs.

Acceptance trace:
`.artifacts/seq2neo-hybrid-pooling-regression/pipeline_info/execution_trace_2026-09-06_00-53-01.txt`.

For the full WES_LL hybrid retry, use
`examples/seqc2/hybrid/config_pooling_fix.yaml` with the existing SEQC2 driver.
It targets `output/seqc2.wes.ll.hybrid.pooling_fix` and
`runs/pooling_fix_state.json` under the hybrid directory. Its completion globs
explicitly require all six caller outputs, both consensus VCFs, and rescue.
The old DNA-only result cannot satisfy these checks.

From the repository root, after any other run using the same resume session has
finished:

```bash
setsid -f nohup python3 -u examples/seqc2/scripts/run_pipeline.py \
  --config examples/seqc2/hybrid/config_pooling_fix.yaml \
  > .artifacts/seqc2-hybrid-pooling-fix.nohup.log 2>&1 < /dev/null
```

The full retry was launched at 2026-09-06 01:09:10 (local time), run name
`clever_bhabha`. Driver PID 2098364 and Nextflow PID 2098370 were verified alive.
This records launch evidence, not full hybrid acceptance; check its separate
state file, trace, and required outputs for the terminal result.
