# Hybrid realignment: missing second pass

Investigated 2026-09-07 on `benchmark/seqc2`, HEAD `7f81c94`.
Status: resource-routing fix implemented in the working tree; full hybrid rerun
has not been launched as part of this fix.

## Implemented correction

The supplied-splice audit now resolves the supplied or generated FAI, checks
single-file cardinality, and rejects missing resources explicitly. Audited and
generated splice tuples, and generated HISAT2 indices, are exposed as reusable
values for multiple RNA samples. GTF-derived splice generation is preserved.
Initial BWA/STAR alignment, ingress, and preprocessing code are unchanged.

Six real Nextflow resource integration cases passed: supplied FAI/splice,
generated FAI/supplied splice, generated HISAT2 index, seq2neo-style GTF-derived splice, disabled
realignment, and missing supplied FAI. Both RNA consumers completed in each
enabled positive case. The tests run the production PREPARE_GENOME workflow and
real splice/FAIDX utilities with small fixtures; they do not perform alignment
against a full human reference. See `tests/seqc2/test_reference_resource_routing.py`
for the opt-in invocation. Existing hybrid validation/runner tests also passed
(9 passed, 1 skipped), as did the existing seq2neo alignment-pooling and input
compatibility suite (13 passed, 1 skipped).

The broader completion-validator observations below remain recommendations,
separate from this resource-routing correction.

## Observed outcome

The attempt recorded in
`examples/seqc2/hybrid/output/seqc2.wes.ll.hybrid.realign.full/pipeline_info/execution_trace_2026-09-07_15-45-07.txt`
finished at approximately 17:18 local time with 102 cached tasks, 17 completed tasks,
and zero failed tasks. No HISAT2 task was submitted. The completion validator
correctly rejected missing HISAT2, second rescue, and VEP.

RNA FASTQ reconstruction completed in
`examples/seqc2/hybrid/work/26/bbd1af1e0267277b05e836a2243e5d/`, producing two
approximately 2.3 GB compressed FASTQs. This investigation checked their presence,
not every FASTQ record. There is no new FASTQ mismatch or HISAT2 executable failure
in this attempt: HISAT2 never received all required inputs.

The configured state file is `examples/seqc2/hybrid/runs/realign_full_state.json`,
which records `failed`. The older `runs/run_state.json` belongs to another launcher
configuration and must not be used to report this attempt.

## Confirmed cause

1. `examples/seqc2/seqc2.shared.config:51` supplies `params.fasta_fai`.
   The latest `.nextflow.log:179` confirms the effective value.
2. `conf/modules/prepare_resources/prepare_genome.config:80` disables
   `SAMTOOLS_FAIDX` when that parameter is supplied. Reusing the existing index is
   intentional.
3. `subworkflows/local/prepare_genome/main.nf:150` nevertheless builds the supplied
   splice-site audit input by combining it with `SAMTOOLS_FAIDX.out.fai` directly.
   That skipped process emits no FAI tuple, so the combined channel is empty.
4. `FILTER_HISAT_SPLICESITES` never executes and its output is empty. The trace
   contains neither FAIDX nor the splice audit.
5. `subworkflows/local/prepare_reference_and_intervals/main.nf:43-44` forwards
   the canonical HISAT2 resources, including the empty splice channel.
6. `modules/nf-core/hisat2/align/main.nf:11-13` requires reads, an index tuple,
   and a splice tuple. FASTQs alone cannot schedule HISAT2. All dependent
   realigned-RNA calling, second rescue, and final VEP work consequently remains
   unscheduled.

This is an empty required resource channel, not a valid zero-candidate outcome.
The RNA candidate extraction and read reconstruction ran successfully.

The direct dependency on the generated FAI was added in `9d68937`.
Commit `1755715` then routed explicitly supplied HISAT2 resources through this
audited producer. Together these changes exposed the missing supplied-FAI path.

## Related findings

- Missing `vep` in the hybrid tool list is not the cause: second rescue passes
  `realignment=true` to `VCF_ANNOTATE` (`second_rescue/main.nf:86-91`), whose gate
  allows `vep` OR realignment (`vcf_annotate/main.nf:26`).
- The filtered splice resource is a process-output queue. Generated HISAT2 index
  and splice resources also use raw process-output queues. Shared reference
  resources should be validated and made reusable across multiple RNA samples.
  This is a related multi-sample hazard, not the cause of this single-RT failure.
- `validate_realign_run.py:validate_trace` initially matches caller and consensus
  names across the whole trace. First-pass tasks therefore satisfy those labels
  in its missing-group message. Its later scope check requires some second-pass
  rows, but does not prove every required caller completed in the second pass.
  The printed three missing groups understate the actual missing second pass.
- Direct HISAT2 tests with generated resources and one sample do not exercise
  the supplied-FAI canonical resource route that failed here.

## Recommended repair and verification

Resolve a single effective FAI from either the supplied file or generated output
before splice auditing. Keep audit enforcement; do not bypass it or rebuild a
valid supplied index merely to supply a channel event. Validate mandatory resource
presence/cardinality and expose shared resources as reusable values, with explicit
errors for missing required resources.

Regression checks should exercise supplied and generated FAI, supplied/generated
HISAT2 resources, and at least two RNA samples sharing references. Include missing
resource failure, the real SEQC2 parameter combination, and the seq2neo raw-read
DN/DT/RT path with realignment disabled and enabled. Use small fixtures and avoid
writing to original inputs or existing non-hybrid SEQC2 outputs.

Require each second-pass caller/consensus/rescue/annotation group in its expected
workflow scope and bind completion evidence to the current attempt. Preserve the
distinct valid zero-candidate contract.

After those checks pass, resume the hybrid run and verify actual HISAT2 submission,
second-pass caller completion, second rescue, and the final annotated VCF. Until
then a further unchanged resume will encounter the same empty resource channel.

## Domain terminology

Existing glossary terms in `CONTEXT.md` apply unchanged: **RNA realignment**,
**Realignment candidate locus**, and **Second-round rescue**. Workflow-engine
completion alone is insufficient evidence that second-round rescue exists.
No new domain term or architectural trade-off is required by this diagnosis.
