# Hybrid second-rescue VEP omission

Verified 2026-09-08 after committing the reference-resource correction as
`997d66b`.

## Observed failure

The full hybrid trace `execution_trace_2026-09-07_19-33-03.txt` contains successful
second-rescue filtering but no `ENSEMBLVEP_VEP` task. The stripped intermediate
`*.filtered.vcf.stripped.vcf.gz` exists. The missing artifact is its VEP-annotated
successor, `*.rescue.filtered.stripped.vep.vcf.gz`. The runner correctly rejects
the incomplete result despite Nextflow exiting successfully.

## Verified cause and correction

The recorded `params_2026-09-07_19-33-51.json` has `vep_genome=null` and
`vep_cache_version=110`. `VCF_ANNOTATE` replaced a missing assembly with
`Channel.empty()`. Because assembly is a required VEP process input, no task could
be scheduled. Second rescue enables annotation through its `realignment=true`
argument even when `vep` is absent from `params.tools`.

The earlier explanation that the first annotation branch drained the VEP cache
was incorrect. `collect()` produces a reusable value, and DSL2 supports multiple
consumers. A production-workflow regression confirms that two annotation branches
each process two samples with the same unmodified collected cache channel.
See [Nextflow channel documentation](https://docs.seqera.io/nextflow/reference/stdlib-types/channel)
and [DSL2 channel reuse](https://docs.seqera.io/nextflow/migrations/dsl1).

The hybrid realignment config now explicitly selects `GRCh38`, `homo_sapiens`,
and cache version `115`, matching seq2neo's existing settings. Read-only inspection
confirmed `homo_sapiens_merged/115_GRCh38/info.txt` under the configured cache root;
its metadata identifies assembly GRCh38 and species homo_sapiens. The inherited
VEP arguments already contain `--merged`. Keeping version 110 would have exposed
a second error after scheduling was repaired.

`VCF_ANNOTATE` now rejects missing assembly, species, or cache-version settings
with an explicit parameter error whenever annotation is enabled. Annotation-off
runs retain their previous behavior. The cache channel and alignment code are
unchanged, as is the seq2neo configuration.

## Regression evidence

Before the fix, this command failed because the missing-assembly case exited
zero and printed `ANNOTATED_SAMPLES=[]`:

```bash
PATH=/t9k/mnt/joey/micromamba/envs/nextflow/bin:$PATH \
RUN_NEXTFLOW_VEP_TESTS=1 .venv/bin/python -m pytest \
  tests/seqc2/test_vep_routing.py -q -k 'missing_genome or seq2neo or hybrid'
```

Both configured positive controls passed before the fix. After the correction,
all six annotation cases pass: hybrid forced annotation, seq2neo-style explicit
VEP with shared cache across both branches, annotation disabled, and missing
assembly/species/version errors. These tests use the production annotation
subworkflows and module command with fake VEP/tabix executables; they verify
scheduling and output routing, not biological annotation accuracy.

The combined annotation/runner/completion-validator suite passes 15 tests with
one skipped. The preceding reference-resource suite passes all six cases when
the required HISAT2 tools are on PATH. The seq2neo alignment-pooling and
realignment-pairing suites pass 21 tests with one skipped.

A real-executable smoke test in `/tmp/hybrid-vep-smoke.LJQJaD` also completed
both production VEP tasks using VEP 115.2, the installed merged 115_GRCh38 cache,
and the configured GRCh38 FASTA. Its input is a 263-record chr1 subset of the
existing stripped rescue VCF, copied with bcftools into the temporary directory.
Both sample outputs were emitted successfully. This smoke uses minimal offline
merged-cache arguments; it does not exercise every optional production annotation
flag or the entire full-run VCF. Source inputs, cache, and published hybrid outputs
were not edited.

Full hybrid resumption is a separate run;
test completion does not imply the final full-run VEP artifact exists.

## Completion-validator hardening

The hybrid validator now matches exact required task names in their intended
preparation, RNA realignment, or second-rescue scopes. First-pass tasks and
downstream caller postprocessing cannot substitute for the actual second-pass
caller. Every matching required task must be COMPLETED or CACHED; ABORTED,
RUNNING, and SUBMITTED rows are rejected. The final output identity must match a
successful second-rescue VEP task, rather than any occurrence in the trace.

Final VCF validation now requires a CSQ INFO format header, eight columns on
every record, and nonempty CSQ annotation for every retained variant. An empty
VCF with the proper annotation header remains valid. Gzip and tabix checks remain
in place. This is content and task-identity validation, not a checksum-based proof
that an artifact came from a specific work directory.

The focused validator/runner suite passes 32 tests with one skipped. In-memory
mutations of the real hybrid trace now reject missing second-pass DeepSomatic
and aborted VEP while accepting the otherwise complete trace with a successful
VEP task added. Both real VEP smoke outputs pass the hardened content/index check.
The full production annotation arguments also passed on the 263-record fixture,
preserving all original sites, FILTER values, and INFO fields. No seq2neo workflow
or configuration changes were needed for these validator corrections.
