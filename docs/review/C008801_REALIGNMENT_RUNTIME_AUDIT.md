# C008801 realignment runtime audit — 2026-09-08

Scope: persistent integration worktree, branch `integration/seqc2-hybrid-realign`, C008801 smoke execution. This is a runtime and integration audit, not a SEQC2 precision/recall result. Original reference files, GTF and input manifest remain unchanged.

## Run evidence

Output root: `/t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/.artifacts/c008801.realign.full`.
Resume trace: `pipeline_info/execution_trace_2026-09-08_08-32-53.txt`.

- Both first-round DeepSomatic tasks are CACHED with exit 0; DNA and RNA consensus and first rescue are also cached.
- `VALIDATE_READ_IDS`, `SORT_MERGE_BED`, enhanced CRAM conversion and interval compression completed with exit 0 after `a56dabf`.
- The read-ID list contains 373,067 rows / 231,627 unique names. These are not verified complete pairs.
- Merged candidates contain 25,480 intervals covering 161,636 bases.
- Converted RNA BAM is 208,465,799 bytes. Samtools reports version 1.21; output quickcheck and indexing succeeded.
- At inspection, first-rescue COSMIC/gnomAD annotation was RUNNING and Picard FilterSamReads was NEW in the submission queue. The defaults request 64 GB and 56 GB respectively; those reservations explain serialization under the approximately 78 GiB cgroup limit. Lowering them requires measurement, not an assumption that the pipeline is hung.

## Reproduced defects and repairs

### RTA-01 — supplied FASTA index silently suppresses splice-site processing (P0)

`conf/modules/prepare_resources/prepare_genome.config` skips `SAMTOOLS_FAIDX` when `params.fasta_fai` is supplied. `PREPARE_GENOME` nevertheless combined supplied splice sites exclusively with `SAMTOOLS_FAIDX.out.fai`. The real C008801 config supplies its FAI; therefore the combination emits nothing and HISAT2 cannot receive its required splice-site tuple.

A Nextflow fixture executing the production combination emitted `COUNT=0` for a supplied FAI. Generated-FAI input emitted one tuple. Repair: use the configured FAI when present, otherwise the generated FAI. The regression requires exactly one tuple for both paths. Missing configured FAI fails through `checkIfExists`.

### RTA-02 — third version-heredoc failure (P0)

`FILTER_HISAT_SPLICESITES` retained the same escaped `\${task.process}` that broke the two helpers fixed in `a56dabf`. Once RTA-01 is repaired, this would become a hard failure. A real module execution with a valid tiny FAI and splice row reproduced Bash `bad substitution` after successful filtering. Repair: interpolate the process name in Nextflow. The test also checks retained data and nonempty version content. A tracked-source search found no further occurrences of this exact escaped expression in local modules/subworkflows after repair.

### RTA-03 — disk estimates measure symlinks rather than input files (P1)

The enhanced conversion log reports a 129-byte CRAM and estimates 258 bytes of required disk. The actual CRAM is 185,150,614 bytes. Its `stat -c%s` statements measure Nextflow's staged symlinks. A fixture executing the production input-size shell assignment measured 153 bytes for a staged 4,096-byte target. Repair: use `stat -Lc%s` for file-size measurements. The remaining two-times-input estimate is still a heuristic, not an upper bound on CRAM-to-BAM expansion.

## Open findings requiring separate acceptance evidence

| ID / prior ticket area | Priority | Evidence and consequence | Acceptance criterion |
| --- | --- | --- | --- |
| RTA-04 / 04–05 | P1 | `VALIDATE_READ_IDS` runs only `test -s`; extraction is `samtools view -L BED ... | cut -f1` with empty extra arguments. Neither proves that a selected name has both usable primary mates. | Validate synchronized nonempty R1/R2 outputs after recovery; report paired, singleton and missing counts; demonstrate a nonempty-ID/no-pair fixture fails before HISAT2. |
| RTA-05 / 05 | P1 | No module/workflow invocation of `recover_rna_read_pairs.py`. Active route converts the already preprocessed RNA CRAM, filters by names, then converts to FASTQ. Original sequence and per-library identity have not been demonstrated at the second alignment. | Wire library-scoped original-read recovery and prove exact sequence/mate/library preservation. Quantify any difference from the current reconstructed FASTQs. |
| RTA-06 / 06–08 | P1 | No production invocation of `verify_rna_nominations.py` or `build_deepsomatic_labels.py` in this branch. Their standalone tests do not establish final-label admission behavior. | Integrate verification and admission with final second-rescue lineage, retaining per-allele evidence and testing unverified nominations cannot become training Somatic labels through downstream annotation. |
| RTA-07 / 04,09 | P1 | `VALIDATED_SAMTOOLS_CONVERT` logs “Skipping validation for non-file object: UnixPath”; it checks `hasProperty('exists')`, then catches validation exceptions as warnings. Outer conversion checks did succeed for this sample. | Exercise actual Nextflow Path objects, missing/unreadable inputs and fail-closed error propagation; do not silently drop bad samples. |
| RTA-08 / 09 | P1 | Completion validator checks caller-name presence globally. It requires some realignment-scope rows but does not require each caller's successful rows in that scope. | Remove each second-round caller from a fixture while retaining first-round callers and assert failure. Verify status, sample identity and artifact lineage for each required stage. |
| RTA-09 / 09 | P2 | gnomAD annotation logs no resource for chromosome M and passes those records through unannotated. | Report unsupported contigs and absent AF evidence explicitly; do not treat absence of annotation as demonstrated rarity. Keep the user-selected seq2neo resources. |
| RTA-10 / 09 | P2 | RNA realignment omits `BAM_GATK_PREPROCESSING.out.versions` from its version aggregation. Toolchain mixes samtools 1.15.1 for extraction, 1.21 for conversion/FASTQ, and 1.20 in the HISAT2 environment. | Capture actual per-stage versions, arguments, model and input/reference identities; assess compatibility from execution evidence. Version differences alone do not prove an error. |

Existing ticket numbers refer to `docs/review/SEQC2_WAITING_PERIOD_PLAN.md`. The original `.scratch/.../drafts/issues` directory was unavailable in this integration worktree during this audit; this tracked table preserves actionable criteria rather than claiming those drafts were updated.

## Tool arguments and operational interpretation

HISAT2 is pinned to 2.2.1 with samtools 1.20. Its paired path sets `--no-mixed --no-discordant` and pipes to `samtools view -bS -F 4 -F 8 -F 256`. These selections can reduce available RNA evidence; their contribution to precision/recall needs paired read accounting and matched completed results. No policy change is made solely from this inspection.

Picard is pinned to 3.3.0 and uses `includeReadList`, `--VALIDATION_STRINGENCY LENIENT`, and `--CREATE_INDEX true`. The FASTQ conversion uses samtools 1.21, separate mapped/unmapped subsets, `-N`, and separate singleton outputs. The active HISAT2 route does not consume singleton FASTQs.

DeepSomatic is configured with four shards and a 64 GB memory reservation. The current successful DNA task's trace reports 120 GB peak RSS despite the approximately 78 GiB cgroup cap, so that trace value must not be interpreted as unique physical memory. Cgroup `oom_kill=4` is cumulative and cannot attribute a particular historical failure without timestamped evidence. At inspection, roughly 63.6 GB of cgroup usage was cache and 19.4 GB RSS. Host `free` reports 503 GiB and is not the container's usable limit.

The output filesystem had approximately 59 GiB free (rounded utilization 100%). This is an observed capacity constraint, not a current ENOSPC failure. No active work or user files were removed.

VEP is deliberately invoked in this integration branch for `realignment=true` even without `vep` in `--tools`; the completion contract expects a stripped eight-column VEP rescue VCF. This is internally consistent here, although the implicit CLI behavior should be made explicit in future interface work.

## Validation and rollout

`tests/seq2neo/test_realignment_helpers.py` executes production Nextflow helper modules, supplied/generated FAI channel selection, and the conversion's shell size assignment. Before fixes: splice version reporting and supplied-FAI cases failed; the symlink fixture separately failed. After the first two fixes all five initial fixtures passed. Final helper plus completion-validator results are recorded below.

Small test artifacts are under `.artifacts/realignment-audit/` in the persistent integration worktree. The running full workflow has already parsed its modules, so these source repairs require a later `-resume` launch to take effect. Leave the current run's useful work intact; do not start a competing run into the same work/output directories. Final realignment, second-round callers and training-label validation remain pending.

Final validation: `micromamba run -n nextflow /t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/.venv/bin/python -m pytest tests/seq2neo/test_realignment_helpers.py tests/seqc2/test_hybrid_realign_validation.py -q --basetemp .artifacts/realignment-audit/verified` completed with **15 passed in 25.78s**. `git diff --check` passed.

The last live process check showed four `bcftools annotate` workers on different chromosomes with active CPU use, including workers launched after earlier chromosomes completed. The gnomAD task is making progress. This does not remove the independently reproduced empty splice-site channel in the already parsed full run.
