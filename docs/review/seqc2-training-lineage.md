# Training-label lineage clarification

The user clarified during the review that downstream training uses only rescue results derived from realigned RNA, because they have observed fewer false positives in those results. A bounded, read-only check of this repository's manifests supports the input-lineage statement. External training executions were not inspected.

## Verified local selection

Paths and line numbers refer to the original repository at `/t9k/mnt/hdd/work/Vax/pipeline/rnadnavar`; processed manifests are local runtime artifacts, not reproduced into the review worktree.

- `examples/seq2neo/build_sample_manifest.py:34–37,160–161` constructs rescue paths under `vcf_realignment/rescue` with `RT_realign` in the identity; RNA caller paths at lines 49–67 also use realigned outputs.
- `examples/seq2neo/data/processed/sample_manifest.tsv` has 66/66 rescue paths and RNA DeepSomatic inputs with realignment provenance.
- `examples/seq2neo/data/processed/sample_manifest_rerun.tsv` retains realigned-RNA lineage for all 66 rescue entries. It has 63 nonempty `training_label_vcf` entries: 56 rerun rescue VCFs and seven QC-cleaned derivatives. Three samples have no training path.
- `examples/seq2neo/docs/data_split_strategy.md:73–82` identifies that rerun manifest's `training_label_vcf` as the source of truth, superseding historical `training_manifest_65.tsv`. The latter also appears in the older runbook state at `docs/RECONSENSUS_RERUN.md:309`.
- `examples/seq2neo/scripts/build_rerun_manifest.py:85–87` publishes rerun labels under `output_reconsensus/.../rescue`, while input lineage and `RT_realign` identities remain. A directory named `rescue` alone does not establish a first-round input.

## Helper policies and their limits

`examples/seq2neo/config/runner.yaml:49–53` permits ordinary or realignment rescue to satisfy completion, and `run_batch_from_json.py:215–221` uses any matching rescue plus failed-trace checks. This completion policy is broader than the stated training endpoint, but is not evidence that first-round labels reached training: the inspected manifest entries retain realigned-RNA lineage.

`build_rerun_manifest.py:93–109,127–138` can fall back to uncleaned rerun labels when QC is absent or a cleaned WARN derivative is missing; FAIL exclusion depends on a hardcoded sample exclusion list. The inspected manifest has populated QC and its FAIL sample excluded. These are prospective provenance/QC-enforcement concerns for later implementation review, not observed first-round training contamination.

## Consequence for this review

Use final rescue labels derived from realigned RNA (including their selected re-consensus/QC derivatives) as the training endpoint and the primary comparison against DeepSomatic DNA. Retain first-round rescue as an intermediate diagnostic. The measured first-round SEQC2 precision loss must not be presented as the precision of the current training-label artifact.

No matched first-versus-second precision/FP comparison was readily found in the scoped local documentation and review evidence. The lower-FP statement is retained as the user's production observation. Quantifying it requires comparable annotation, selection, truth and regions, reporting both absolute FP changes and precision/recall. The failed SEQC2 realignment attempt supplies no final second-round score.
