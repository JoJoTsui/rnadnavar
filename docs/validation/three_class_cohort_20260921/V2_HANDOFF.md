# v2-only audit and three-sample handoff

All derivatives use the already completed `separated_three_class_v2` final
realignment-rescue export. No mapping, calling or new consensus is performed.
The original 66-sample export remains untouched and is **not** an inclusion
list for training: it contains three historically excluded samples.

## Sources and outputs

Shared root:
`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919`.

Original v2 inventory:

- `variant_parquet/three_class_candidates.parquet`: 15,862,893 candidate rows.
- `variant_parquet/refined_three_class_manifest.parquet`: 66 sample rows.
- `exports/manifest.tsv`: the same 66 sample IDs and v2 rescue paths.

New review-only handoff under `handoff_review_20260921/`:

- `cohort63.review.parquet` and `cohort63.review.tsv`.
- `three_samples.review.parquet` and `three_samples.review.tsv`.
- `audit.json`: input/output hashes, per-sample counts and 48 source-record checks.

| Subset | Samples | Somatic | Germline | Reference | Total |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original v2 | 66 | 32,933 | 1,585,619 | 14,244,341 | 15,862,893 |
| Quarantined review subset | 63 | 28,241 | 1,522,528 | 12,791,828 | 14,342,597 |
| Direct three-sample stage | 3 | 1,017 | 63,350 | 912,310 | 976,677 |

The one-sample 4007 stage is skipped. The three-sample stage still contains
`PRJNA298376_4007`, `PRJNA298376_4060`, and `PRJNA298376_4072` together.
They are all `train_pool` in the existing shared `sample_split.tsv`. Historical
reserved pools must remain reserved; the 63-sample file is not a declaration
that all 63 samples belong to training. No new train/validation/test assignment
is made by this handoff.

| Sample | Somatic | Germline | Reference |
| --- | ---: | ---: | ---: |
| PRJNA298376_4007 | 569 | 21,855 | 50,534 |
| PRJNA298376_4060 | 213 | 19,759 | 278,568 |
| PRJNA298376_4072 | 235 | 21,736 | 583,208 |

## Exclusions and current status

| Historically excluded sample | Current v2 Somatic / Germline / Reference | Handoff |
| --- | --- | --- |
| PRJNA298330_4032 | 397 / 19,194 / 68,009 | Entire sample excluded |
| PRJNA298376_4081 | 286 / 18,241 / 1,298,064 | Entire sample excluded |
| PRJNA298376_4255 | 545 / 25,656 / 86,440 | Entire sample excluded |

All three have completed v2 outputs; none has a new re-admission decision.
The historical causes and exclusion decision are in
`docs/RERUN_LABEL_ANOMALY_ANALYSIS.md` and
`examples/seq2neo/docs/data_split_strategy.md`. Current v2 counts are not the
old count anomaly; successful regeneration and native evidence spot checks do
not clear possible biological/sample issues. Their rows are absent from both
handoff Parquets and TSVs.

Within the remaining 63 samples, 3,464 Somatic rows with explicit annotation,
native-negative or DNA-verification review conflicts are quarantined rather
than relabeled. This is a **proposed conservative training subset**, not a
changed benchmark policy or a claim that every quarantined variant is false.
The same operation withholds 96 Somatic rows across the three-stage samples.
Native Germline/Reference candidates are retained under v2's actual rules;
the stricter 299-read/zero-ALT pilot gate is not imposed. Negative sample volume
is large, but coverage filters and balanced sampling still determine usable
model examples; quantity does not certify label accuracy.

## Audit evidence and limits

- Bound source Parquet, manifest and summary to completed export verification.
- Checked 66 state/sample/policy/final-path bindings.
- Checked native-negative class/provenance consistency across the source rows.
- Reproduced 48 native negative nominations from indexed DeepSomatic records:
  four per class for each of the three stage samples and three excluded samples.
  Checked tumor sample headers, GQ/DP/AD/PL under the existing native rules.
  This deterministic prefix check tests plumbing, not representative accuracy.
- Independently reread both output Parquets to verify every sample/class count
  against the audit and confirm zero rows from excluded samples.
- The sample-level original Parquet and TSV agree for all 66 IDs and v2 paths.

All source approval fields remain unchanged (`training_eligible=false`,
`TRAINING_ELIGIBLE=NO`), and the new manifests have no `training_label_vcf`.
The row-selected Parquet must not be replaced with its unfiltered source VCF.
Status is `review_ready_not_training_approved`: this is a reviewable proposal
for a workflow-derived weak-label release, not independently established truth.
An explicit scope approval and EvoSomatic loader/cache/split check are needed
before launch. No model training or remote manuscript modification was performed.
The full 299-read biological validation is **not** required for that scoped
approval. No audit can guarantee zero biological label errors.

## Reproduction

Choose a fresh destination:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 prlimit --as=17179869184 -- \
  .venv/bin/python -u examples/seq2neo/scripts/prepare_v2_training_handoff.py \
  --source-root /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919 \
  --outdir /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/handoff_review_20260921
```

No original input/output is rewritten; fresh review outputs are bounded-batch
Parquets. A failed run's files must not be consumed: require the passing audit
status and its matching output hashes.
