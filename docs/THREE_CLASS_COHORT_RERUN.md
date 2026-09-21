# Candidate-only three-class cohort rerun

Prepared after the [completed three-dataset validation](validation/separated_three_class_v2_20260919/README.md).
The three-sample pilot completed on 2026-09-19. The full 66-sample candidate
cohort finished on 2026-09-21 with no failed samples. See the
[cohort/export review](validation/three_class_cohort_20260921/README.md) for
verification status and output paths, and the
[pilot review](validation/three_class_pilot_review_20260919/README.md) for pilot evidence.
The new policy is `separated_three_class_v2`; the old cohort configuration and
production workflow defaults remain unchanged.

## What runs

```text
Existing DNA caller VCFs ----------------> established refined DNA consensus
Existing realignment RNA callers/rescue -> established gated realignment rescue
Existing DNA caller VCFs ----------------> independent native negative nominations
                                               |
                         separated three-class union, independently per stage
                                               |
                         structural audit + source/code/output integrity
                                               |
                         candidate VCFs and candidate manifest (not truth)
```

Only VCF normalization, consensus, rescue and candidate labeling run. There is
no Nextflow launch, mapping, variant calling or cache deletion. Existing caller
VCFs, original rescues, BAMs, references, indexes and prior cohort outputs are
read-only. No fresh VEP run is performed: available annotations are inherited;
missing annotations/evidence stay missing. The cohort's final rescue uses the
**realignment** RNA panel and realignment rescue, not first-round rescue.

The established Somatic algorithms are run without the experimental three-class
flag. A separate invocation uses that flag only to obtain native negative
nominations. The separated adapter preserves each stage's exact Somatic set.
Do not substitute the old `apply_refined_rescue.py --experimental-three-class`
experiment: its Somatic-retention failure is documented in the validation archive.

## Gate and limits

Preparation verifies the committed validation summary, its archived evidence
hashes, and the policy implementation hashes. Missing stages, stale code,
changed evidence or non-candidate mode fail closed. All 66 input sets are checked
for readable indexed VCFs, reference/header compatibility and native DNA
DeepSomatic GQ schema. These are header/first-record checks; full parsing and
source checksumming happen at execution.

All candidates remain `TRAINING_ELIGIBLE=NO`. The published benchmark protects
Somatic performance but is not a Germline/Reference truth benchmark. The pilots
provided provisional Germline support and **zero supported Reference sites**
under the approved 1%/95% per-sample rule. Negative indels still require
haplotype-aware validation. No cohort-wide BAM evidence is collected or inferred
from those pilots by this runner. This rerun prepares candidates for that review,
not an approved three-class training set.

`--approve-pilot` acknowledges a reviewed execution pilot; it does **not** grant
biological training approval. FILTER-only model consumers must not be pointed at
these candidate VCFs as if they were approved truth.

## Paths and resources

Config: `examples/seq2neo/config/separated_three_class_v2_cohort.json`.
Code may run from the current repo or its rsynced shared copy. In both cases:

- Outputs: `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/`
- Work: `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/nf_work/three_class_v2_cohort_20260919/`
- Logs: output root's `logs/prepare.*.log`, `logs/pilot.*.log`, or `logs/execute.*.log`.

The wrapper prints the exact log path before redirecting stdout/stderr.
After the completed two-worker pilot, the execution configuration uses **three
workers**, each with a 16 GiB address-space limit, with 8 GiB reserved in the
cgroup resource check (56 GiB budget total). The old two-class config remains
unchanged. This scheduling-only change preserves the completed pilot identity;
validated policy code, output paths, source hashes and approval gates do not
change. It does not increase a single sample's speed, and a 1.5x throughput gain
is only an ideal upper bound before I/O contention and scheduling imbalance.
The larger disk estimate accounts for separate baselines,
native candidates and disk-backed unions; filesystem free space is not a quota
guarantee. Outputs and work are disjoint from all original sources.

## Run sequence

From either repo root, prepare without generating sample results:

```bash
bash examples/seq2neo/run_three_class_cohort.sh --prepare
```

For the shared copy, first use:

```bash
cd /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar
```

Then generate the three configured execution pilots:

```bash
bash examples/seq2neo/run_three_class_cohort.sh --pilot
```

After reviewing those new-policy pilots, run the cohort:

```bash
bash examples/seq2neo/run_three_class_cohort.sh --execute --approve-pilot
```

Old-policy pilots do not count. Completed candidates are reused only with the
same identity, source hashes and verified output hashes. Changed code, manifest,
reference or validation gate requires a new output namespace; failed attempts
are retained. The global lock prevents two cohort drivers sharing the namespace.

## Products and export

Within `SAMPLE/attemptNNN/`:

| File | Meaning |
| --- | --- |
| `refined.vcf.gz`, `refined.rescue.vcf.gz` | Established Somatic baselines; not the new final three-class outputs |
| `native.vcf.gz` | Independent native nomination intermediate; not the final consensus |
| `three_class.consensus.vcf.gz` | Final separated DNA candidate consensus |
| `three_class.rescue.vcf.gz` | Final separated realignment-rescue candidates |
| `consensus.audit.json`, `rescue.audit.json` | Hash-bound structural checks on the final candidates |
| `completion.json` | Commands, source/output hashes, identity and final artifact paths |

The new `candidate_manifest.tsv` points to the final artifacts and leaves both
`truth_vcf` and `training_label_vcf` empty. It inventories completed samples only;
a three-sample pilot cannot imply completion of all 66.

After all 66 finish, export all three candidate classes (not only Somatic):

```bash
.venv/bin/python examples/seq2neo/scripts/build_refined_rerun_artifacts.py \
  --source-manifest /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/data/processed/sample_manifest.tsv \
  --prior-summary examples/seq2neo/data/processed/refined_native_v2_20260917/summary.json \
  --output-root /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919 \
  --outdir /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/exports \
  --parquet /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/variant_parquet/three_class_candidates.parquet
```

The exporter checks final artifact/audit hashes, refuses a fallback to old
baselines, preserves per-caller evidence and both classification traces, and
keeps `training_eligible=false`. Use fresh export destinations. The manifest
preserves original sample/input paths and leaves training-label approval empty.
