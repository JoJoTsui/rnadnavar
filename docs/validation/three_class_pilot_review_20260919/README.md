# Three-class cohort pilot review and execution tuning

The three candidate-generation pilots finished on 2026-09-19 at 08:24 Asia/Shanghai.
They ran from the working repo and published to the shared repo. No full
66-sample execution, mapping, variant calling, fresh VEP, or biological training
approval was performed by this review.

## Results and provenance

Canonical output root:
`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919`.
Each sample uses `attempt001/three_class.consensus.vcf.gz` and
`attempt001/three_class.rescue.vcf.gz`; the latter uses realignment RNA evidence.
Source VCFs, previous outputs and work caches remain unmodified.

| Sample | DNA Somatic | Rescue Somatic | Germline, each stage | Reference, each stage | DNA / rescue total records |
| --- | ---: | ---: | ---: | ---: | ---: |
| PRJNA298376_4255 | 541 | 545 | 25,656 | 86,440 | 1,115,754 / 1,313,569 |
| PRJNA298376_4278 | 523 | 528 | 26,137 | 136,125 | 1,536,518 / 1,883,593 |
| PRJNA298330_4032 | 397 | 397 | 19,194 | 68,009 | 407,917 / 719,240 |

The remaining records are NoConsensus, not training negatives. The production
reports give zero Somatic membership mismatches; all six structural audits have
no issues. The independent reviewer checks **every exact full allele and all
classes**, including missing/duplicate output records, native negative
provenance, FILTER/UNIFIED_FILTER, baseline/native class traces, and eligibility.
It uses a disk-backed union rather than an in-memory set of millions of records.
All **6,976,591 final records across six VCFs passed**. The
[machine-readable summary](summary.json) records the final verdict, detailed
timings and file hashes. Source/output hashes were rechecked after the replay.

The Somatic allele sets also match the previous cohort's corresponding
`output_refined_native_v2_20260916/SAMPLE/attempt001/refined[.rescue].vcf.gz`
exactly. Sorted `CHROM,POS,REF,ALT` SHA-256 checks agree on both sides:

| Sample | Stage | Somatic allele-set SHA-256 |
| --- | --- | --- |
| 4255 | DNA | f0a4f93f62db130e0344187b04dc6ff74e828d6c1d693c6afc6297025cf722a3 |
| 4255 | Rescue | 1204c3d5cf50c6be6399b3e3a434afd16b04d95621202185785c2c51aa6ad2b2 |
| 4278 | DNA | a9b6957ced979ea31c978859922470625608a7c81d435050231a21bf910daee8 |
| 4278 | Rescue | eaef53d41ba52a22200c9e772f2ad2f5e47f79aa57807172aa53283b8267efa1 |
| 4032 | Both | 3a800553050df0f6fbb9996006ddea80eb756629a54bb500392b132f1df3dd16 |

This is membership preservation, not a fresh truth benchmark or proof of
Germline/Reference precision. All candidates retain `TRAINING_ELIGIBLE=NO`.
The cohort runner does not collect paired-BAM negative-evidence validation;
the prior benchmark pilots' zero supported Reference yield and unresolved
negative indels remain limitations. Do not point FILTER-only training consumers
at these candidates as approved truth.

## Why the new pilot took longer

| Sample | Previous sample runtime (minutes) | Three-class sample runtime (minutes) |
| --- | ---: | ---: |
| 4255 | 25.83 | 67.05 |
| 4278 | 33.22 | 86.60 |
| 4032 | 9.53 | 24.80 |

The observed total pilot span was about 92 minutes with two workers; the third
sample waited for a worker. The earlier 09:10–10:40 completion estimate was too
conservative: the third sample was substantially smaller.

The new path intentionally adds a complete native negative-nomination pass
and two disk-backed full-allele assemblies. Nomination took approximately
6–24 minutes/sample, and assemblies together took about 7–21 minutes/sample.
The pilot samples' total runtime is consistently about 2.6 times their previous
Somatic-only path. These are different workloads; a slowdown is not evidence
of a stalled job. Individual stage estimates come from local artifact mtimes,
not a CPU profiler, and include nearby indexing/publication overhead.

## Applied optimization

Only the new three-class cohort config changes from **two to three workers**.
The validated production driver, classification/rescue/assembly code, command
arguments, source and output namespaces, and validation summary are unchanged.
Consequently the completed pilot identity is preserved and verified pilot
outputs remain reusable. The previous cohort config and Nextflow configs are
untouched.

Three workers at 16 GiB each plus an 8 GiB reserve budget 56 GiB, below the
observed 78 GiB cgroup limit and 46-CPU quota. Cgroup usage was approximately
13 GiB during review (including cache/other activity); three workers retain more
headroom than a four-worker/72-GiB budget. Existing resource checks still reject
over-budget configurations; this is not an aggregate hard memory reservation.

The theoretical throughput ceiling is 1.5x relative to two workers, **not a
measured speedup**, and individual samples do not become faster. I/O contention
and load imbalance reduce this gain. Do not assume the expanded three-class
66-sample run will finish in the previous workload's 24–26 hours.

No policy-code micro-optimization is promoted here. Sharing parsing between the
two consensus passes or replacing per-record SQLite operations requires a
separate profiler experiment and full-output differential validation; editing
the frozen driver/policy would invalidate the current pilot identity.
Reusing the previous cohort's intermediate baseline also needs
source/reference/policy provenance equivalence, not just matching filenames or
Somatic counts on three samples. No checks are skipped to make this run faster.

### Isolated assembly compression experiment (not promoted)

The existing adapter was profiled on 4278 realignment-rescue/native records in
`chr22:20000000-21000000`, yielding 3,276 union records. This is a fixed interval,
not a truth-selected variant subset. Of 2.021 profiled seconds, per-record zlib
compression used 0.525 seconds (26%); 18,174 SQLite execute calls used 0.114
seconds (6%). SQLite call-count optimization is not the largest hotspot in this
small interval.

Three alternating unprofiled repeats compared default internal compression
level 6 with level 1, patched **only in the isolated experiment process**:

| Internal SQLite compression | Seconds, three repeats | Median | SQLite size |
| --- | --- | ---: | ---: |
| Level 6 (current production) | 1.538, 1.634, 1.602 | 1.602 s | 11,661,312 bytes |
| Level 1 (experimental) | 1.238, 1.502, 1.305 | 1.305 s | 11,784,192 bytes |

All six decompressed VCFs (including headers and INFO) have the same SHA-256:
`8df48b33ecec139b69095a64d884f49853655b2b3e6ad9fe0f55132dbf871965`.
Counts are 3,100 NoConsensus, 150 Reference, 25 Germline, 1 Somatic. Level 1
reduced median assembly wall time by 18.5% on this small local fixture, with
about 1.1% extra SQLite space. It did not change final BGZF compression.

This is **not a measured whole-cohort gain or a full-data equivalence result**.
Assembly is only part of the new runtime, so the implied overall gain is much
smaller. The immediately usable scheduling change preserves the validated
pilot cache; changing the assembly implementation would require refreshed
full-data equivalence evidence and a new validated identity before promotion.
Do not use experimental output as a cohort result or validation-gate input.

Ignored experimental files are under
`examples/seqc2/comparison/three_class_pilot_review_20260919/profile_experiment/`.
They were moved intact from `/tmp/three-class-assembly-profile.X1UNsgj5/` after
execution; paths inside their temporary experiment reports retain that old
prefix. No experiment modified production Python source or original inputs.

To reproduce, use a fresh temporary directory, subset the two indexed pilot
VCFs with `bcftools view -r chr22:20000000-21000000 -Oz -o ...`, then run the
adapter with their actual SHA-256 values under `python -m cProfile -o ...`.
For the unprofiled comparison, the isolated-process override was:

```python
import apply_three_class_labels as m
original = m.zlib.compress
# Run levels in order 6,1,1,6,6,1; each run uses a fresh output directory.
m.zlib.compress = lambda value, level=level: original(value, level)
report = m.run(baseline, native, outdir,
               m.digest(baseline), m.digest(native), 'realignment')
```

Hash decompressed VCF bytes, compare counts, and include the larger SQLite
footprint when evaluating this separately from the current release.

## Reproduction and next action

Independent review, using a **new** scratch/report destination each time:

```bash
.venv/bin/python examples/seq2neo/scripts/review_three_class_pilot.py \
  --config examples/seq2neo/config/separated_three_class_v2_cohort.json \
  --outdir examples/seqc2/comparison/three_class_pilot_review_20260919
```

The review directory contains ignored SQLite scratch and `review.json`. It is
disjoint from the shared cohort output, which is read-only during review.
Lightweight summary/provenance live with this document; VCFs and SQLite do not.

Somatic membership reproduction, for each previous/current stage pair:

```bash
set -o pipefail
bcftools query -i 'FILTER="Somatic"' -f '%CHROM\t%POS\t%REF\t%ALT\n' "$VCF" | LC_ALL=C sort | sha256sum
```

Read-only 66-sample preflight with the new resource setting:

```bash
.venv/bin/python examples/seq2neo/scripts/run_refined_cohort.py \
  --config examples/seq2neo/config/separated_three_class_v2_cohort.json \
  --plan examples/seqc2/comparison/three_class_pilot_review_20260919/three_worker_preflight.json
```

This preflight passed for all 66 samples with three workers and sufficient
reported disk headroom. Existing older-than-VCF index warnings remain; these
checks establish header/first-record readability, not complete index accuracy
or an independent filesystem-quota guarantee. No original index was modified.

Verification: **87 focused tests passed**, covering independent review failures,
synthetic execution/cache reuse after changing worker count, frozen policy code,
cgroup-budget rejection, separated labels and candidate export:

```bash
.venv/bin/python -m pytest tests/seqc2/test_three_class_pilot_review.py \
  tests/seqc2/test_refined_cohort_preparation.py \
  tests/seqc2/test_refined_rerun_artifacts.py \
  tests/seqc2/test_separated_three_class.py -q
```

After execution-pilot review, the operator may launch the candidate cohort:

```bash
bash examples/seq2neo/run_three_class_cohort.sh --execute --approve-pilot
```

If launching from the rsynced repo, sync the updated configuration first.
The driver rechecks resource limits, validation identity and source/output
hashes; it should reuse the three completed pilots and generate the remaining
63 samples. `--approve-pilot` is **execution approval only**, never biological
training-label approval.
