# Training-readiness evaluation

> Scope correction after reviewing the manuscript's current integrated plan:
> the paired-read pilot below is a stringent diagnostic, **not** a mandatory
> v2 weak-supervision release gate. Reference-specific research is deferred.
> See [the v2-only handoff](V2_HANDOFF.md) for the subsequent efficient audit,
> historical sample exclusions and direct three-sample-stage preparation.

The verified Parquet and manifest are suitable for candidate exploration,
feature extraction and pipeline integration tests, not automatically approved
supervised training or truth-based performance evaluation. FILTER is a candidate
class; `training_eligible=false` remains authoritative. Consumers that ignore
eligibility must not be pointed at this dataset as training truth.

## Scope

This follow-up first audits all exported candidate class/review-reason counts,
then samples 32 Germline and 32 Reference SNPs per configured execution pilot:
`PRJNA298330_4032`, `PRJNA298376_4278`, `PRJNA298376_4255`. Selection is the
smallest SHA256 of CHROM:POS:REF:ALT separately per sample/class, independent of
truth or read evidence. These three samples are an execution pilot, not a
representative cohort accuracy study. Both native-negative/Somatic conflicts
are inspected separately and must not be mixed into pilot retention estimates.

The script reuses the existing BAM pileup implementation and strict negative
evidence assessment. MAPQ/BQ >=20, BAQ, duplicate/secondary/supplementary filtering,
overlapping-mate filtering and orphan exclusion are retained. Reference requires
zero ALT/other observations and >=299 quality-filtered observations in **each**
DNA sample at the approved idealized 1%/95% criterion. The older exploratory
60-read Reference categorization in the pileup helper is **not** used.

RNA measurements use the manifest's `bam_rt` (original mapped RNA), not an
assumed realignment BAM. They are exploratory and do not approve a DNA negative.
No coverage is not evidence of a Reference genotype. Negative indels remain
outside this SNP pileup evaluation. No labels, original VCFs, BAMs, indexes,
mapping/calling caches, or model configurations are changed.

## Reproduce

Use a fresh output directory; existing results are refused:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 prlimit --as=17179869184 -- \
  .venv/bin/python -u examples/seq2neo/scripts/evaluate_three_class_training_readiness.py \
  --config examples/seq2neo/config/separated_three_class_v2_cohort.json \
  --outdir /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/training_readiness_20260921
```

Report: `training_readiness_20260921/report.json` under the shared cohort root.
Log: `logs/training_readiness_20260921.log`. Execution uses one worker and a
16 GiB virtual-address-space cap; the observed resident memory is below 0.5 GiB.

## Input caveats

Existing indexes are older than their BAMs for the DN/DT of
`PRJNA298330_4032` and DN/DT/RT of `PRJNA298376_4112`. HTSlib warns but can
query them. This does not prove corruption; index provenance/validity remains
unresolved. Do not treat read support as certification of these inputs. Original
indexes were not rebuilt or modified. The other inspected DN/DT/RT indexes for
4252, 4255 and 4278 are newer than their BAMs; timestamp order alone is not an
index-integrity proof either.

The evaluator binds the verified export and helper code by SHA256, and checks
BAM/reference sizes and modification times. Whole BAMs and indexes are not
rehashed. This is a pilot evidence screen, not a full input-integrity audit.

## Completed pilot findings

Status: `pilot_complete_not_training_approved`. All 194 sites (192 pilot + two
targeted conflicts) were measured. Export counts and code/input integrity checks
passed. See [lightweight results](training_readiness_summary.json) for the exact
report path/hash, outcomes and conflict-site counts. All 33 focused tests passed.

| Class | Pilot sites | Provisionally supported | Withheld |
| --- | ---: | ---: | ---: |
| Germline | 96 | 52 | 44 |
| Reference | 96 | 0 | 96 |

Germline withholding: 29 other-allele contexts, eight insufficient paired
evidence, seven large allele-fraction shifts. Reference withholding: 56 with
some ALT evidence, 39 with insufficient depth, one other-allele context. These
are evidence-gate outcomes, not measured FP counts. In particular, the native
Reference nomination permits limited ALT observations, whereas this training
evidence gate requires zero; nomination confidence is not training eligibility.

Supported Germline counts per pilot are 22/32 (4032), 15/32 (4255), 15/32 (4278).
The 4032 results retain the index-provenance caveat above. Do not extrapolate
these counts into a cohort-wide approved-label yield.

The Germline/Somatic conflict at chr6 has normal ref/alt=8/1 (depth 9), tumor
6/18 (24), and RNA 2/16 (18). Normal coverage is insufficient to adjudicate it.
The Reference/Somatic conflict at chrX has normal 29/0 (29), tumor 40/1 (41),
and RNA 0/289 (289). Its strong RNA signal alone establishes neither a somatic
DNA mutation nor a reliable Reference label. Both remain withheld for training.

Among the 192 pilot sites, 106 have zero quality-filtered RNA observations;
34 have RNA depth >=20 (23 Germline, 11 Reference). Twenty-one Germline sites
have >=3 RNA ALT observations. This shows uneven RNA evidence availability,
not a measured benefit or harm to a model. The three-class dataset is **not**
approved for training: no Reference pilot cleared the selected criterion,
negative indels are not validated, and Somatic conflicts/QC remain unresolved.

## Historical proposed gates (superseded as release prerequisites)

The following were proposed before checking the updated manuscript plan.
They remain possible research/quality investigations, not requirements that all
Reference labels meet a 299-read criterion before workflow-derived training.

1. Adjudicate Somatic annotation/native-class conflicts and complete label QC.
2. Extend negative read-evidence validation beyond the pilot; validate indels
   with an appropriate haplotype-aware method. Resolve index provenance.
3. Make an explicit biological label approval decision; export an approved
   subset separately, never relabel missing/withheld evidence as Reference.
4. Only then evaluate DNA-only versus DNA+RNA models on identical labels and
   patient-disjoint splits, stratified by RNA coverage. Model training has not
   been launched or claimed by this evaluation.
