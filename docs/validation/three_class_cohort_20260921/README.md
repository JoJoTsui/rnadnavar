# Completed three-class candidate cohort and export review

Policy: `separated_three_class_v2`. All 66 samples completed; no failed sample.
Full execution began 2026-09-19 13:52:06 and the last sample finished
2026-09-21 01:46:15 Asia/Shanghai (about 35 h 54 min). Three completed pilots
were reused. No mapping, variant calling or fresh VEP annotation was run.

## Canonical paths

Shared output root:
`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919`.

- Final DNA candidates: `SAMPLE/attempt001/three_class.consensus.vcf.gz`.
- Final realignment-rescue candidates: `SAMPLE/attempt001/three_class.rescue.vcf.gz`.
- Completion inventory: `candidate_manifest.tsv` (66 samples).
- Export manifest: `exports/manifest.tsv`.
- Export summary: `exports/summary.json`.
- Variant Parquet: `variant_parquet/three_class_candidates.parquet`.
- Sample-level Parquet: `variant_parquet/refined_three_class_manifest.parquet`.

The final VCF is the **realignment** rescue, not first rescue. Original caller
VCFs, BAMs/CRAMs, reference files, work caches and previous cohort outputs are
read-only. Heavy files are stored in the shared output tree, not in Git.

## Audit review

All 132 structural audits report no issues. All 132 separated-label reports
declare zero Somatic membership mismatches and unchanged inputs. Their JSON
hashes match completion records, final VCF hashes agree between state/report/
audit, and all 66 completion records match the frozen policy identity and
current production code. The candidate manifest has exactly the expected 66
unique sample IDs and empty truth/training-label approval fields.

Counts match the previous refined rerun for **Somatic in every sample and both
stages**. The full export independently rehashes final VCFs and audits and checks
each exported class count against its rescue audit. A separate full-allele
Somatic comparison checks equality of variant membership, not just counts.
Export and verification are complete: all 66 samples and all 15,862,893 Parquet
rows passed. All 132 old/new Somatic stage comparisons have zero added/removed
alleles. These are integrity and membership checks, not biological approval.

| Stage | Somatic | Germline | Reference | NoConsensus |
| --- | ---: | ---: | ---: | ---: |
| Previous refined DNA | 32,725 | 15,839 | 38,480 | 105,963,208 |
| New separated DNA | 32,725 | 1,585,620 | 14,244,342 | 90,264,266 |
| Previous refined rescue | 32,933 | 986,181 | 7,197,503 | 114,226,688 |
| New separated rescue | 32,933 | 1,585,619 | 14,244,341 | 109,339,082 |

The previous DNA also contained 76,701 Artifact records; previous rescue had
2,090,509 Artifact and 668,161 RNAedit records. New separated outputs contain
only Somatic/Germline/Reference/NoConsensus by the declared union contract.
Previous counts are the recorded prior summary, not freshly rescanned counts
for all previous non-Somatic classes. New counts come from completion audits
and are cross-checked by the exporter.

The candidate Parquet includes all three selected classes from final rescue:
**15,862,893 rows verified**, not just the Somatic SNP/indel subset.
NoConsensus is deliberately excluded; it must not be interpreted as Reference.
Full ALT fields, including multiallelic records, are retained.

| Exported class | SNP | Indel |
| --- | ---: | ---: |
| Somatic | 31,571 | 1,362 |
| Germline | 1,479,633 | 105,986 |
| Reference | 14,244,341 | 0 |

The two negative-to-Somatic conflicts were checked in both final VCFs:
`PRJNA298376_4112 chr6:158621604 G>A` (Germline) and
`PRJNA298376_4252 chrX:13379221 T>G` (Reference). Both pass the established
`dna_nominated_rna_supported` rescue branch and retain
`THREE_CLASS_REVIEW_REASON=native_negative_conflict`. Somatic precedence
preserves the prior result; it does not establish the correct biological class.

The exported review reasons include 3,607 Somatic common-population-AF conflicts,
two native-negative conflicts, 29,324 remaining Somatic candidates needing label
QC, and 15,829,960 negative candidates needing paired evidence. No row is approved.

There are no sample-level increases in Germline or Reference between the new
DNA consensus and new rescue. The two negative classes each decrease by one
record in the cohort totals, consistent with stage-specific Somatic precedence.
Rescue adds 208 Somatic records overall relative to DNA consensus; that is the
same Somatic behavior as the established refined baseline.

## Interpretation limits

This validates candidate generation and preservation of the established Somatic
results, **not three-class biological accuracy**. Every exported row must retain
`training_eligible=false`, `TRAINING_ELIGIBLE=NO`, and the separated policy ID.
The exported manifest is marked `NOT_APPROVED`; no training-label VCF is set.

The larger Germline/Reference counts reflect independent native DeepSomatic
nominations and the separate negative-class rules, not newly established truth.
Reference candidates are caller-emitted candidate sites, not all genomic
reference bases. Cohort-wide paired-BAM evidence was not collected by this run.
The approved 299-observation-per-DNA-sample Reference evidence rule was not
silently relaxed or automatically granted to native candidates. Negative indels
still require haplotype-aware assessment. Do not use this Parquet as approved
training truth without a separate label-QC/evidence/biological approval step.

## Reproduction

From the current repo or a synchronized copy, export to **fresh** destinations:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 prlimit --as=17179869184 -- \
  .venv/bin/python -u examples/seq2neo/scripts/build_refined_rerun_artifacts.py \
  --source-manifest /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/data/processed/sample_manifest.tsv \
  --prior-summary examples/seq2neo/data/processed/refined_native_v2_20260917/summary.json \
  --output-root /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919 \
  --outdir /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/exports \
  --parquet /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/variant_parquet/three_class_candidates.parquet
```

The exporter refuses existing destinations. It uses bounded 50,000-row batches,
not an in-memory concatenation of the cohort. Export log:
`logs/export.XEdaHPK7.log` under the shared output root.

Exact Somatic membership comparison:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 prlimit --as=17179869184 -- \
  .venv/bin/python -u examples/seq2neo/scripts/compare_three_class_somatic_sets.py \
  --previous-root /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_refined_native_v2_20260916 \
  --current-root /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919 \
  --manifest /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/data/processed/sample_manifest.tsv \
  --report /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/review_20260921/somatic_membership.retry16g.json
```

This compares full `CHROM,POS,REF,ALT` sets for 132 old/new pairs. It records
declared prior VCF hashes but does not independently rehash all prior VCF bytes;
current output hashes are verified by the exporter. Neither check proves truth
accuracy. Log: `logs/somatic_comparison.9pQz2jom.log`.

An initial verifier-only attempt with a 4 GiB virtual-address-space limit hung
before opening its first VCF. The environment reproduced this with
`prlimit --as=4294967296 -- bcftools --version` (timeout), while the same command
at 16 GiB succeeded. Only that verification child was terminated; no cohort or
exporter was stopped. Its failed `review_20260921/somatic_membership.json` and
`logs/somatic_comparison.9wvyFTqv.log` are retained as diagnostics, **not canonical
validation results**. Use the `retry16g` report. The exact linked-library cause
was not established; this is a startup/address-space limitation, not evidence
of VCF corruption or a need for 16 GiB resident RAM.

After successful export, verify every Parquet row with bounded Arrow batches:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 prlimit --as=17179869184 -- \
  .venv/bin/python examples/seq2neo/scripts/verify_three_class_export.py \
  --summary /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/exports/summary.json \
  --manifest /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/exports/manifest.tsv \
  --report /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/review_20260921/parquet_verification.json
```

Verification scripts also refuse existing report paths. Tests:

```bash
.venv/bin/python -m pytest tests/seqc2/test_three_class_cohort_comparison.py \
  tests/seqc2/test_three_class_export_verification.py \
  tests/seqc2/test_refined_rerun_artifacts.py -q
```

All 16 focused tests passed, including equal-count/different-allele failure,
duplicate Somatic alleles, wrong policy, false training approval, count/hash
mismatches, and all-three-class export behavior.
