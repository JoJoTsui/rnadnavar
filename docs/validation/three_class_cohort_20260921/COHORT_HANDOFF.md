# Full-cohort v2 handoff

## Approved release

The user approved this exact weak-label selection on 2026-09-21.
[RELEASE_APPROVAL.json](RELEASE_APPROVAL.json) is the authoritative, hash-bound
approval sidecar: 58 samples for training and 5 reserved for held-out evaluation.
It includes the 3,464 Somatic exclusions and native-evidence negative labels.
This is not independent biological truth certification. Historical review files
retain false/NO flags and their hashes; consumers must explicitly acknowledge
this release sidecar rather than globally bypass eligibility checks.

Training will use another user's code on another host, not the inspected local
EvoSomatic copy. The full cohort supersedes the three-sample execution plan.

## Selection

Use the reviewed 63-sample `separated_three_class_v2` subset, excluding
PRJNA298330_4032, PRJNA298376_4081 and PRJNA298376_4255 entirely.
The approved selection also quarantines 3,464 Somatic conflicts (3,462 common-AF
annotations and two native-class conflicts); these are not proven false
positives. Native Germline/Reference labels are retained without a new
299-read requirement. Approval remains separate from technical verification.

Expected counts: Somatic 28,241; Germline 1,522,528; Reference 12,791,828;
total 14,342,597. These are workflow-derived supervision, not independent
genotype truth. Historical review flags remain false/NO; release approval is
recorded separately above and does not approve the original 66-sample inventory.

The split-table join confirms 58 train-pool and 5 reserved samples. Reserved:
PRJNA298310_3812, PRJNA298376_4166, PRJNA298376_4214,
PRJNA298376_4231 and PRJNA298376_4242. These must not be silently included
in the training pool. All 189 DN/DT/RT paths and their indexes were found;
availability does not establish BAM integrity or sample identity by itself.

## Artifacts

Root: `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/`.

- `handoff_review_20260921/cohort63.review.parquet`: annotation-rich variants.
- `handoff_review_20260921/cohort63.review.tsv`: manifest and alignment paths.
- `handoff_review_20260921/audit.json`: selection provenance and hashes.
- `cohort63_model_review_20260921/`: indexed minimal review VCFs, JSON manifests
  and report. Require a successful report; directory existence is insufficient.

Completed export on 2026-09-21: status
`review_bridge_verified_not_training_approved`, 63 VCFs and 63 tabix indexes,
14,342,597 records, 58 train-pool / 5 reserved, 189 alignment entries.
All ordered allele/class hashes matched the reviewed Parquet; source checksums
remained unchanged. This status is technical verification, not training approval.

Verification report SHA-256:
`eb6f11ef47e8c7c2e9c32853de2843894c98fa68d640f6dbb8d4b96a3cb554b8`.

Reviewed Parquet SHA-256:
`8e61c3c73b51cbca095b0da26e34f76b5adfe80c594723acaaedf4d740662a6d`.

Reviewed TSV SHA-256:
`05602cf063e96e15d9c4c9485a854446c79d0d8d39c8fc68e9e7bf3664ac4cdc`.

Split TSV SHA-256:
`08776aa101ee79056d4c12869d3d5d4d8f8398c1404af6232e811225acc200df`.

The bridge verifies full ordered CHROM/POS/REF/ALT/FILTER roundtrips, class
counts, sample membership and unchanged source hashes. Alignment/index
inventory verifies availability and reports BAM sizes, not BAM content hashes.
No original outputs, mapping, calling, classification or training are changed.

`samples.review.json` inventories all samples, not just training samples.
`train_pool.review.json` and `reserved.review.json` preserve existing pools.
Only sample IDs/pools are imported from the old split TSV; its old VCF paths,
QC approvals and counts are not reused.

## Destination checks

The destination host uses the same shared paths (confirmed by the user).
No data transfer or remapping is required. Keep the existing absolute paths;
the user accepted the completed shared-data checks; repeating them on the
destination is not required. The checker remains available for troubleshooting.

1. Use the shared selected Parquet, manifest, audit, VCF/index files and report.
   Label SHA-256 checks have passed; retain their report and release sidecar.
2. Resolve DN/DT/RT alignment/index paths using full sample IDs and matching
   reference/contig names. Supply the matching FASTA/index for the actual model.
3. Preserve reserved samples. Apply the documented chromosome split separately:
   train chr2–20, validation chr21/22, test chr1. Do not silently add chrX/Y/M.
4. Test the actual loader: Reference=0, Germline=1, Somatic=2; VCF POS is
   one-based. Account for normalization, multiallelic and long-indel scope.
   Do not convert Artifact/NoConsensus/RNAedit or unselected sites to Reference.
5. Use the approved named release, including its Somatic quarantine.
   A loader ignoring approval fields alone does not establish authorization.
6. Use a fresh cache identity bound to label content, splits, alignment/reference
   identities, model code and feature settings, not merely filenames.

## Reproduction

Use a fresh destination; no GPU, feature cache or training is initialized:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 prlimit --as=17179869184 -- \
  .venv/bin/python examples/seq2neo/scripts/build_v2_model_review_bridge.py \
  --cohort --handoff /path/to/handoff_review_20260921 \
  --split-manifest /path/to/data/processed/sample_split.tsv \
  --outdir /path/to/fresh/cohort63_model_review
```

The exporter streams 32,768-row batches without Parquet reader threads; this
command caps virtual address space at 16 GiB. Focused bridge/handoff/export
tests passed (11 tests). Heavy artifacts remain outside Git.

## Shared-host preflight (no workflow or training run)

The standalone checker uses Python's standard library; no model environment
or local development checkout is required. Optional reproduction only (not a
required repeat on the destination host):

```bash
python3 /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/cohort63_handoff_tools_20260921/check_v2_cohort_handoff.py --bridge /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/cohort63_model_review_20260921 --expected-report-sha256 eb6f11ef47e8c7c2e9c32853de2843894c98fa68d640f6dbb8d4b96a3cb554b8
```

Success is `shared_data_preflight_passed_not_training_approved`. The command
only reads inputs and prints JSON. Source code paths in historical provenance
belong to the exporter host; they are not required on the destination. Shared
source data and all exported VCF/index/manifest hashes are rechecked.
BAMs are checked for readable nonempty files and original sizes, not fully
hashed or revalidated. The report's class counts are trusted only after its
recorded hash and exact output hashes match the previous full roundtrip audit.

Destination model code, FASTA configuration, label acceptance, chromosome
splits and cache settings still require that user's loader preflight. The local
EvoSomatic implementation is not used. Do not train on `samples.review.json`:
it includes reserved samples; use the appropriate pool with RELEASE_APPROVAL.json.
