# Cohort 63 approved labels and BAM evidence

This record integrates the approved 63-sample label release with its separately
generated BAM-derived variant evidence sidecar for the same downstream training
dataset. The approved release files remain immutable; this document is a
companion record rather than a modification of the hash-bound approval report.

The approved cohort report is
`examples/seq2neo/output_three_class_v2_20260919/cohort63_model_review_20260921/report.json`.
Its SHA-256 is
`eb6f11ef47e8c7c2e9c32853de2843894c98fa68d640f6dbb8d4b96a3cb554b8` and the
approval sidecar SHA-256 is
`ce21605f5bbc24f2166b36010637c6fdbbf83be7c900b60facb5b6dbf6ad6f80`.
The report describes the approved weak-label release and remains
`review_bridge_verified_not_training_approved` with `training_approved=false`;
the downstream consumer must apply its own training approval gate.

The verified evidence sidecar is
`examples/seq2neo/info_v1_20260922/`, generated at extractor revision
`3d20646a2b7c9a536a71ffa2c2d4f605ee58d238`. It contains separate
`development/` and `reserved/` partitions, retains the approved `FILTER`, and
keys every row by the full
`sample_id,CHROM,POS,REF,ALT` identity. It recomputes DN/DT/RT evidence from
the registered BAMs; it does not copy historical measurements or classification
annotations.

## Reconciled release

The approved cohort report and evidence verifier agree exactly:

| Quantity | Approved release | Verified evidence |
|---|---:|---:|
| Samples | 63 | 63 |
| Records | 14,342,597 | 14,342,597 |
| Reference | 12,791,828 | 12,791,828 |
| Germline | 1,522,528 | 1,522,528 |
| Somatic | 28,241 | 28,241 |
| Train-pool samples | 58 | 58 |
| Reserved samples | 5 | 5 |

The evidence split totals are train 11,722,395, validation 525,937, test
1,288,502, and reserved 805,763. Reserved samples remain isolated across all
chromosomes. Training-pool records use chr1 as test, chr21/chr22 as validation,
and the remaining approved training-pool records as train, as documented in
`TRAINING_DATA_GUIDE.md`.

The strict verifier completed with `full_release_complete=true`, zero remaining
approved rows, and PASS for full-key/FILTER equality, output hashes, pool and
split assignment, and numerical invariants. Export runtime was 88,550.43
seconds (~24.6 hours) with two workers; verification runtime was 401.75 seconds.
The full-run verifier did not inspect held-out model performance and did not
run model training, variant calling, alignment, or quantification.

## Reproduction and audit files

The exact commands and frozen extraction settings are in
`info_v1_20260922/RUN.md` and `examples/seq2neo/docs/VARIANT_EVIDENCE.md`.
The authoritative generated reports are:

- `info_v1_20260922/manifest.json`
- `info_v1_20260922/schema.json`
- `info_v1_20260922/qc.json`
- `info_v1_20260922/verification.json`
- `info_v1_20260922/verification.md`

At integration time their SHA-256 values were:

| File | SHA-256 |
|---|---|
| `manifest.json` | `55ce30aee1b122266b120127e79eb797800825d6af7b71d272de7acd5469dbbb` |
| `schema.json` | `808e3988d01a64a6b39f138da6a5b860246f644957af8abbe4438af183f5679c` |
| `qc.json` | `4198273141736a29ecc81ea794a44a34082c259460219c5b8ca5ad7329221ccd` |
| `verification.json` | `1caadb08084f42322923c41c4efebdc32d2db98ac9f2fb2be11a5c3b779b97ae` |
| `verification.md` | `9a19b371bca9b2d3794de30c266e5ef95f97863b079cd0b910ff3aa3c82d683e` |

The generated parquet data and per-sample receipts remain in the versioned
output directory and are intentionally not copied into the approved release.
