# v2 training sample exclusions

Decision confirmed by the user on 2026-09-21: retain the following historical
training exclusions for the `separated_three_class_v2` handoff.

| Sample ID | Historical concern | Current v2 status | Training disposition |
| --- | --- | --- | --- |
| PRJNA298330_4032 | TruthQC-suspect; historical RNA-only common-AF leakage | Completed; 397 Somatic, 19,194 Germline, 68,009 Reference | Excluded |
| PRJNA298376_4081 | TruthQC-suspect, RefCall-heavy inputs; historical label anomalies | Completed; 286 Somatic, 18,241 Germline, 1,298,064 Reference | Excluded |
| PRJNA298376_4255 | TruthQC-suspect, RefCall-heavy inputs; historical label anomalies | Completed; 545 Somatic, 25,656 Germline, 86,440 Reference | Excluded |

The historical investigation is recorded in
[RERUN_LABEL_ANOMALY_ANALYSIS.md](../../RERUN_LABEL_ANOMALY_ANALYSIS.md).
It distinguishes old baseline/version effects from unresolved sample concerns;
do not assume every historical failure signature still occurs in v2.
Neither successful v2 completion nor bounded source-caller evidence checks
constitute re-admission. The 66-sample export intentionally remains a complete
candidate inventory, **not a training inclusion manifest**.

## Enforced handoff boundary

The derived `handoff_review_20260921/cohort63.review.parquet` and
`cohort63.review.tsv` exclude all records from these three samples. The direct
three-sample-stage files also contain none of them. The shared root is:

`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919`.

Do not use the original 66-sample Parquet or fall back to an unfiltered source
VCF as a replacement for the selected handoff. The exclusions concern training
and evaluation admission; no original result is deleted or relabeled.

The next training stage comprises **PRJNA298376_4007, PRJNA298376_4060 and
PRJNA298376_4072 together**. Only the separate 4007-only stage is skipped.
Existing reserved pools remain reserved. Sample inclusion does not approve
every variant: the handoff remains review-only pending its scoped release
decision and model loader/cache/split verification.

Any future re-admission requires a fresh, version-specific review and an
explicit recorded decision; a caller PASS or pipeline completion alone is not
sufficient.
