# HG008 frozen-policy benchmark results — 2026-09-16

Source: `examples/seqc2/comparison/current_native_audit_20260916_v2/hg008_frozen_full_domain_v1/validation.json`.
All 12 scoring runs completed; input, code and manifest integrity checks passed.
Rules remain frozen. This is evaluation, not HG008-guided threshold selection.

## UKB

| Method | Type | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| DNA DeepSomatic | snp | 408 | 15 | 22 | 0.964539 | 0.948837 | 0.956624 |
| DNA DeepSomatic | indel | 64 | 9 | 179 | 0.876712 | 0.263374 | 0.405063 |
| DNA DeepSomatic | records | 472 | 24 | 197 | 0.951613 | 0.705531 | 0.810300 |
| Optimized consensus | snp | 410 | 16 | 20 | 0.962441 | 0.953488 | 0.957944 |
| Optimized consensus | indel | 90 | 4 | 153 | 0.957447 | 0.370370 | 0.534125 |
| Optimized consensus | records | 500 | 20 | 169 | 0.961538 | 0.747384 | 0.841043 |
| Optimized first rescue | snp | 410 | 21 | 20 | 0.951276 | 0.953488 | 0.952381 |
| Optimized first rescue | indel | 90 | 4 | 153 | 0.957447 | 0.370370 | 0.534125 |
| Optimized first rescue | records | 500 | 25 | 169 | 0.952381 | 0.747384 | 0.837521 |
| Optimized realignment rescue | snp | 410 | 19 | 20 | 0.955711 | 0.953488 | 0.954598 |
| Optimized realignment rescue | indel | 90 | 4 | 153 | 0.957447 | 0.370370 | 0.534125 |
| Optimized realignment rescue | records | 500 | 23 | 169 | 0.956023 | 0.747384 | 0.838926 |
| Workflow first rescue | snp | 410 | 148 | 20 | 0.734767 | 0.953488 | 0.829960 |
| Workflow first rescue | indel | 80 | 2 | 163 | 0.975610 | 0.329218 | 0.492308 |
| Workflow first rescue | records | 490 | 150 | 179 | 0.765625 | 0.732436 | 0.748663 |
| Workflow realignment rescue | snp | 410 | 111 | 20 | 0.786948 | 0.953488 | 0.862250 |
| Workflow realignment rescue | indel | 81 | 2 | 162 | 0.975904 | 0.333333 | 0.496933 |
| Workflow realignment rescue | records | 491 | 113 | 178 | 0.812914 | 0.733931 | 0.771406 |

## MedExome

| Method | Type | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| DNA DeepSomatic | snp | 135 | 4 | 19 | 0.971223 | 0.876623 | 0.921502 |
| DNA DeepSomatic | indel | 19 | 2 | 39 | 0.904762 | 0.327586 | 0.481013 |
| DNA DeepSomatic | records | 154 | 6 | 57 | 0.962500 | 0.729858 | 0.830189 |
| Optimized consensus | snp | 135 | 4 | 19 | 0.971223 | 0.876623 | 0.921502 |
| Optimized consensus | indel | 17 | 0 | 41 | 1.000000 | 0.293103 | 0.453333 |
| Optimized consensus | records | 152 | 4 | 59 | 0.974359 | 0.720379 | 0.828338 |
| Optimized first rescue | snp | 135 | 5 | 19 | 0.964286 | 0.876623 | 0.918367 |
| Optimized first rescue | indel | 17 | 0 | 41 | 1.000000 | 0.293103 | 0.453333 |
| Optimized first rescue | records | 152 | 5 | 59 | 0.968153 | 0.720379 | 0.826087 |
| Optimized realignment rescue | snp | 135 | 4 | 19 | 0.971223 | 0.876623 | 0.921502 |
| Optimized realignment rescue | indel | 17 | 0 | 41 | 1.000000 | 0.293103 | 0.453333 |
| Optimized realignment rescue | records | 152 | 4 | 59 | 0.974359 | 0.720379 | 0.828338 |
| Workflow first rescue | snp | 136 | 49 | 18 | 0.735135 | 0.883117 | 0.802360 |
| Workflow first rescue | indel | 17 | 0 | 41 | 1.000000 | 0.293103 | 0.453333 |
| Workflow first rescue | records | 153 | 49 | 58 | 0.757426 | 0.725118 | 0.740920 |
| Workflow realignment rescue | snp | 136 | 28 | 18 | 0.829268 | 0.883117 | 0.855346 |
| Workflow realignment rescue | indel | 17 | 0 | 41 | 1.000000 | 0.293103 | 0.453333 |
| Workflow realignment rescue | records | 153 | 28 | 58 | 0.845304 | 0.725118 | 0.780612 |

## Interpretation and decision

- UKB optimized consensus versus DeepSomatic: net +28 TP, −4 FP, −28 FN;
  aggregate F1 improves from 0.810300 to 0.841043. Most of the TP gain is indels
  (+26 TP, −5 FP). SNPs gain 2 TP with 1 additional FP, so SNP precision is
  slightly lower despite a higher SNP F1.
- MedExome optimized consensus versus DeepSomatic: net −2 TP, −2 FP, +2 FN.
  SNP metrics tie; the indel trade-off lowers aggregate F1 from 0.830189 to
  0.828338 despite improved precision. This is not a universal win.
- Relative to optimized consensus, UKB first rescue adds 5 FP and realignment
  rescue adds 3 FP, with no net TP/FN improvement. MedExome first rescue adds
  1 FP; realignment rescue ties consensus. No threshold is changed in response.
- Both optimized rescues greatly improve F1 over their historical workflow
  counterparts, but that does not establish a benefit over optimized DNA
  consensus. The workflow MedExome rescues have one additional TP, accompanied
  by substantially more FP.
- Current results do not support enabling the optimized rescue by default on
  the claim that it always improves consensus, or a three-dataset SOTA claim.
  Keep experimental opt-in status and separate structural correctness from
  biological quality. Further policy development belongs on SEQC2; preserve
  HG008 outcomes rather than optimizing away these failures.

These are som.py metrics from HG008 v0.3 truth, its HC BED (-R), and the named
BED (-T), using the configured GRCh38 FASTA and no -P. All supplied DNA caller
records entered consensus; truth and target intervals were used only for
scoring. This does not imply unrestricted WGS caller coverage. Aggregate
records are som.py's native aggregate, not a recomputed SNP+indel sum; truth
representation can produce different totals across those categories.

## Authoritative artifacts

### Full-output label scope and unresolved conflicts

The optimized consensus contains 724,384 records, including 604 Somatic
records. First rescue contains 2,076,392 records and realignment rescue
974,609; each retains all 604 baseline Somatic records and admits five
additional Somatic records. These whole-output admission counts are not
interchangeable with HC/target-restricted benchmark TP/FP counts (for example,
realignment adds five records overall but only three scored UKB FP).

Both rescue reports flag 48 baseline population/editing annotation conflicts.
Baseline-versus-historical-rescue-negative conflicts are 18 for first rescue
and 15 for realignment. These counters identify retained conflicting evidence;
they are not automatic FP labels, nor evidence of structural corruption.
No conflict was silently relabeled. They remain reasons to withhold a
biological training-label release even if all structural checks pass.

Consensus structural checks passed across all 724,384 records, with unchanged
source checksum. First rescue (2,076,392 records) and realignment rescue
(974,609 records) also completed with zero structural issues and unchanged
source checksums. All three audits passed the implemented structural checks;
this does not resolve the biological conflicts described above.

### Paths

All new artifacts are under the source report directory above:

- `refined.vcf.gz`: full-input optimized consensus.
- `first/refined.rescue.vcf.gz`: optimized first rescue.
- `realignment/refined.rescue.vcf.gz`: optimized realignment rescue.
- `validation.json`: exact source paths, SHA-256, commands, all metrics.
- `*.query.vcf.gz`: PASS-only benchmark copies, not training-label VCFs.
- `consensus.structural_audit_v1.json`, `first.structural_audit_v1.json`,
  `realignment.structural_audit_v1.json`: independent label-structure audits
  (all complete with `structural_pass_not_training_approval`).

The historical workflow output remains under the shared repo's
`examples/seqc2/hybrid/output/hg008.wgs.hybrid/`. The manifest and source hashes
retain exact provenance for both rescue rounds. No input or workflow cache was
rewritten, and no workflow default or cohort training label was promoted.
