# Full-input validation and retained-conflict review

Continuation after checkpoint `c4e92fdb`, authorized by the user. No policy
threshold, default configuration, cohort label or caller cache is changed.
HG008 results remain evaluation outcomes rather than threshold-tuning inputs.

## Full-input SEQC2 evaluation

The same `validate_frozen_hybrid_policy.py` used for HG008 was run against both
verified SEQC2 caller bundles. The candidate selection sees all supplied DNA
caller records; HC and target BEDs restrict scoring only. This checks whether
the earlier region-scoped baseline concealed candidate-selection or protected
negative-label differences. All first/rescue RNA panels and annotated sources
match the previous SEQC2 assays. WES first and realignment outputs retain their
distinct historical provenance; this is not a new complete hybrid execution.

Manifests:

- `examples/seqc2/hybrid/frozen_seqc2_wes_validation.json`
- `examples/seqc2/hybrid/frozen_seqc2_wgs_validation.json`

These are exact-path evaluation manifests, not relocatable workflow defaults.
Both passed the 16-input preflight. Commands from the repository root:

```bash
.venv/bin/python examples/seqc2/scripts/validate_frozen_hybrid_policy.py \
  --manifest examples/seqc2/hybrid/frozen_seqc2_wes_validation.json \
  --outdir examples/seqc2/comparison/current_native_audit_20260916_v2/seqc2_wes_frozen_full_domain_v1

.venv/bin/python examples/seqc2/scripts/validate_frozen_hybrid_policy.py \
  --manifest examples/seqc2/hybrid/frozen_seqc2_wgs_validation.json \
  --outdir examples/seqc2/comparison/current_native_audit_20260916_v2/seqc2_wgs_frozen_full_domain_v1
```

Both evaluations completed: 24 benchmark comparisons in total, with source,
code and manifest integrity checks passing. Each evaluates consensus,
DeepSomatic, and workflow/optimized first and realignment rescues on UKB and
MedExome. All eight optimized-rescue comparisons reproduce the region-scoped
assays exactly for SNP, indel and aggregate metrics. Full-input consensus also
reproduces the earlier metrics in all four cells. No old result was replaced
or deleted. Machine-readable metrics and all commands are in each output's
`validation.json`.

| Dataset/target | DeepSomatic F1 | Optimized consensus F1 | First rescue F1 | Realignment rescue F1 |
|---|---:|---:|---:|---:|
| WES/UKB | 0.619019 | 0.621199 | 0.624301 | 0.625294 |
| WES/MedExome | 0.798016 | 0.799433 | 0.802539 | 0.803383 |
| WGS/UKB | 0.966347 | 0.968346 | 0.968346 | 0.968346 |
| WGS/MedExome | 0.913468 | 0.916070 | 0.916070 | 0.916070 |

These are aggregate record F1 values, not evidence of superiority in every
variant type: WES MedExome optimized indel F1 remains 0.750000 versus
DeepSomatic's 0.779661. The full-input check supports the earlier SEQC2 aggregate
gains without truth/target preselection of DNA candidates. It does not remove
the HG008 MedExome regression or establish a universal rescue benefit.

Independent structural audits of all six new full-input output VCFs completed
with zero detected structural issues and unchanged source checksums. They are
recorded separately as
`consensus.structural_audit_v1.json`, `first.structural_audit_v1.json`, and
`realignment.structural_audit_v1.json` in each dataset's output directory.

| Dataset | Output | Records audited | Structural issues |
|---|---|---:|---:|
| WES | Consensus | 188716 | 0 |
| WES | First rescue | 999099 | 0 |
| WES | Realignment rescue | 441973 | 0 |
| WGS | Consensus | 466878 | 0 |
| WGS | First rescue | 1213602 | 0 |
| WGS | Realignment rescue | 660799 | 0 |

All launched evaluations and structural audits in this follow-up are complete.
Structural passes do not certify biological label correctness or authorize
workflow-default promotion or cohort training-label regeneration.

Metric parity is not whole-VCF identity: full-input WES consensus contains
2,504 Somatic records, and each optimized rescue contains 2,518 (14 additional
records overall). Full-input WGS consensus contains 6,074 Somatic records;
first rescue contains 6,076 and realignment rescue 6,079.
These whole-output counts differ from the earlier region-scoped baselines;
the protected-negative and baseline-retention rules now see the complete
supplied DNA candidate universe. Identical HC/target scoring does not certify
labels outside those scoring regions.

## HG008 retained-baseline annotation conflicts

`audit_frozen_baseline_conflicts.py` reviews both rounds without changing
labels. It reconstructs flagged baseline alleles from read-only union SQLite
files and requires their counts to match the adapter reports. It replays the
frozen full consensus query with retained som.py scratch partitions, requiring
exact SNP/indel/aggregate metric parity before assigning any partition labels.
Unmatched alleles are marked `unscored_or_representation_unresolved`, not FP.

The audit compares exact gnomAD exomes v4.1 allele AFs and retains database
record filters. Caller tumor/normal FORMAT observations are recorded per caller;
missing measurements are not filled with zero. Missing Strelka indel AD is
unavailable in this helper, not inferred from SNP fields. AF confirmation does
not establish a somatic/germline diagnosis or approve the database record.

```bash
.venv/bin/python examples/seqc2/scripts/audit_frozen_baseline_conflicts.py \
  --validation examples/seqc2/comparison/current_native_audit_20260916_v2/hg008_frozen_full_domain_v1/validation.json \
  --outdir examples/seqc2/comparison/current_native_audit_20260916_v2/hg008_baseline_conflicts_v1
```

Reports retain input hashes and large-file size/mtime checks; database record
hashes identify matched observations without hashing the full gnomAD database.
Original indexes are not regenerated. This audit was launched after its helper
tests passed. Its result remains `not_training_approval`, even when read-only
integrity and annotation confirmation pass. Any unresolved biological conflicts
are documented for release decisions; they are not silently corrected or used
to tune the HG008 policy.

## Completed HG008 conflict evidence review

The new audit completed with input hashes and large-file stat checks unchanged.
It reconstructed exactly the same 48 baseline conflicts in both HG008 rescue
rounds, all flagged for population AF. All 48 source AF values match exact
gnomAD exomes v4.1 alleles. Both benchmark replays passed exact metric parity.

| Scoring target | TP | FP | Unscored or representation-unresolved |
|---|---:|---:|---:|
| UKB | 34 | 6 | 8 |
| MedExome | 2 | 1 | 45 |

These are partition labels, not 48 false positives. Of the UKB TP, five are
SNPs and 29 are indels/other alleles; all six FP are SNPs. Mutect2 normal
alternate counts are available for 30 TP observations, including five with
positive normal alternate counts, and two FP observations, including one
positive. Neither population AF nor a blanket nonzero-normal-alt rule is a
clean separator in this selected subset. Missing evidence stays unavailable;
these observations do not establish the biological cause of each conflict.

Several local gnomAD indexes have older modification times than their VCFs.
Indexed lookups succeeded and exact matching alleles were recovered; no original
index was regenerated. Database FILTER values are retained in the audit, and
confirmation of a matching AF is not endorsement of a database record's quality.

Decision: retain the frozen labels and thresholds. The AF evidence review is
complete and finds no AF lookup mismatch for these sites, but broader biological
training-label approval remains unresolved. Do not silently veto baseline
records or tune a rule against these HG008 labels. Full evidence is in
`hg008_baseline_conflicts_v1/audit.json` under the comparison root above.

Combined regression verification passed 273 tests with 62 existing NumPy
deprecation warnings. No policy or label changes were made in response to these
findings, and no original inputs, indexes or caches were modified.
