# Rescue-gate domain audit — 2026-09-14

Follow-up to the [nomination-gate investigation](SEQC2_RESCUE_FP_INVESTIGATION.md),
committed in `1a3b6174`. The candidate policy is frozen for this audit: no
threshold changes, candidate additions, production changes, or new workflow runs.

## Question and result

Does the apparent FP reduction hide true-variant losses outside the originally
scored targets? Partition all historical rescue additions by the same truth HC
BED and UKB/MedExome targets, then attribute each exclusion to the fixed rule.

| Dataset | Historical additions | HC additions retained | HC additions excluded | Outside-HC retained | Outside-HC excluded |
| --- | ---: | ---: | ---: | ---: | ---: |
| WES-LL | 303 | 11 TP + 2 FP | 2 FP | 64 | 224 |
| WGS-IL | 610 | 0 | 3 FP | 131 | 476 |

All HC additions are already inside UKB. There are **no additional HC sites
outside UKB** in these historical candidate sets, including no MedExome-only
HC additions. Therefore expanding this particular audit beyond the target
does not provide an independent test of the rule or expose additional scored TPs.
This is a limitation of the historical candidate universe, not a general claim
about WGS or MedExome coverage.

Outside HC, excluded records break down as follows:

- WES: 224 have `GNOMAD_AF > 0.001`.
- WGS: 474 have `GNOMAD_AF > 0.001`; two others lack DNA nomination.
- The editing veto removes no additional site beyond the nomination veto in
  these datasets. Its independent benefit has **not** been demonstrated here.

The biological veto's independent scored gain is one common-AF WES FP. Thus
the eight-run comparison supports nomination plus a population-frequency gate;
it does not establish a separate empirical gain from the editing component.

Outside-HC records remain **unassessed**, not FP. Common population frequency
provides a selection rationale but is not a substitute for truth validation.
The retained 64 WES and 131 WGS outside-HC records are likewise unvalidated;
the candidate VCF must not be presented as validated training labels globally.

## Method and safeguards

The audit reads the immutable verified native/gated VCFs, original annotated
second-rescue VCF, and saved nomination/combined exclusion lists. It:

1. Derives the full rescue-only allele set without truth-based selection.
2. Queries HC membership with bcftools `-R`; UKB/MedExome membership adds `-T`.
3. Checks exact SNP truth presence only within HC. These same-domain counts
   agree with the previously completed som.py rescue deltas; this is not a
   new indel/haplotype benchmark.
4. Reconstructs exclusion reasons and asserts exact agreement with every
   frozen selection decision, preserving overlapping reasons.
5. Records source hashes, a row per allele, and mutually exclusive domain counts.

No outside-HC record is assigned TP/FP even if a target contains it. The
existing candidate gate and historical benchmark results are not overwritten.

Evidence:

`examples/seqc2/comparison/rescue_fp_investigation_20260914/domain_audit/domain_audit.json`

Reproduce from the repo root into a fresh output directory:

```bash
.venv/bin/python examples/seqc2/scripts/audit_rescue_gate_domain.py --outdir examples/seqc2/comparison/rescue_domain_NEW
.venv/bin/python -m pytest tests/test_rescue_gate_domain_audit.py tests/test_rescue_nomination_experiment.py tests/test_historical_native_gate_replay.py -q
```

Fourteen focused tests passed. This is a fixed-bundle forensic assay, not a
portable launch script or an independent validation cohort.

## Next validation boundary

Do not continue tuning these two SEQC2 datasets to eliminate individual errors.
Freeze the candidate and validate against an independent completed dataset,
with matching truth/reference and a declared target domain, without retuning.
Before treating it as a workflow default, also resolve unrestricted candidate
discovery, generic sample roles, the evidence-based indel rule, and annotated
FILTER/rationale parity. No mapping or variant calling is needed for these
VCF-level checks when the necessary completed caller outputs already exist.

This audit has not inspected or tuned HG008, changed defaults, or modified
original inputs, workflow publications, archives, or caches.
