# SEQC2 cross-dataset policy runbook

This runbook describes the reproducible benchmark and policy-freeze boundary. It does not certify a policy before the four-cell matrix has been evaluated.

## Decision flow

```text
caller VCFs
  ├─ DNA callers ──> DNA consensus ──> raw/filtered audit
  ├─ RNA callers ──> RNA consensus
  ├─ realigned RNA callers ──> RNA-realign consensus
  └─ DNA + RNA evidence ──> first rescue ──> realignment rescue ──> final annotation

Each arrow preserves caller/sample/round provenance.
Missing RNA evidence is unavailable evidence; it is not a negative vote.
Indels are evaluated under their own evidence rule; the SNP native override does not admit them.
```

## Matrix execution

Use existing caller-ready output directories. The command performs VCF preparation and `som.py` benchmarking only:

```bash
bash examples/seqc2/scripts/run_seqc2_policy_matrix.sh \
  <comparison-root> \
  <wes-ll-output> WES_LL_T_1_vs_WES_LL_N_1 \
  <wgs-il-output> WGS_IL_T_1_vs_WGS_IL_N_1
```

The output preserves `wes_ll/ukb`, `wes_ll/medexome`, `wgs_il/ukb` and `wgs_il/medexome`. The benchmark normalizes truth and query by default; set `NORMALIZE_ALL=0` only to reproduce a historical command. Optional RNA and rescue controls are passed through `EXTRA_QUERIES=name=vcf[,name=vcf]` and rescue through `RESCUE_VCF`.

## Acceptance

Run the matrix validator after all four cells exist:

```bash
python3 examples/seqc2/scripts/validate_policy_matrix.py \
  --matrix <comparison-root> --candidate consensus \
  --out <comparison-root>/policy_validation.json
```

The candidate must have higher F1 than each declared DNA comparator and precision at least matching DeepSomatic for SNP, indel and records in every cell. A failed gate is a valid investigation outcome; it does not authorize threshold changes.

## HG008 handoff

Freeze the exact policy, thresholds, command provenance and matrix result before evaluating HG008. Use the existing shared wrapper after the HG008 run is complete:

```bash
bash examples/seqc2/hybrid/run_hg008_wgs_hybrid_benchmark_shared.sh
```

This handoff must use the completed final realignment-rescue VEP artifact and its declared truth/reference contract. Do not modify the active HG008 workflow or use its truth-derived errors for policy tuning.
