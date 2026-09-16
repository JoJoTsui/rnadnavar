# Optimized consensus and rescue comparison

Date: 2026-09-16  
Branch: `seqc2-consolidated`

This record supersedes comparisons that used the unrestricted historical
native-plus-gated replay. The optimized rescue policy is the frozen
DNA-nomination plus biological-veto gate from
`SEQC2_RESCUE_FP_INVESTIGATION.md`:

```text
retain every optimized native-consensus record
add only SNVs with RNA support >= 2 and a positive DNA nomination
reject common gnomAD alleles and the DNA-unsupported canonical editing case
keep indels on the validated threshold-consensus path
```

All values below are `som.py -N` record metrics using the same truth, HC BED,
target BED, and GRCh38 reference within each dataset. TP/FP/FN are records;
SNP and indel metrics remain in the source JSON files.

| Dataset / method | TP | FP | FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| SEQC2 WES-LL optimized consensus + gated rescue | 1062 | 38 | 1238 | 0.9655 | 0.4617 | **0.6247** |
| SEQC2 WES-LL optimized consensus | 1051 | 36 | 1249 | 0.9669 | 0.4570 | 0.6206 |
| SEQC2 WES-LL DNA DeepSomatic | 1048 | 38 | 1252 | 0.9650 | 0.4557 | 0.6190 |
| SEQC2 WES-LL realignment rescue | 1011 | 314 | 1289 | 0.7630 | 0.4396 | 0.5578 |
| SEQC2 WES-LL first-round rescue | 1011 | 515 | 1289 | 0.6625 | 0.4396 | 0.5285 |
| SEQC2 WGS-IL optimized consensus + gated rescue | 2169 | 19 | 131 | 0.9913 | 0.9430 | **0.9666** |
| SEQC2 WGS-IL optimized consensus | 2169 | 19 | 131 | 0.9913 | 0.9430 | 0.9666 |
| SEQC2 WGS-IL DNA DeepSomatic | 2168 | 19 | 132 | 0.9913 | 0.9426 | 0.9663 |
| SEQC2 WGS-IL realignment rescue | 2144 | 25 | 156 | 0.9885 | 0.9322 | 0.9595 |
| SEQC2 WGS-IL first-round rescue | 2141 | 25 | 159 | 0.9885 | 0.9309 | 0.9588 |
| HG008 WGS optimized consensus | 499 | 20 | 170 | 0.9615 | 0.7459 | **0.8401** |
| HG008 WGS DNA DeepSomatic | 472 | 24 | 197 | 0.9516 | 0.7055 | 0.8103 |
| HG008 WGS first-round rescue | 490 | 150 | 179 | 0.7656 | 0.7324 | 0.7487 |

HG008 realignment rescue is not included: the second-round artifact is not
available. Consequently, HG008 has not yet validated optimized rescue, only
optimized consensus and the unoptimized first-round rescue.

## Ranking

1. WES-LL optimized consensus + gated rescue (F1 0.6247).
2. WGS-IL optimized consensus and gated rescue (tie, F1 0.9666).
3. HG008 optimized consensus (F1 0.8401).
4. DNA DeepSomatic is below the optimized method in all three completed
   comparisons.

Neither rescue round currently beats the optimized DNA baseline when used
   without the optimized gate. No standalone Mutect2 or Strelka2 result
   exceeds the primary methods under this contract.

## Provenance

- SEQC2 optimized gate: `examples/seqc2/comparison/rescue_fp_investigation_20260914/gate_tests/{wes_ll,wgs_il}/nomination_biological/`.
- SEQC2 optimized consensus and DeepSomatic controls:
  `examples/seqc2/comparison/historical_manual_replay_20260914/`.
- WGS first/realignment rescue benchmark-only metrics:
  `/tmp/wgs-optimized-rounds-pass-20260916/`.
- HG008 metrics: `/tmp/hg008-rescue-benchmark-20260916/bench/`.

The benchmark copies are derived PASS queries. Original VCFs, workflow
outputs, mapping/calling caches, and input files were not modified.

## Next implementation boundary

Apply the frozen optimized gate to HG008's available first-round caller
evidence in a separate output namespace. Do not tune thresholds on HG008 and
do not treat realignment as an independent vote. Indel rescue remains disabled
until one shared indel rule improves all completed datasets without an FP
regression.
