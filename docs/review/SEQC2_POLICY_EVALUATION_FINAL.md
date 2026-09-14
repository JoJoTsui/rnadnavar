# SEQC2 policy evaluation result

The four-cell VCF-only experiment completed for ordinary 2-of-3 consensus, native-SNV consensus and strict 3-of-3 consensus. All 36 policy/class rows were evaluated against Mutect2, Strelka2 and DeepSomatic with precision at least DeepSomatic and F1 strictly above every DNA comparator required in every cell.

The result is **no qualifying shared policy**. The native policy is the strongest SNP/record candidate in most cells, while ordinary consensus is strongest for indels, but neither satisfies every cell. The WES-LL cells remain the limiting cases; the WGS-IL cells are closer but still fail the strict DeepSomatic F1 gate.

The combined policy therefore remains an evaluation result, not a default change:

```text
SNV: native evidence candidate
indel: ordinary caller threshold candidate
rescue: retained as a separate stage; no broad promotion adopted
release: rejected/inconclusive until a new evidence rule passes all cells
```

The cached rescue transitions explain why rescue is not promoted automatically. First and realignment rescue add more false than true records in both datasets and can remove true baseline records. Filtering and final annotation must remain separate from voting attribution.

Machine-readable outputs:

- `examples/seqc2/comparison/common_policy_20260914/policy_experiments/evaluation.json`
- `examples/seqc2/comparison/common_policy_20260914/*/*/*.transitions.json`

No mapping, variant calling, source VCF or cache was changed. HG008 is not used for tuning; it remains the next independent validation stage after a qualifying policy exists.
