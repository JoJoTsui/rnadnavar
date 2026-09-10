# SEQC2 consensus and rescue FP/FN investigation

Status: interview and evidence collection; no implementation approved by this document.

## Scope and preserved constraints

Inspect DNA consensus, initial RNA consensus, realigned RNA consensus, first rescue and second-round rescue from the latest WES_LL hybrid output. Compare against DNA DeepSomatic using the same truth, reference and high-confidence/target regions. Retain SNP, indel and records metrics separately.

Existing user constraints remain: biological FILTER=Somatic defines selection; obsolete RaVeX_FILTER is not evidence; original inputs and variant-calling cache remain unchanged; the original seq2neo FASTQ workflow must remain runnable. Population AF, caller evidence and sample-specific AD/DP/VAF require provenance and missingness checks before use.

## Verified starting evidence

Source: examples/seqc2/hybrid/comparison/realignment_rescue_clair_latest/WES_LL_T_1_vs_WES_LL_N_1/WES_LL_T_1_vs_WES_LL_N_1.benchmark_comparison.csv.

| Callset | TP records | FP records | FN records | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| DNA consensus | 973 | 97 | 1327 | 0.9093 | 0.4230 | 0.5774 |
| Final rescue | 1011 | 314 | 1289 | 0.7630 | 0.4396 | 0.5578 |
| DNA DeepSomatic | 1048 | 38 | 1252 | 0.9650 | 0.4557 | 0.6190 |

Final rescue has net +38 TP and +217 FP relative to DNA consensus. Net differences do not establish incremental rescue yield: paired variant transitions must identify retained, gained and lost TP/FP separately. The table does not include initial RNA consensus, realigned RNA consensus or first rescue; their provenance and comparable benchmark availability need verification.

## Investigation tree

1. Verify artifact provenance and matching semantics for all five outputs.
2. Build paired error transitions, including DeepSomatic TP lost by each output.
3. Attribute errors to caller/class agreement, tumor and normal evidence, population AF, RNA editing, alignment round, variant class and genomic context.
4. Distinguish absent caller candidates from candidates rejected by downstream rules. RNA-wide FN includes unexpressed and unselected loci; retain the common-domain score and separately report ascertainment limitations.
5. Test proposed rules offline against frozen caller evidence, with independent evaluation partitions before judging improvements.

## Open decisions — round 1

- Optimization acceptance: recommend higher TP with no increase in FP relative to DNA DeepSomatic on held-out evaluation; report SNP and indel separately, and do not claim success from F1 alone.
- Validation: recommend fixed chromosome/block development and held-out evaluation partitions, documenting that a split within one cell line is not independent biological validation.
- Evaluation policy for RNA: recommend keep the same HC/target primary domain for comparability and add secondary expression/coverage/candidate-domain diagnostics without substituting an easier primary denominator.

No thresholds or causal conclusions have been established by these aggregate counts.
