# Shared SEQC2 consensus and rescue policy investigation

Status: interview and evidence audit, 2026-09-14. No production changes or workflow runs performed.

## Accepted decisions

- Use one evidence-based policy across WES-LL and WGS-IL, without dataset-specific thresholds or truth-site exceptions.
- Require higher F1 than standalone DNA callers and precision at least matching DNA DeepSomatic on each dataset. Report SNP, indel and combined-record TP/FP/FN, precision, recall and F1 separately for consensus, first rescue and realignment rescue.
- Genomic holdouts are exploratory: SEQC2 assays share a biological sample and WES was previously tuned. Freeze the policy before independent HG008 evaluation.
- Reuse caller and annotation artifacts for isolated experiments after the interview; preserve original inputs, outputs and caches. Leave running HG008 unchanged.
- Missing RNA coverage is unavailable evidence. Investigate SNP and indel rules separately.

## Confirmed comparison gaps

1. The Indel candidate audit in `SEQC2_OPTIMIZATION_FOLLOWUP.md` records experimental indels anchored to the strongest baseline caller. Its implementation section leaves indels threshold-based. The final WES report mixes experimental metrics with the deployed description.
2. WGS metrics record `som.py -N` (normalize both inputs); authoritative WES native/gated commands omit it. Numerical effects remain unmeasured.
3. WGS uses HC plus UKB targets; authoritative WES uses HC plus MedExome. These are distinct benchmark domains.
4. Historical gated rescue retained every native baseline record. WGS compares raw DNA consensus with final filtered/annotated rescue whose DNA input was already filtered. Stage attrition remains unresolved.
5. The WGS table omits first rescue and RNA branch controls.

Sources: `docs/SEQC2_CONSENSUS_RESCUE_FINAL_REPORT.md`; `docs/review/SEQC2_OPTIMIZATION_FOLLOWUP.md`; `examples/seqc2/hybrid/comparison/authoritative_medexome_bed_20260910/`; `examples/seqc2/hybrid/comparison/WGS_IL_T_1_vs_WGS_IL_N_1/`.

WGS records TP/FP/FN: DeepSomatic 2168/19/132; consensus 2165/22/135; final rescue 2144/25/156. SNP consensus gains 9 TP and adds 10 FP versus DeepSomatic; indels lose 12 TP and remove 7 FP. These are net counts, not paired allele transitions.

## Accepted comparison and indel scope

Q4 and Q5 accepted: benchmark WES-LL and WGS-IL in both UKB and MedExome regions and preserve every result. Use identical truth/reference, normalization and selection semantics within each domain. Reproduce deployed policies separately from historical experiments and trace raw/filtered consensus, first rescue, realignment rescue and final annotation before tuning.

Examine indel rules across datasets: allele representation, insertion/deletion length, repeats, caller rejection, tumor/normal evidence and missingness. Select one shared evidence-based rule against the agreed precision and F1 criteria, retaining all four dataset/domain results. Strongest-caller indel copying is not an accepted replacement.

Preserve attempted rules, commands, input identities, thresholds, paired TP/FP/FN transitions and negative results in separate experiment outputs. Do not use pooled scores to hide a dataset or region regression. If no candidate meets the criteria, report that outcome. Freeze the policy before HG008 evaluation.

Interview scope is settled. Next: common-contract baseline reproduction and paired error attribution using existing artifacts only.
