# Combined native-consensus + gated-rescue policy (available SEQC2 evidence)

This is the frozen non-HG008 combined-policy record. It combines native SNV consensus, the DNA-vetoed gated rescue (DNA nomination plus RNA corroboration), and the retained threshold-consensus indel policy. It does not change production defaults or launch a workflow.

| Cell | TP | FP | FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| wes_ll/ukb | 1062 | 38 | 1238 | 0.9655 | 0.4617 | 0.6247 |
| wes_ll/medexome | 570 | 19 | 259 | 0.9677 | 0.6876 | 0.8039 |
| wgs_il/ukb | 2169 | 19 | 131 | 0.9913 | 0.9430 | 0.9666 |
| wgs_il/medexome | 702 | 6 | 127 | 0.9915 | 0.8468 | 0.9135 |

Source: `examples/seqc2/comparison/rescue_fp_investigation_20260914/gate_tests/evaluation.json` (SHA-256 `eee0c7f0478f4dc26bdcd359fee6083a305b13af49d584dc3a0435da8fbb770f`).

The four available cells complete the bounded combined-policy evaluation. Indels remain threshold consensus because the shared indel comparator found no candidate satisfying TP/FP non-regression in all four cells. HG008 independent validation, production adoption, and any claims beyond these SEQC2 domains remain deferred.
