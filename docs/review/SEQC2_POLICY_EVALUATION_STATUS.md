# SEQC2 policy evaluation status

Ticket 02 stage attribution is complete for cached WES-LL and WGS-IL rescue artifacts. Ticket 03–05 experiments use `run_consensus_policy_experiments.sh` and existing caller VCFs only. Candidate outputs are isolated by policy name (`ordinary`, `native`, `strict`) and include source provenance.

Candidate interpretation:

```text
ordinary: 2-of-3 caller threshold for SNVs and indels
native:   native-evidence SNV policy; indels remain 2-of-3
strict:   3-of-3 caller threshold for SNVs and indels
```

Each candidate is evaluated independently per UKB and MedExome target. The candidate gate requires higher F1 than every declared DNA comparator and precision at least matching DeepSomatic for SNP, indel and records in every cell. No candidate is adopted when any cell fails.

The current deployed native consensus already fails this gate in all four cells; this baseline decision is preserved. Rescue remains evaluated separately because stage transitions show false additions exceeding true additions in the cached WES and WGS outputs.
