# 09: README — cohort + workflow sections

**What to build:** The root README gains two sections, in English, surgically inserted: (a) **Cohort (seq2neo)** — 66 eligible patients from SRA projects PRJNA298376/298330/298310, three modalities per patient (DN/DT/RT), WES + RNA-seq on GRCh38, partitioned into 4 disease-exclusive cross-validation folds, with the wrapper-script entry points; (b) **Workflow** — the actual run flow (manifest parsing → per-patient 3-modality pipeline runs → statistics), the utilized caller set (Mutect2, Strelka2, DeepSomatic), WES mode, and the consensus/rescue stages. The stale caller list (SAGE) is corrected, and the EnsembleVar naming is introduced with attribution to the nf-core/rnadnavar origin.

**Blocked by:** None (can start immediately).

**Status:** ready-for-agent

- [ ] Cohort section matches the facts in the seq2neo examples (patient counts, projects, modalities, folds)
- [ ] Workflow section reflects the real invocation path, not the stock nf-core template
- [ ] Caller list reads Mutect2, Strelka2, DeepSomatic; SAGE/Manta noted as unused
- [ ] EnsembleVar naming introduced with fork attribution
- [ ] No unrelated README churn
