# Context / Glossary

Domain terms for this repository. Glossary only — no implementation details or decisions.

## Tools and projects

- **EnsembleVar** — this repository: the ensemble somatic variant-calling pipeline forked from nf-core/rnadnavar. Combines multiple variant callers (Mutect2, Strelka2, DeepSomatic) over matched DNA/RNA data. (Resolved: "EnsembleVar" names the pipeline/repo, not the downstream model.)
- **TruthQC** — external standalone tool (outside this repo) that audits truth-label VCFs for mislabelled germline-like sites, using rules R1–R6 over VCF INFO/FILTER fields. Reference implementation for label QC.

## Samples and data

- **DN / DT / RT** — the three modalities of one patient: DNA normal (sample status 0), DNA tumor (status 1), RNA tumor (status 2).
- **seq2neo cohort** — the 66-patient WES + RNA-seq cohort assembled from SRA projects PRJNA298376, PRJNA298330, PRJNA298310, partitioned into 4 disease-exclusive cross-validation folds. Its pipeline outputs are the training labels for the downstream model.

## Pipeline concepts

- **Consensus** — within-modality merging of per-caller VCFs using caller-support thresholds.
- **Rescue** — cross-modality (DNA ↔ RNA) recovery of variants that failed or were missed in one modality.
- **Truth label** — a final Somatic-filtered consensus/rescue VCF record used as a supervised training label for the downstream model. The model is a label consumer: label precision matters more than recall.

- **Re-consensus rerun** — forward-only regeneration of consensus and rescue VCFs for the whole cohort from the existing, read-only per-caller VCF outputs, using fixed pipeline code, written to a new output location. Original outputs are never modified.

## Label-quality failure classes (observed)

An **abnormal sample** is one whose truth labels fail QC at scale. Three distinct classes observed so far:

- **Label self-contradiction** — VCF `FILTER=Somatic` while the record's own unified per-caller evidence (e.g. `UNIFIED_FILTER`) says Reference/Germline. Seen massively in samples 4081 and 4255.
- **Normal contamination** — the paired normal carries alt support at "somatic" truth sites (alt-VAF ≥ 5% in normal at most FN sites). Also seen in 4081/4255.
- **RNA-only germline leakage** — records supported only by RNA callers (no DNA-caller entry) at common population frequencies (high gnomAD AF, heterozygous-like VAF) labelled Somatic. Seen in PRJNA298330_4032.
