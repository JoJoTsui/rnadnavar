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
- **Reference-compatible BAM** — an external BAM whose sequence dictionary has exactly the selected reference's contig names, lengths, and order. Off-reference contigs are outside the caller-input domain.
- **Reference normalization** — the audited, provenance-preserving transformation of a safely compatible external BAM into a reference-compatible BAM before any caller.
- **Caller input pair** — the single reference-compatible tumor/normal alignment pair shared by Mutect2, Strelka2, and DeepSomatic.
- **Expected caller panel** — the complete caller set required for one sample/modality consensus. Missing, duplicate, or unexpected inputs make the panel incomplete rather than changing its size.
- **Caller-detection support** — the `ENS_SUPPORT=k/n` count of expected callers whose records support the variant's existence. It is distinct from biological-class agreement and calibrated label confidence.
- **DNA-only run** — a matched DNA tumor/normal run through DeepSomatic, Mutect2, Strelka2, normalization, and within-DNA consensus, with no RNA rescue branch.

- **Re-consensus rerun** — forward-only regeneration of consensus and rescue VCFs for the whole cohort from the existing, read-only per-caller VCF outputs, using fixed pipeline code, written to a new output location. Original outputs are never modified.

## Split and downstream-task concepts

- **Working cohort** — the 63 samples eligible for dataset splitting: rows of the rerun sample manifest with a non-empty training label VCF (56 PASS + 7 WARN). The 3 `useless` samples are never split.
- **Selected variant** — a truth-label VCF record whose FILTER is in the whitelist {Somatic, Germline, Reference}; only selected variants enter dataset manifests. (NoConsensus, Artifact, RNAedit are excluded.)
- **Sample pool** — the sample-level assignment: *reserved* (held out wholly for downstream evaluation) or *train pool*. For train-pool samples the real split is per-variant, by chromosome.
- **Reserved pool** — 5 hardcoded non-digestive PASS samples spanning the disease folds, used for downstream evaluation (zero-shot and tag-defined sub-pools). Never used in training.
- **Chromosome split (deepsomatic)** — per-variant assignment for PASS and WARN train-pool samples: chr1 → test, chr21–22 → val, chr2–20 → train, any other chromosome → train.
- **Verdict routing** — PASS and WARN samples use the same chromosome split for non-reserved variants: chr1 → test, chr21–22 → val, and chr2–20/other chromosomes → train. PASS is the primary evaluation stratum; WARN is retained as a separate sensitivity stratum.
- **Sub-pool tags** — per-variant booleans defining downstream evaluation subsets: is_zero_shot, is_low_vaf_a/b, is_low_dp, is_rescued, is_non_rescued, is_indel. Tags are derived from the truth-label VCF's own per-modality depth/VAF/rescue INFO fields.
- **Split manifest** — the pair of artifacts recording pools, per-variant splits, tags, and label verdicts; the contract for downstream dataset validation and path tracking.

## Label-quality failure classes (observed)

An **abnormal sample** is one whose truth labels fail QC at scale. Three distinct classes observed so far:

- **Label self-contradiction** — VCF `FILTER=Somatic` while the record's own unified per-caller evidence (e.g. `UNIFIED_FILTER`) says Reference/Germline. Seen massively in samples 4081 and 4255.
- **Normal contamination** — the paired normal carries alt support at "somatic" truth sites (alt-VAF ≥ 5% in normal at most FN sites). Also seen in 4081/4255.
- **RNA-only germline leakage** — records supported only by RNA callers (no DNA-caller entry) at common population frequencies (high gnomAD AF, heterozygous-like VAF) labelled Somatic. Seen in PRJNA298330_4032.
