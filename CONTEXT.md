# Context / Glossary

Domain terms for this repository. Glossary only — no implementation details or decisions.

## Tools and projects

- **EnsembleVar** — this repository: the ensemble somatic variant-calling pipeline forked from nf-core/rnadnavar. Combines multiple variant callers (Mutect2, Strelka2, DeepSomatic) over matched DNA/RNA data. (Resolved: "EnsembleVar" names the pipeline/repo, not the downstream model.)
- **TruthQC** — external standalone tool (outside this repo) that audits truth-label VCFs for mislabelled germline-like sites, using rules R1–R6 over VCF INFO/FILTER fields. Reference implementation for label QC.

## Samples and data

- **DN / DT / RT** — the three modalities of one patient: DNA normal (sample status 0), DNA tumor (status 1), RNA tumor (status 2).
- **DN/DT/RT FASTQ triplet** — one paired-end raw-read input for each of DN, DT, and RT: three biological input rows and six FASTQ files. _Avoid_: “three FASTQs,” which confuses sample roles with file count.
- **Pooled RT sample** — one logical RNA-tumor sample formed from explicitly identified RNA library repeats while retaining each repeat's provenance. It is not evidence that the libraries came from the same extraction or aliquot as the DNA samples.
- **RNA library repeat** — an independently prepared RNA-seq library for the same biological sample. It remains distinct from a sequencing lane so duplicate handling and provenance preserve the library boundary.
- **Cell-line-matched hybrid benchmark** — an integrated DNA/RNA analysis whose inputs represent the same cell-line identities but are not proven to share a specimen, extraction, or aliquot. Its rescue output is benchmark evidence, not specimen-matched training truth.
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
- **Hybrid-ingress run** — one patient analysis whose modalities enter from different preparation stages, such as caller-ready DN/DT alignments together with raw RT reads, and converge before variant calling. _Avoid_: mixed-input run, which does not say that the inputs begin at different stages.
- **Caller-ready alignment** — an indexed, coordinate-sorted alignment whose sample identity, preparation state, and reference compatibility are sufficient for direct variant calling without implicit remapping or duplicate marking.
- **Input stage** — the declared preparation state at which an input enters the pipeline: raw reads, an alignment intended for remapping, or a caller-ready alignment. It is independent of whether the containing file is FASTQ, BAM, or CRAM.
- **Hybrid input manifest** — a samplesheet in which every row explicitly declares its input stage so multiple preparation stages can coexist without inference. _Avoid_: partially staged manifest.
- **Ingress provenance** — the auditable record of how each source input reached the caller-ready boundary, including its identity, preparation stage, reference decision, and library relationship.
- **Shared callable region** — the genomic region in which DNA and RNA evidence are intentionally compared for integrated consensus and rescue labels.

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
