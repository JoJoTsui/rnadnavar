# EnsembleVar Optimization Guidelines

**Date:** 2026-08-29
**Status:** living document
**Scope:** future optimization of the EnsembleVar pipeline (consensus / rescue / filtering / classification) and its handoff to the downstream DN+DT+RT model. Companion documents: the adversarial audit (`dev_docs/audit/2026-08-28_adversarial_review.md`), the label-QC design (`dev_docs/implementation/label_qc_design.md`), the rerun runbook (`docs/RECONSENSUS_RERUN.md`).

This document carries (a) fix specs for findings deliberately deferred from the 2026-08 fix round, (b) the pipeline↔model contract rules learned from the two model tryouts, and (c) the known-unknowns posture.

---

## 1. Deferred findings — fix specs

### 1.1 Annotation-stage germline leakage (M6) — highest priority deferred item

**Problem.** In the annotation-stage classifier, Rule 1 reclassifies a variant to Somatic when ≥2 DNA callers carry the label "Somatic" (>50% of DNA callers) with **no gnomAD check**, taking priority over the germline rule; and the germline-demotion rule additionally requires a Germline-labeled caller in *both* DNA and RNA, which is nearly unreachable because RNA Mutect2 is sparse. A caller's "Somatic" label is just its own `FILTER == PASS`; if Mutect2 runs without `--germline-resource`/panel-of-normals, common SNPs PASS and get labeled Somatic.

**Fix spec.**
1. First confirm configuration reality: check whether the cohort's Mutect2 invocations used `--germline-resource` (af-only-gnomad) and `--panel-of-normals` (1000g PON is present in the seq2neo shared config — verify it reaches the process). If they did, severity drops; if not, this is the structural leak behind the 4032-class failures on the DNA side.
2. Reorder the rules: evaluate the population-frequency (gnomAD) demotion **before or jointly with** Rule 1, so ≥2-caller agreement cannot override common-AF evidence.
3. Drop the RNA-germline requirement for demotion: DNA-only germline evidence at common AF must suffice.
4. Regression-test with a fixture: two DNA callers PASS a variant at gnomAD AF 0.3 → must demote to Germline.

### 1.2 RNA-editing over-masking (M7)

**Problem.** Any variant with an exact REDIportal match and ≥2 RNA callers is relabeled `RNAedit` — including MEDIUM tier (has DNA presence) and LOW tier (non-canonical transition at a known site). This removes genuine high-VAF DNA-supported A>G somatic mutations at editing sites, and the step runs *after* COSMIC/gnomAD reclassification, overwriting it.

**Fix spec.** Restrict FILTER changes to the VERY_HIGH/HIGH tiers (no DNA support); for MEDIUM/LOW tiers emit INFO annotation only (e.g. an `RNA_EDIT_TIER` field) and leave FILTER alone. Depends on the ticket-02 tumor-VAF fix being in place (it is), because the tier logic keys off DNA VAF. Regression-test: DNA-supported A>G at a REDIportal site with tumor DNA VAF 0.3 stays Somatic with an INFO annotation.

### 1.3 Multi-sample cross-product rescue (M10)

**Problem.** Rescue crosses DNA×RNA consensus by *patient* while individual caller VCFs are grouped by patient — for patients with >1 tumor or >1 RNA aliquot, each pair's rescue receives other samples' calls.

**Fix spec.** Group caller VCFs by pair id (`meta.id`), not patient; add an nf-test or fixture-level check with a two-tumor patient proving no cross-sample leakage. Low urgency for the current cohort (single tumor per patient) but mandatory before multi-sample patients are processed.

### 1.4 Minor findings (from the audit — address opportunistically)

- **QUAL averaging** across incomparable caller scales (Strelka EVS vs Mutect2 QUAL vs DeepSomatic) — drop or replace with a defined composite.
- **Homopolymer filter never runs**: the filtering config never passes `--ref`; also its deletion-context handling treats the anchor base as SNV context, and `alt == "-"` is dead MAF notation. Either wire `--ref` and fix the context logic, or remove the dead code.
- **Threshold drift**: gnomAD 0.0001 (rescue filter path, MAF path) vs 0.001 (annotation/filtering path) — pick one germline-frequency threshold family and document the rationale.
- **MAF-path indel rescue** matches any positional overlap regardless of allele (legacy `run_consensus.R` path) — require allele equality.
- **Filename-derived caller identity** is fragile (dots in sample names break parsing; duplicate caller names silently overwrite) — pass caller identity through metadata instead.
- **Header-copy logic** in the union writer strips other callers' INFO/FORMAT definitions — audit which fields are silently dropped for downstream consumers.
- **`--normal-sample` propagation**: tumor resolution currently infers from VCF sample names/ordering (ticket 02); the robust follow-up is a CLI flag carrying the normal sample name into the consensus/rescue scripts.
- **Model-side dead channel precedent**: `neo_var`'s `supports_alt` channel was dead (constant) due to a gating bug — when regenerating tensors after the rerun, verify channel liveness before training.

---

## 2. Pipeline → model contract rules

Learned from the two model tryouts (`neo_var`, `set_somatic`). Treat these as hard constraints on pipeline output changes:

1. **FILTER is the label.** Both tryouts read the VCF FILTER column directly, whitelisting `{Somatic, Germline, Reference}` and silently dropping everything else (NoConsensus/Artifact/RNAedit become *silent sample loss*, never errors). Any FILTER-vocabulary change is a three-repo coordinated migration (see §3).
2. **Class-index order differs between consumers** (neo_var: Germline/Somatic/Reference = 0/1/2; set_somatic: Reference/Germline/Somatic = 0/1/2). Shared tooling must read class order from config, never assume.
3. **The pipeline VCF is the candidate universe.** Neither model proposes candidate sites; precision/recall numbers are conditioned on the rescue VCF's prevalence. Keep the rescue VCF deterministic (dedup by chrom/pos/ref/alt, decompose multi-allelics) and treat any change to candidate generation as a label-contract change.
4. **Reference negatives come from the same VCF's FILTER=Reference rows** — Reference-label purity matters as much as Somatic precision (both tryouts' dominant failure was Reference bleed / score inflation). The label-QC gate must therefore check negative classes too, not just Somatic.
5. **AD/VAF features are circular with consensus labels.** The pipeline should keep emitting per-caller tumor AD/DP/VAF reliably in INFO (done — tickets 02/03/07), but models must not train on them (set_somatic excludes them deliberately; neo_var's DeepSomatic-style candidate filter uses them only for filtering).
6. **Patient-holdout CV is available but unused.** The manifest's `set_number` (4 disease-exclusive folds) is parsed but both tryouts split by chromosome *within* all patients — train/test share patient genomes (optimistic generalization). Wire `set_number` into the split stage for the final model.
7. **Depth gates exist downstream**: set_somatic drops sites with <2 reads per required modality (~4.2% of sites); neo_var tags low-DP pools at BAM_DT_DP < 20. Stable depth reporting in pipeline outputs avoids split-count drift.

---

## 3. Future standards migration: FILTER → PASS + INFO/VC

SOTA VCF hygiene would put PASS in FILTER and the biological class in INFO (the `VC` header field is already declared but never written). This is **blocked on the model contract** (§2.1): both tryouts read FILTER. Migration sequence when the time comes: (1) both model repos gain a config switch to read the class from INFO; (2) pipeline emits dual (PASS + INFO class, keeping FILTER classes one release behind a flag); (3) cut over after a full cohort rerun + label-QC pass. Until then, FILTER keeps the biological classes and `CLASSIFICATION_RATIONALE` (ticket 07) carries the audit trail.

---

## 4. Label-quality operations (standing procedure)

1. **Never train on ungated labels.** Run `bin/label_qc.py` (Tier A always; Tier B `--verify-bam` for any sample entering a training fold) over every new/re-run cohort. Training inclusion = PASS; WARN requires manual review; FAIL is excluded or remediated.
2. **Calibration anchors**: the three 2026-08 abnormal samples are the standing regression set for the QC gate — 4081/4255 (self-contradiction + normal contamination) and PRJNA298330_4032 (RNA-only common-AF leakage) must FAIL; 4278/4252/PRJNA298330_4096 must PASS. If a pipeline change alters these verdicts unexpectedly, stop and investigate.
3. **Re-consensus rerun** (this round's deliverable): regenerate labels via the consensus+rescue-only rerun path (`docs/RECONSENSUS_RERUN.md`), never by re-aligning or re-calling; originals stay read-only; checksum guards verify that.
4. **After any classification-logic change**, diff the label VCFs on a fixed sample subset (one clean + the three anchor abnormals) before any cohort-wide regeneration.

---

## 5. Known unknowns

- The 2026-08 fix set (C1–C3, M1–M5, M8a, M9 + plumbing) is **a subset of all issues** — it is what one adversarial review plus one external label audit surfaced. There is no proof the label logic is now correct; there is only a stronger net. The QC gate (`label_qc.py`) is the standing mitigation: keep expanding its rule set whenever a new failure class is diagnosed in model evaluation.
- **Normal contamination is a cohort provenance problem, not a pipeline bug.** Tier B detects it, but the source (sample swap? low tumor purity? shared variants?) needs wet-lab/data-provenance follow-up for 4081/4255.
- **Low-coverage FN (~28% of evaluation FNs, tumor DP < 50)** is a genuine model/data limitation, not label noise — address via candidate generation or model design, not label QC.
- **Mutect2 germline-resource/PoN usage in the cohort configs is unverified** (§1.1 step 1) — until confirmed, treat the germline-leak fix as speculative in impact.
- The **MAF-based legacy path** (`run_consensus.R`, `enable_maf_workflow`) was only lightly audited; if it is ever re-enabled, audit it first (indel positional-overlap rescue is known-bad).
- Both model tryouts diagnosed **score inflation vs score separation** as the dominant residual failure even on cleaned labels — expect the next failure class to come from label *ambiguity* (low-VAF subclonal vs artifact), not the gross corruption fixed here.
