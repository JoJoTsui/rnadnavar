# Spec: Consensus/Rescue Label Quality — Fixes, QC Gate, and Re-Consensus Rerun

**Status:** ready-for-agent
**Date:** 2026-08-28
**Tracker:** local files (GitHub Issues disabled on the fork)

## Problem Statement

EnsembleVar (this repo, forked from nf-core/rnadnavar) produces the consensus/rescue VCFs that serve as **truth labels** for training a downstream DN+DT+RT deep-learning somatic caller (`neo_var`, `set_somatic` tryouts). Final evaluation during model training surfaced three **abnormal samples** — `4081_rnadnavar`, `4255_rnadnavar`, `PRJNA298330_4032` — and an external audit (TruthQC) showed label cleaning alone moves model F1 from 0.633 to 0.908.

Root causes are in the pipeline, confirmed by adversarial review with empirical reproduction:

- **Label self-contradiction** (4081/4255): records carry `FILTER=Somatic` while their own unified per-caller evidence says Reference/Germline (74–80% of records in those samples). Genotype/DP/VAF aggregates are extracted from the **normal** sample, not the tumor; caller-support counts include caller-rejected records; the RaVeX filtering flags fire on 100% of records because they read a FORMAT column the consensus VCF does not have; the gnomAD germline guard is silently inert due to a field-name case mismatch.
- **Normal contamination** (4081/4255): the paired normal carries alt support (alt-VAF ≥5%) at most "somatic" truth sites — undetectable from VCFs alone.
- **RNA-only germline leakage** (PRJNA298330_4032): records supported only by RNA callers at common population frequencies get labelled Somatic; the rescue stage can never actually rescue (cross-modality agreement returns NoConsensus), RNA can override a DNA Artifact label, and `RESCUED`/`PASSES_CONSENSUS_*` flags mean "present in the union file", not "passed".

Meanwhile both model tryouts consume the VCF **FILTER column directly as labels** and never propose their own candidate sites — so label corruption translates one-to-one into model failure (observed as Reference bleed / score inflation).

## Solution

Three coordinated workstreams, all forward-only (existing cohort outputs stay read-only):

1. **Fix the consensus/rescue label logic at the source**: tumor-sample-aware genotype extraction, truthful filter flags, working gnomAD guard, correct support counting, a rescue contract that can actually rescue (cross-modality promotion, DNA-Artifact veto, passed-not-present flags), FILTER vocabulary preserved (model contract) with enriched INFO.
2. **A standalone label-QC gate** (`label_qc.py`): a port-and-expansion of TruthQC's proven R1–R6 rules plus cohort-adaptive sample-level gates and variant-level checks (including `COSM_` INFO fields and normal-contamination verification via samtools), emitting a 4-part report contract and per-sample PASS/WARN/FAIL verdicts for training inclusion.
3. **A re-consensus rerun path**: regenerate consensus and rescue VCFs for the whole cohort from the existing, read-only per-caller VCFs using the fixed code — consensus and rescue steps only; FASTQ→BAM alignment and per-caller variant calling are never re-triggered — into a new output location, followed by the QC gate to select high-confidence samples/variants for training.

## User Stories

1. As a model developer, I want truth labels whose FILTER agrees with the record's own per-caller evidence, so that the Somatic class is not poisoned by self-contradictory labels.
2. As a model developer, I want VAF/DP/GT aggregates computed from the tumor sample, so that quantitative label features and tiering logic reflect tumor biology.
3. As a model developer, I want a per-sample PASS/WARN/FAIL verdict before training, so that abnormal samples (4081/4255/4032 class) never reach the training set.
4. As a model developer, I want per-variant confidence flags within retained samples, so that negative sampling (Reference/Germline classes) is as clean as the positives.
5. As a model developer, I want normal-contamination detection, so that samples whose paired normal carries tumor signal are excluded — the models literally read DN reads.
6. As a model developer, I want RNA-only Somatic calls at common population frequencies flagged, so that germline polymorphisms do not leak through the RNA rescue branch.
7. As a model developer, I want the QC verdicts machine-readable, so that training-split scripts can consume them without parsing prose.
8. As a model developer, I want the label FILTER vocabulary unchanged, so that both model tryouts keep working without a coordinated migration.
9. As a pipeline maintainer, I want caller-support counts to include only non-Artifact caller records, so that consensus thresholds mean what they say.
10. As a pipeline maintainer, I want the rescue stage to actually rescue variants when DNA and RNA callers agree, so that cross-modality support is real signal.
11. As a pipeline maintainer, I want DNA Artifact labels to veto RNA overrides, so that RNA-specific error modes cannot flip DNA-flagged artifacts to Somatic.
12. As a pipeline maintainer, I want `RESCUED`/`PASSES_CONSENSUS_*` flags to mean "passed", so that downstream tiering and stats are not inflated.
13. As a pipeline maintainer, I want the gnomAD germline guard to actually fire, so that common polymorphisms cannot survive rescue filtering.
14. As a pipeline maintainer, I want filtering-stage flags derived from real alt counts, so that a 100%-bogus flag column never misleads downstream consumers again.
15. As a pipeline maintainer, I want CLI consensus thresholds to reach the rescue classifier, so that configuration is not silently ignored.
16. As a pipeline maintainer, I want the consensus entry step to guarantee normalized input, so that indel representation mismatches cannot silently erase caller agreement.
17. As a pipeline maintainer, I want an audit report with evidence for every finding, so that future changes can be checked against documented failure modes.
18. As a pipeline maintainer, I want a rerun driver that regenerates only consensus+rescue outputs, so that the cohort gets clean labels without re-paying alignment and calling compute.
19. As a pipeline maintainer, I want the rerun to treat existing outputs as strictly read-only, so that the original cohort results remain reproducible and comparable.
20. As a pipeline maintainer, I want regression tests reproducing each fixed bug on synthetic fixtures, so that the bugs cannot silently return.
21. As a cohort operator, I want QC thresholds that adapt to cohort statistics with absolute biological floors, so that the gate works on future cohorts without retuning.
22. As a cohort operator, I want optional BAM-based verification that does not require pysam, so that deep checks stay fast at cohort scale.
23. As a cohort operator, I want the QC tool to check COSMIC (`COSM_`) annotations against population frequency, so that hotspot-coincidence germline leakage is caught.
24. As a cohort operator, I want cleaned VCFs produced only under an explicit apply flag, so that QC never mutates inputs by default.
25. As a repo reader, I want the README to document the seq2neo cohort and the actual workflow, so that the repo explains what it is really used for.
26. As a future contributor, I want an optimization guidelines document covering the deferred findings (germline-leak rule ordering, RNA-editing over-masking, multi-sample rescue contamination, FILTER→PASS+INFO migration), so that known issues are not rediscovered by accident.

## Implementation Decisions

- **Naming**: "EnsembleVar" is this pipeline/repo. The downstream model is a separate consumer; nothing in the model repos changes in this work.
- **FILTER vocabulary preserved**: FILTER keeps the biological classes `{Somatic, Germline, Reference, Artifact, NoConsensus, RNAedit}` because both model tryouts read FILTER as the label contract. The FILTER-vs-INFO contradiction is resolved by making FILTER derivable from the unified classification, not by migrating to PASS+INFO (that migration is documented as a future, coordinated three-repo change in the guidelines). INFO is enriched instead: per-caller tumor AD/DP/VAF, classification rationale.
- **Tumor-sample resolution**: genotype/AD/DP/VAF extraction resolves the tumor sample per caller convention (Strelka: sample named TUMOR; Mutect2: via the normal-sample metadata; DeepSomatic: first sample). Strelka's max-read-count row selection is replaced by explicit tumor resolution.
- **Filter-flag provenance**: consensus output carries per-caller tumor alt-count INFO so downstream filtering stages never need to read a FORMAT column that consensus VCFs do not have; biological-class FILTER values are exempt from "caller rejected it" checks.
- **gnomAD guard**: population-frequency field lookup is case-insensitive / standardized on one canonical field name, with a regression test on an annotated VCF.
- **Support counting**: caller support counts only records the caller itself did not reject (non-Artifact), matching the documented consensus rules.
- **Rescue contract**: cross-modality promotion is enabled — when at least one DNA caller and at least one RNA caller agree on Somatic (thresholds configurable), the variant is rescued as Somatic and tagged in INFO; a DNA Artifact label vetoes RNA-driven overrides; `RESCUED`/`PASSES_CONSENSUS_*` flags are computed from records that passed consensus as Somatic, not from presence in the union file; CLI `--snv_thr/--indel_thr` reach the rescue classifier.
- **Rerun scope**: the rerun consumes existing per-caller VCFs (DNA and RNA branches) and BAMs read-only and executes only consensus → rescue → rescue post-processing (with normalization guaranteed at the entry step). FASTQ→BAM alignment and per-caller variant calling are structurally unreachable from the rerun path. Outputs go to a new location.
- **QC tool**: standalone Python CLI in `bin/` named `label_qc.py`, following repo conventions. Two tiers: Tier A is fast and VCF-only (ported TruthQC R1–R6 plus expanded sample gates — somatic-count outlier, self-contradiction rate, RNA-only fraction, modality completeness, spectrum sanity, coverage floors — and variant-level checks including strand bias context, low-complexity context, clustered variants, RNA-editing overlap, and COSM_ field cross-checks; sample-level rules also get variant-level counterparts where meaningful). Tier B is optional BAM-based verification (normal contamination, strand bias) via samtools subprocess — no pysam; a Rust/PyO3 backend next to the existing stats_core precedent is the designated extension point if profiling demands it.
- **QC thresholds**: adaptive per-cohort statistics (median/MAD outlier gates) with absolute biological floors as backstops; all thresholds in a bundled, overridable config.
- **QC output contract** (four parts): human report; machine-readable summary + per-sample verdict table; per-site flagged detail; cleaned VCFs only under an explicit apply flag. Inputs are never modified in place.
- **Docs**: README gains "Cohort (seq2neo)" and "Workflow" sections and a corrected caller list (Mutect2, Strelka2, DeepSomatic); guidelines live in `docs/`; the audit report lives in `dev_docs/audit/`; the QC design doc lives in `dev_docs/implementation/`.
- **Acknowledged incompleteness**: the fixed set is known to be a subset of all issues; the QC gate is the standing safety net, and the guidelines carry a known-unknowns section.

## Testing Decisions

- **What makes a good test**: tests assert on external behavior of the invoked artifact (output VCF records, INFO/FILTER content, verdict tables), never on internal call graphs. Each fixed bug gets a regression test that reproduces the original failure on a minimal synthetic fixture — the review already demonstrated this technique catches C1/C2-style bugs.
- **Primary seam — Python CLI**: the bin scripts (`run_consensus_vcf.py`, `run_rescue_vcf.py`, `filter_vcf.py`, `filter_rescue_vcf.py`, `label_qc.py`) invoked as CLIs on small synthetic caller-VCF fixtures under pytest. Prior art: the existing `tests/vcf_utils/` and `tests/vcf_stats/` pytest suites. This single seam covers the tumor-extraction fix, filter-flag fixes, gnomAD guard, support counting, the rescue contract, and every QC rule.
- **Secondary seam — nf-test at module level**: only where Nextflow wiring changes (normalization gating at the consensus entry step; rescue workflow contract). Prior art: `modules/local/*/tests/`.
- **Tertiary — smoke verification**: a scripted (non-CI) check running the rerun driver on 1–2 real samples, asserting new consensus+rescue VCFs are produced and input checksums are unchanged.
- No new seams are introduced; the Python CLI seam is the highest seam that exercises the label logic.

## Out of Scope

- Modifying, moving, or deleting any existing cohort output (all 66 samples' results remain read-only).
- Re-running FASTQ→BAM alignment or per-caller variant calling for any sample.
- Changes to the `neo_var` / `set_somatic` model repos (their contracts are inputs to this spec, not targets).
- The FILTER→PASS+INFO standards migration (documented as future work; requires coordinated change across three repos).
- The deferred findings M6 (annotation-stage germline rule ordering), M7 (RNA-editing over-masking), M10 (multi-sample cross-product rescue) — spec'd in the guidelines, not fixed in this round.
- SAGE/Manta module work (unused by the active workflow).
- Porting TruthQC itself into the repo (it remains the external reference implementation; `label_qc.py` is a fresh, expanded implementation).

## Further Notes

- Evidence for every finding, with reproduction transcripts, is deliverable of the audit-report ticket and is the authoritative reference for the fix tickets.
- The external TruthQC run (`truth_qc_out`, run id 20260826_204134_29d8e75350) is the baseline the QC gate must reproduce and exceed: 66 samples, rules R1–R6, F1 0.633→0.908 after cleaning; the three abnormal samples' signatures (self-contradiction rate, RNA-only fraction, normal alt-VAF) are the calibration set for the new gates.
- The two model tryouts differ in class-index order (neo_var: Germline/Somatic/Reference = 0/1/2; set_somatic: Reference/Germline/Somatic = 0/1/2) — any shared tooling must read class order from config.
- The 4 disease-exclusive folds (`set_number`) exist in the sample manifest but are unused by both model tryouts (chromosome splits share patient genomes across train/test); the guidelines recommend wiring patient-holdout CV.
