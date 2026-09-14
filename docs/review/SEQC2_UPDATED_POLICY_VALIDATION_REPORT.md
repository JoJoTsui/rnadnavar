# SEQC2 updated-policy validation report

Date: 2026-09-14. Validation implementation commit: pending. This report is
generated from the frozen 20260914 evidence bundle and completed VCF-only
assays; it does not represent independent HG008 validation.

## Executive verdict

The frozen nomination-plus-biological rescue candidate reproduces its recorded
metrics and passes the available VCF-level contract checks. The current
workflow is **not equivalent** to that candidate: its enabled native/rescue
flags implement older semantics, including native majority fallback, different
indel handling, and no nomination/biological gate wiring. Therefore production
defaults remain unchanged and release approval is blocked pending parity work,
independent validation, and an explicit indel policy.

Machine-readable output:

`examples/seqc2/comparison/updated_policy_validation_20260914/validation_report.json`

The validator completed with status `pass_with_known_gates`; this means the
bounded checks passed while known release gates remain open. It does not mean
the workflow is approved for training-label generation.

## Checks completed

- Four candidate benchmark cells reproduced exactly: WES-LL/WGS-IL ×
  UKB/MedExome, with record counts 1062/38/1238, 570/19/259, 2169/19/131,
  and 702/6/127 respectively.
- Domain audit accounted for all 303 WES and 610 WGS rescue additions and
  found no truth-present exclusion within the evaluated HC domain.
- Root defaults were verified as native SNV on, rescue promotion on, and
  rescue floors DNA=1/RNA=2.
- The experimental nomination/biological gate remains isolated and marked
  exploratory; no production module contains its policy name.
- Re-consensus configs use `step: consensus`, consensus/rescue tools, no caller
  tools, and dedicated checksum namespaces.
- Seventeen focused tests passed, covering historical replay, rescue evidence,
  domain attribution, multiallelic/sample-role handling, and validator output.

## Material findings

1. **Policy mismatch (release blocker).** Production native consensus is
   additive and can fall back to ordinary majority; the frozen replay starts
   from the selected baseline and applies its exact veto ordering. Production
   indels remain threshold-consensus, while the historical benchmark retains a
   DeepSomatic-derived indel baseline.
2. **Rescue semantics mismatch (release blocker).** Production rescue can
   reclassify a DNA Somatic baseline during cross-modality disagreement. The
   frozen experiment retains the entire native baseline and filters only rescue
   additions.
3. **Vote-semantics defect risk.** The production rescue comment refers to DNA
   Somatic support, but the implementation counts eligible caller lists, which
   can include Reference observations.
4. **Annotation-order gap.** Population/editing annotations are produced in
   downstream rescue post-processing. The experimental biological veto uses
   final annotated rescue fields, so equivalent production placement is not
   demonstrated.
5. **Evidence limitations.** Outside-HC additions are unassessed; read-level
   evidence and independent HG008 validation are not yet available. RNA
   realignment is not independent biological confirmation.

## Input modes and cache safety

The static routing audit confirms that the same generic native/rescue flags are
forwarded for DNA-only BAM, hybrid caller-ready DNA plus RNA FASTQ, and DN/DT/RT
FASTQ-triplet paths. DNA-only has no RNA branch and therefore cannot produce a
cross-modality rescue. Consensus/rescue-only rerun configs structurally omit
caller names and use a consensus entry step; they write checksums and outputs
to dedicated namespaces. This is static reachability evidence, not live cache
reuse certification.

No mapping, variant calling, cohort rerun, HG008 execution, source modification,
or cache cleanup was performed by this validation.

## Required follow-on tickets

- Add a production seam that can represent the frozen policy without changing
  legacy behavior until explicitly selected.
- Align observed/eligible/Somatic support semantics and annotation timing, with
  regression fixtures for Reference, missing, duplicate, and multiallelic data.
- Define and independently validate an evidence-based indel policy; do not
  silently copy DeepSomatic indels as a consensus innovation.
- Run the frozen policy unchanged on completed HG008 WGS data and a matching
  truth/reference contract when available.

These findings are implementation backlog, not changes made in this report.
