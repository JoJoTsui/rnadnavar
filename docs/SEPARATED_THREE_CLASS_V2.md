# Separated three-class candidate policy v2

Frozen before revalidation on 2026-09-18. Opt-in; not a training-approved policy.
Implementation: `bin/apply_three_class_labels.py` and negative-evidence gate v2.

Completed execution validation (2026-09-19): all three datasets, consensus and
both rescue rounds passed exact Somatic membership, UKB/MedExome SNP/indel/record
metric parity, structural and recorded integrity checks. See the
[full comparison and evidence archive](validation/separated_three_class_v2_20260919/README.md).
The 128-SNP-per-class pilots supported 57/126 usable Germline sites in WES,
51/126 in SEQC2 WGS and 71/128 in HG008, with zero supported Reference sites.
These are evidence-retention counts, not negative-class accuracy. Candidate-only
cohort preparation is permitted; biological training approval remains open.

## Diagnosis

The first three-class experiment reused old rescue biological labels and
population-AF annotations as new vetoes against an established Somatic baseline.
It also required confirmed DNA verification for every new rescue addition,
whereas the frozen benchmark policy treated missing verification differently.
This changed Somatic admission while ostensibly adding negative classes.

Both rescue rounds have completed for all three datasets. WES UKB loses 50/40
TPs (first/realignment), WGS and HG008 also lose TPs. WES allele-level replay
attributes 46/36 losses to inherited-label-only vetoes and four to population
AF vetoes, removing 11 FPs in each round. Population frequency is not this
person's genotype; failed Somatic classification is not independent negative
truth. The v1 runs are retained as counterexamples, never overwritten.

## Two separate decisions

```text
Existing caller results
  |
  + Established refined Somatic policy, independently per output stage
  |    + DNA consensus baseline
  |    + first-rescue baseline
  |    + realignment-rescue baseline
  |
  + Native DNA negative nominations (native_three_class_v1)
       + explicit DeepSomatic GERMLINE / RefCall and native confidence
       + compatible paired caller read evidence; missing stays missing
  |
  v
Separated union, for each exact allele
  + in selected stage's Somatic baseline --> Somatic candidate, unchanged
  + otherwise native Germline nominated --> Germline candidate
  + otherwise native Reference nominated -> Reference candidate
  + otherwise ---------------------------> NoConsensus
  |
  + annotation/class conflicts --> explicit review INFO, not silent relabeling
  |
  v
Independent negative-evidence assessment and biological approval
  + supported / withheld / conflict annotations
  + all records remain TRAINING_ELIGIBLE=NO
```

There is no truth-dependent inclusion, caller shortcut, score fitting or new
Somatic admission. The exact Somatic allele set must equal the declared,
hash-bound baseline at each stage, including its existing indel decisions and
rescue additions. A native-negative conflict is recorded for adjudication,
not silently called truth. This preserves benchmark membership, **not** an
assurance that every retained Somatic candidate is safe for training.

The established baseline is `seqc2_refined_v2` for DNA and
`seqc2_refined_gate_v1` for rescue. Run those algorithms from existing caller
VCFs, or use their checksum-verified existing results. Do not substitute an
inferior workflow output or the v1 three-class rescue as the baseline.

## Negative nomination versus evidence support

Native negative nomination remains independent of Somatic PASS failures:
DeepSomatic native FILTER, GQ >=30, DP and allele-read total >=20, consistent
class PL when supplied; Germline ALT >=3; Reference ALT <=2 and AF <=0.05.
GQ is model confidence, not a normal-sample genotype confidence. DeepSomatic
GERMLINE GT=0/0 recoding must not be interpreted as Reference. Mutect2/Strelka
can corroborate or contradict read evidence but cannot nominate a negative
class merely by failing Somatic admission. See the v1 review for full nomination
details. Native-only candidates are not high-confidence training labels.

Separate paired-BAM evidence support requires MAPQ/BQ >=20, BAQ, nonduplicate
primary alignments, proper-pair filtering and overlap suppression. Exclude
unmapped, secondary, QC-failed, duplicate **and supplementary** reads (`0xF04`).
The evidence report must explicitly attest these settings; old reports without
that provenance are rejected, not silently reused.

- Reference SNP: zero ALT and other-allele counts, at least 299 observations in
  **each DNA sample**. This is the approved 1% detection limit at 95% confidence
  per sample, under an idealized independent-observation model. Do not sum
  callers or samples, or relax the limit to improve retention.
- Germline SNP: normal depth >=20, ALT >=5, AF >=0.2; tumor depth >=20 and ALT
  >=3. Other alleles, a normal/tumor AF difference >0.3, or unresolved context
  withhold support. Normal depth >=60 with zero ALT is a conflict. These are
  provisional corroboration thresholds, not calibrated Germline precision.
- Indels: keep native nomination provenance but withhold negative evidence
  approval pending haplotype-aware validation. Somatic indel behavior is unchanged.
- Missing/low-depth evidence: withhold, never invent Reference evidence.

Even SUPPORTED is not training approval. Somatic-truth absence cannot prove
Germline/Reference. A zero supported-Reference yield is a coverage limitation,
not a license to lower the threshold or label all non-Somatic positions negative.

## Validation gates frozen before execution

1. Preserve exact Somatic membership for consensus and both rescue rounds on
   all three datasets. Require matching SNP, indel and aggregate som.py metrics
   in UKB and MedExome against the same truth/reference/domain on both sides.
   HG008 uses the recommended tumorvariants truth, not the historical variant.
2. All negative outputs must be traceable to native DNA nomination, not legacy
   rescue labels. Preserve input hashes and avoid mapping/calling entirely.
3. Repeat truth-blind, deterministic negative SNP pilots with corrected read
   filters; report support, abstention and Somatic-truth collisions separately.
   Use available independent normal evidence as corroboration, not germline truth.
4. Test VCF/export eligibility and missingness. No FILTER-only training manifest
   may be produced. Candidate-generation readiness is distinct from biological
   validation, which cannot be claimed without suitable class truth/context.
5. Only after candidate-generation gates pass may a candidate-only 66-sample
   wrapper be prepared, using new shared output paths. Do not automatically
   execute the cohort or enable production workflow defaults. Biological approval
   remains required for model-training labels.

All original callsets, configs and Nextflow calling/cache definitions are unchanged.
