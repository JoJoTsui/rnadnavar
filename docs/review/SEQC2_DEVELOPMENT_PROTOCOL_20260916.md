# SEQC2 development and subsequent HG008 evaluation

User-approved split, 2026-09-16: select rules using SEQC2 WES-LL and WGS-IL;
evaluate frozen candidates on HG008. Preserve UKB and MedExome comparisons.
Earlier HG008-informed optimization prevents calling HG008 untouched.
Subsequent HG008 failures are evaluation outcomes, not tuning inputs.

## Baseline audit

The winning historical SEQC2 replay retained DeepSomatic-derived indels.
The combined-policy and comparison summaries incorrectly described these as
threshold-consensus indels; those descriptions are now corrected.

Read-only exact allele checks of `examples/seqc2/verified/20260914/` found
45 WES and 96 WGS native length-changing alleles, all present in source
DeepSomatic PASS/. records. This is a membership check, not normalized scoring.
The separate threshold-indel experiment lost 6/4/12/3 TP versus DeepSomatic
in WES-UKB/WES-MedExome/WGS-UKB/WGS-MedExome respectively.
Historical aggregate gains therefore do not validate threshold-indel consensus.

## Rejected rescue hypothesis

Evidence: `examples/seqc2/comparison/rescue_fp_investigation_20260914/evidence.json`,
restricted to scored WES additions retained by the `nomination_biological`
selection in `gate_tests/evaluation.json`.

Requiring positive RNA realigned Mutect2 alternate counts in both F1R2 and
F2R1 would remove the remaining FP chr17:76353649 C>T and chr19:33206715 T>C.
However, nine of the eleven retained TP also have only one alternate
orientation; one has both and one has no Mutect2 record. The requirement loses
at least nine TP, or ten if missing evidence fails closed, to remove two FP.
Reject this rule. These are attributed allele counts from existing evidence,
not a fresh som.py run. Read orientation and forward/reverse strand are
different measurements and must not be conflated.

## Next experiments

1. Reproduce one current consensus policy on both SEQC2 datasets from verified
   callers, recording flags, hashes and allele differences from historical replay.
2. Attribute FN to absent caller candidates, representation changes, consensus
   rejection and rescue rejection. Missing RNA VCF evidence does not establish
   absent RNA coverage; that remains unknown without coverage measurement.
3. Inspect indel candidates rejected by voting and remaining rescue FP using
   caller-specific tumor/normal evidence; test shared rules in all four cells.
4. Preserve per-type TP/FP/FN, precision, recall and F1, including failed
   experiments. Freeze code and arguments before subsequent HG008 evaluation.

Historical WGS replay has a restricted candidate universe. Portable experiments
must select without truth membership and apply HC/target intervals at scoring.
Original inputs, workflow outputs and caller caches remain unchanged.
The initial audit changed documentation only. The subsequent standalone
reproduction and candidate experiments are recorded in
[the current-policy experiment](SEQC2_CURRENT_NATIVE_EXPERIMENT_20260916.md).
Production defaults remain unchanged.

The next [indel refinement and rescue-template audit](SEQC2_INDEL_RESCUE_FOLLOWUP_20260916.md)
preserves rejected hypotheses and identifies the selected development comparator.
