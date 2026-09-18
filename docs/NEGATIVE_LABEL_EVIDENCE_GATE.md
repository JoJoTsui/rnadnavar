# Native negative candidates and training eligibility

Status: opt-in evidence assessment, not a training-approved label policy.
The 2026-09-18 validation found Reference candidates overlapping low-VAF somatic
truth and a HG008 indel with mixed germline/somatic context. Candidate nomination
must therefore remain separate from biological training approval.

Current evidence implementation: `negative_evidence_gate_v2`, compatible with
[separated three-class v2](SEPARATED_THREE_CLASS_V2.md). The original combined
v1 rescue remains a diagnostic experiment, not the recommended cohort policy.

## Approved Reference detection-limit requirement

The user selected a 1% allele-fraction detection limit at 95% confidence.
With zero ALT observations, the idealized one-sided binomial upper bound is
`1 - (1 - confidence)^(1/n)`. To make that bound <=0.01 requires
`ceil(log(0.05)/log(0.99)) = 299` observations **in each DNA sample**.

The gate requires zero target ALT and zero other-allele observations, and at
least 299 usable observations in both normal and tumor. Depth is never summed
between samples or callers. Evidence must have MAPQ/BQ >=20 and BAQ enabled.
The BAM pilot explicitly excludes unmapped, secondary, QC-failed, duplicate
and supplementary alignments (`flag_filter=0xF04`), suppresses overlapping mates
and ignores orphan/improper pairs. Reports missing these settings are rejected.
Nevertheless, read independence is an approximation. This statistical
bound does not exclude mapping bias, systematic error, CNV or correlated reads.
Confidence is per sample, not a claimed joint two-sample 95% guarantee.

Lower-depth/ALT-observed sites remain Reference **candidates**, but evidence
support is withheld. Germline/Reference FILTERs are not silently reassigned to
Somatic, PASS, or another training class.

## Evidence gate

`bin/assess_negative_label_evidence.py` consumes a biological candidate VCF and
a paired-BAM report from `pilot_three_class_bam_evidence.py`. The report is
bound to the exact VCF SHA256. Duplicate allele evidence, class mismatches,
insufficient evidence-generation quality settings and invalid counts fail closed.

```text
Native three-class candidate VCF + hash-bound paired DNA evidence
    |
    + missing evidence / unsupported indel ----> WITHHELD
    + conflicting normal evidence ------------> CONFLICT
    + unresolved allele / haplotype context ---> WITHHELD
    + Reference meets 1% zero-ALT depth bound --> SUPPORTED (evidence only)
    + Germline meets paired-read checks ------> SUPPORTED (evidence only)
    |
    + all records: unchanged FILTER, TRAINING_ELIGIBLE=NO
```

Germline SNP read corroboration requires normal depth >=20, ALT >=5, AF >=0.2,
tumor depth >=20 and tumor ALT >=3. A >0.3 normal/tumor AF shift is conservatively
withheld for context review: it can reflect purity/CNV/LOH or somatic-on-germline,
not necessarily an error. These are provisional review thresholds, not calibrated
Germline precision. Non-SNPs are withheld pending haplotype-aware validation;
the tool does not claim that all indels are false or irrelevant.

Statuses are annotations, not new biological FILTER labels. The output records
the reason and actual normal/tumor ref:alt:other:depth counts when measured.
Missing counts remain missing. No training approval option exists in this tool.
It creates no training manifest. **Do not give its candidate VCF directly to
FILTER-only training consumers**; they do not enforce eligibility annotations.
The cohort candidate Parquet exporter continues to mark all rows unapproved.

The current evidence interface is a bounded pilot JSON report. Unmeasured
records are withheld, not extrapolated from the pilot. This is not a substitute
for collecting cohort-wide evidence before training-dataset generation.

## Usage and verification

Use new output directories; neither source VCFs nor BAMs are modified:

```bash
.venv/bin/python bin/assess_negative_label_evidence.py --vcf /absolute/path/candidate.vcf.gz --bam-evidence /absolute/path/bam_pilot.json --outdir /absolute/path/new_evidence_gate --max-reference-af 0.01 --confidence 0.95
```

Outputs are `candidate.evidence.vcf.gz` and `report.json`. Failed runs leave a
clearly named partial VCF and failed report, not a published candidate output.
Repeated evidence gating is rejected. An existing Number=1/String eligibility
field is accepted only when every value is absent or NO; an existing positive
eligibility decision is never overwritten or silently converted.

Current separated-policy validation uses `validate_separated_three_class.py`
for consensus and both rescue rounds, against each stage's established Somatic
baseline. See its policy document for the frozen gates and fresh-output usage.

The following commands reproduce the **historical v1 experiment**, which has
known Somatic-retention failures. They must not be used as cohort approval.
Native nomination validation uses `validate_three_class_policy.py`; the v1
rescue experiment is handled by `validate_three_class_rescue.py`:

```bash
.venv/bin/python examples/seqc2/scripts/validate_three_class_rescue.py --consensus-validation /absolute/path/completed_consensus/validation.json --outdir /absolute/path/new_rescue_validation --execute
```

Without `--execute`, it prepares and verifies a plan only. Use a different output
directory for execution. It rejects stale classifier provenance and checks the
consensus checksum. HG008 should specify `--truth` with the recommended
`HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz`. First and realignment
rounds use their own RNA panels and annotated rescue sources. UKB and MedExome
metrics stay separate. Germline/Reference matches to somatic truth remain
collision screens, never estimates of negative-class precision.

These commands do not launch Nextflow, alignment, calling or the 66-sample cohort.
Production workflow defaults and variant-calling task/cache definitions are unchanged.
