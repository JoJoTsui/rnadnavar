# Refined-policy cohort pilot review — 2026-09-16

Policy and runner: `c7405f4e`, frozen `seqc2_refined_v2+seqc2_refined_gate_v1`.
The three-sample pilot completed; the remaining cohort was not launched during
this review. No policy thresholds, VCF labels or original inputs were changed.

## Execution and integrity

Shared output root:
`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_refined_native_v2_20260916`.
Each sample's artifacts are in `<sample_id>/attempt001/`; `state.json` and
`completion.json` identify commands, inputs, hashes and output paths.
Run log: the shared `examples/seq2neo/refined_native_v2.pilot.session.log`.

All six consensus/rescue structural audits passed. A fresh read-only comparison
confirmed zero baseline Somatic alleles lost, and SHA-256 of every recorded source
VCF/index and output still matched completion state after report-only QC.
No mapping, variant calling or fresh VEP was run.

| Sample | Execution minutes | DNA Somatic | Rescued Somatic | Added | DNA QC* | Rescue QC |
|---|---:|---:|---:|---:|---|---|
| PRJNA298330_4032 | 9.5 | 397 | 397 | 0 | PASS | WARN |
| PRJNA298376_4255 | 25.8 | 541 | 545 | 4 | PASS | PASS |
| PRJNA298376_4278 | 33.2 | 523 | 528 | 5 | PASS | PASS |

Durations are per-sample wall time, not isolated CPU benchmarks. Two workers ran
concurrently. The nine additions are SNPs, marked `RESCUE_PROMOTED=YES` with
`dna_nominated_rna_supported` rationale. Available gnomAD AF values were below
0.001; three additions lacked GNOMAD_AF. Missing AF is not evidence of rarity.
The T>C addition at chr13:25097966 in 4255 has REDI_CANONICAL=YES but no
REDI_ACCESSION and REDI_EVIDENCE=NONE: a canonical transition alone is not an
editing-database match under the frozen gate.

## Report-only QC and schema limitation

QC root under the shared output root: `pilot_qc.Lhj4wuT0/`.
For each sample, `consensus/` and `rescue/` retain the original default QC reports;
`consensus_dna_scope/` retains corrected DNA-only reports. All calls omitted
`--apply` and `--verify-bam`, used `--threads 1` and an 8 GiB address-space limit.
Thus these are VCF-only checks, not BAM-supported biological verification.

The default QC assumes DNA-prefixed caller names. DNA consensus instead uses
`deepsomatic`, `mutect2`, `strelka`; it falsely marked every sample FAIL at S5
(no DNA participation) and every Somatic call R2. The corrected checks used:

```json
{"variant_rules": {"dna_prefix": ""}}
```

This override was applied **only to known DNA-only consensus outputs**, never to
mixed-modality rescue. It changes identifier interpretation, not thresholds.
The temporary config is `.artifacts/pilot_dna_only_label_qc.json` in the development
repo. Reproduce each corrected report into a new output directory with:

```bash
prlimit --as=8589934592 -- .venv/bin/python bin/label_qc.py \
  --truth-vcf /absolute/sample/attempt001/refined.vcf.gz \
  --sample-id SAMPLE --threads 1 \
  --config .artifacts/pilot_dna_only_label_qc.json --out /absolute/new-qc-directory
```

For rescue use `refined.rescue.vcf.gz` and omit the prefix override. Each sample
was evaluated separately, so cohort-adaptive outlier checks are not a substitute
for final all-cohort QC. Consensus lacks the inherited rescue annotations;
its PASS is not evidence that its baseline calls are population-conflict-free.
Legacy rescue caller fields also do not fully describe every refined baseline
record; adapter provenance and the direct baseline comparison remain necessary.

## Biological caveats and rollout decision

| Sample | Baseline annotation conflicts | Rescue QC high | Rescue QC mid | High+mid fraction |
|---|---:|---:|---:|---:|
| 4032 | 71 | 15 | 52 | 16.88% |
| 4255 | 62 | 11 | 47 | 10.64% |
| 4278 | 46 | 8 | 36 | 8.33% |

Baseline annotation-conflict counts equal the report's R1 common-AF hits here.
4032 exceeds the configured 15% flag-burden warning threshold. These conflicts
were retained by the frozen baseline-preservation contract, not introduced by
the nine rescue additions. No automatic common-AF veto or label cleaning was
applied. PASS at sample level does not eliminate individual high-priority flags.

Recommendation: proceed with the remaining samples **only as isolated candidate
generation**, retaining both DNA consensus and rescue plus per-sample audit/QC
reports. Require explicit acknowledgment of this warning before full launch.
Do not use these as approved training labels, claim improved TP without truth,
or replace existing outputs. Training release still requires biological review,
annotation-coverage decisions (no new VEP CSQ), and cohort-level QC.
