# Separated three-class v2: completed candidate validation

**Execution checks passed; biological training approval did not occur.**

All nine stage outputs preserve their declared established Somatic allele sets exactly.
SNP, indel and aggregate metrics match the corresponding baseline in both regions.
Consumed VCF/reference/benchmark inputs and validated code passed the recorded hash checks.
BAM pilots record source paths and read filters, not newly computed whole-BAM checksums.
No cohort was executed.

See [policy](../../SEPARATED_THREE_CLASS_V2.md) and
[loss diagnosis](../three_class_postfix_20260918/RESCUE_LOSS_DIAGNOSIS.md).

## Somatic results

These restore the established refined policy; they are not a newly tuned improvement over it.
HG008 uses the recommended tumorvariants truth. Comparisons are paired within each domain,
not across different truths or target regions. HG008 has already been inspected and is not
an untouched holdout. Negative-class queries against Somatic truth are collision screens,
not Germline/Reference precision estimates.

### ukb: snp

| Dataset | Stage | TP | FP | FN | Precision | Recall | F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| seqc2_wes_ll | consensus | 1012 | 34 | 1193 | 0.967495 | 0.458957 | 0.622578 |
| seqc2_wes_ll | first | 1021 | 37 | 1184 | 0.965028 | 0.463039 | 0.625804 |
| seqc2_wes_ll | realignment | 1023 | 36 | 1182 | 0.966006 | 0.463946 | 0.626838 |
| seqc2_wgs_il | consensus | 2091 | 8 | 114 | 0.996189 | 0.948299 | 0.971654 |
| seqc2_wgs_il | first | 2091 | 8 | 114 | 0.996189 | 0.948299 | 0.971654 |
| seqc2_wgs_il | realignment | 2091 | 8 | 114 | 0.996189 | 0.948299 | 0.971654 |
| hg008_wgs | consensus | 411 | 15 | 13 | 0.964789 | 0.969340 | 0.967059 |
| hg008_wgs | first | 411 | 20 | 13 | 0.953596 | 0.969340 | 0.961404 |
| hg008_wgs | realignment | 411 | 18 | 13 | 0.958042 | 0.969340 | 0.963658 |

### ukb: indel

| Dataset | Stage | TP | FP | FN | Precision | Recall | F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| seqc2_wes_ll | consensus | 40 | 1 | 55 | 0.975610 | 0.421053 | 0.588235 |
| seqc2_wes_ll | first | 40 | 1 | 55 | 0.975610 | 0.421053 | 0.588235 |
| seqc2_wes_ll | realignment | 40 | 1 | 55 | 0.975610 | 0.421053 | 0.588235 |
| seqc2_wgs_il | consensus | 81 | 6 | 14 | 0.931034 | 0.852632 | 0.890110 |
| seqc2_wgs_il | first | 81 | 6 | 14 | 0.931034 | 0.852632 | 0.890110 |
| seqc2_wgs_il | realignment | 81 | 6 | 14 | 0.931034 | 0.852632 | 0.890110 |
| hg008_wgs | consensus | 92 | 2 | 134 | 0.978723 | 0.407080 | 0.575000 |
| hg008_wgs | first | 92 | 2 | 134 | 0.978723 | 0.407080 | 0.575000 |
| hg008_wgs | realignment | 92 | 2 | 134 | 0.978723 | 0.407080 | 0.575000 |

### ukb: records

| Dataset | Stage | TP | FP | FN | Precision | Recall | F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| seqc2_wes_ll | consensus | 1052 | 35 | 1248 | 0.967801 | 0.457391 | 0.621199 |
| seqc2_wes_ll | first | 1061 | 38 | 1239 | 0.965423 | 0.461304 | 0.624301 |
| seqc2_wes_ll | realignment | 1063 | 37 | 1237 | 0.966364 | 0.462174 | 0.625294 |
| seqc2_wgs_il | consensus | 2172 | 14 | 128 | 0.993596 | 0.944348 | 0.968346 |
| seqc2_wgs_il | first | 2172 | 14 | 128 | 0.993596 | 0.944348 | 0.968346 |
| seqc2_wgs_il | realignment | 2172 | 14 | 128 | 0.993596 | 0.944348 | 0.968346 |
| hg008_wgs | consensus | 503 | 17 | 147 | 0.967308 | 0.773846 | 0.859829 |
| hg008_wgs | first | 503 | 22 | 147 | 0.958095 | 0.773846 | 0.856170 |
| hg008_wgs | realignment | 503 | 20 | 147 | 0.961759 | 0.773846 | 0.857630 |

### medexome: snp

| Dataset | Stage | TP | FP | FN | Precision | Recall | F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| seqc2_wes_ll | consensus | 543 | 17 | 252 | 0.969643 | 0.683019 | 0.801476 |
| seqc2_wes_ll | first | 548 | 19 | 247 | 0.966490 | 0.689308 | 0.804699 |
| seqc2_wes_ll | realignment | 549 | 19 | 246 | 0.966549 | 0.690566 | 0.805576 |
| seqc2_wgs_il | consensus | 675 | 2 | 120 | 0.997046 | 0.849057 | 0.917120 |
| seqc2_wgs_il | first | 675 | 2 | 120 | 0.997046 | 0.849057 | 0.917120 |
| seqc2_wgs_il | realignment | 675 | 2 | 120 | 0.997046 | 0.849057 | 0.917120 |
| hg008_wgs | consensus | 136 | 3 | 17 | 0.978417 | 0.888889 | 0.931507 |
| hg008_wgs | first | 136 | 4 | 17 | 0.971429 | 0.888889 | 0.928328 |
| hg008_wgs | realignment | 136 | 3 | 17 | 0.978417 | 0.888889 | 0.931507 |

### medexome: indel

| Dataset | Stage | TP | FP | FN | Precision | Recall | F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| seqc2_wes_ll | consensus | 21 | 1 | 13 | 0.954545 | 0.617647 | 0.750000 |
| seqc2_wes_ll | first | 21 | 1 | 13 | 0.954545 | 0.617647 | 0.750000 |
| seqc2_wes_ll | realignment | 21 | 1 | 13 | 0.954545 | 0.617647 | 0.750000 |
| seqc2_wgs_il | consensus | 29 | 2 | 5 | 0.935484 | 0.852941 | 0.892308 |
| seqc2_wgs_il | first | 29 | 2 | 5 | 0.935484 | 0.852941 | 0.892308 |
| seqc2_wgs_il | realignment | 29 | 2 | 5 | 0.935484 | 0.852941 | 0.892308 |
| hg008_wgs | consensus | 17 | 0 | 41 | 1.000000 | 0.293103 | 0.453333 |
| hg008_wgs | first | 17 | 0 | 41 | 1.000000 | 0.293103 | 0.453333 |
| hg008_wgs | realignment | 17 | 0 | 41 | 1.000000 | 0.293103 | 0.453333 |

### medexome: records

| Dataset | Stage | TP | FP | FN | Precision | Recall | F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| seqc2_wes_ll | consensus | 564 | 18 | 265 | 0.969072 | 0.680338 | 0.799433 |
| seqc2_wes_ll | first | 569 | 20 | 260 | 0.966044 | 0.686369 | 0.802539 |
| seqc2_wes_ll | realignment | 570 | 20 | 259 | 0.966102 | 0.687575 | 0.803383 |
| seqc2_wgs_il | consensus | 704 | 4 | 125 | 0.994350 | 0.849216 | 0.916070 |
| seqc2_wgs_il | first | 704 | 4 | 125 | 0.994350 | 0.849216 | 0.916070 |
| seqc2_wgs_il | realignment | 704 | 4 | 125 | 0.994350 | 0.849216 | 0.916070 |
| hg008_wgs | consensus | 153 | 3 | 58 | 0.980769 | 0.725118 | 0.833787 |
| hg008_wgs | first | 153 | 4 | 58 | 0.974522 | 0.725118 | 0.831522 |
| hg008_wgs | realignment | 153 | 3 | 58 | 0.980769 | 0.725118 | 0.833787 |

## Provisional paired-read evidence

| Dataset | Stage | Germline supported / usable / selected SNPs | Reference supported / usable / selected SNPs |
| --- | --- | ---: | ---: |
| seqc2_wes_ll | consensus | 57 / 126 / 128 | 0 / 128 / 128 |
| seqc2_wes_ll | first | 57 / 126 / 128 | 0 / 128 / 128 |
| seqc2_wes_ll | realignment | 57 / 126 / 128 | 0 / 128 / 128 |
| seqc2_wgs_il | consensus | 51 / 126 / 128 | 0 / 128 / 128 |
| seqc2_wgs_il | first | 51 / 126 / 128 | 0 / 128 / 128 |
| seqc2_wgs_il | realignment | 51 / 126 / 128 | 0 / 128 / 128 |
| hg008_wgs | consensus | 71 / 128 / 128 | 0 / 128 / 128 |
| hg008_wgs | first | 71 / 128 / 128 | 0 / 128 / 128 |
| hg008_wgs | realignment | 71 / 128 / 128 | 0 / 128 / 128 |

These are deterministic, truth-blind SNP pilot retention rates, not class accuracy.
Usable means valid, nonzero paired-DNA read counts, not adequate depth for class support.
Stages reuse many sites and must not be pooled as independent samples. Unassessed
records remain withheld; negative indels require haplotype-aware validation.
No supported negative in the targeted exact-allele challenge overlapped known Somatic truth.
Rescue stages use the same native nominations and retain every DNA Somatic baseline allele;
their negative sets are therefore subsets of the challenged consensus negatives.
This bounded challenge is not a genome-wide guarantee or a haplotype-equivalence test.

Reference uses zero ALT/other alleles and at least 299 observations in each DNA sample
(1% detection limit, 95% confidence per sample under the idealized independence model).
Zero supported Reference yield does not justify relaxing this approved threshold.
HG008 N-P normal gVCF corroboration is different-tissue evidence, not Germline truth
or paired tumor/normal Reference approval. Every output remains TRAINING_ELIGIBLE=NO.

## Artifacts and reproduction

Heavy outputs, retained outside Git: `/t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/examples/seqc2/comparison/separated_three_class_v2_20260919`.
Each dataset contains `consensus`, `first` and `realignment` directories:
`candidates.vcf.gz` is the separated candidate output; `evidence_gate/candidate.evidence.vcf.gz`
adds pilot evidence status without changing FILTER. Neither is an approved training VCF.
Lightweight evidence, source/output hashes, commands and all metrics are archived here.
Earlier Sept18 attempts are failed/interrupted development artifacts, not canonical results.

Run each validation in a **fresh** destination; existing outputs are never overwritten:

```bash
.venv/bin/python examples/seqc2/scripts/validate_separated_three_class.py --native-validation examples/seqc2/comparison/three_class_postfix_20260918/seqc2_wes_ll/validation.json --samplesheet examples/seqc2/hybrid/csv/seqc2_wes_ll_hybrid.csv --outdir /path/to/fresh_validation/seqc2_wes_ll
.venv/bin/python examples/seqc2/scripts/validate_separated_three_class.py --native-validation examples/seqc2/comparison/three_class_postfix_20260918/seqc2_wgs_il/validation.json --samplesheet examples/seqc2/hybrid/csv/seqc2_wgs_il_hybrid.csv --outdir /path/to/fresh_validation/seqc2_wgs_il
.venv/bin/python examples/seqc2/scripts/validate_separated_three_class.py --native-validation examples/seqc2/comparison/three_class_postfix_20260918/hg008_wgs/validation.json --samplesheet examples/seqc2/hybrid/csv/hg008_wgs_hybrid.csv --outdir /path/to/fresh_validation/hg008_wgs --truth /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/data/giab/data/giab_benchmarks/smvar_v0.3/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz
.venv/bin/python examples/seqc2/scripts/check_negative_truth_collisions.py --vcf /path/to/fresh_validation/seqc2_wes_ll/consensus/candidates.vcf.gz --truth /t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/examples/seqc2/comparison/common_policy_20260914/wgs_il/ukb/benchmark_truth.vcf.gz --samplesheet examples/seqc2/hybrid/csv/seqc2_wes_ll_hybrid.csv --fasta /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta --out /path/to/fresh_validation/seqc2_wes_ll/negative_collision_check.json
.venv/bin/python examples/seqc2/scripts/check_negative_truth_collisions.py --vcf /path/to/fresh_validation/seqc2_wgs_il/consensus/candidates.vcf.gz --truth /t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/examples/seqc2/comparison/common_policy_20260914/wgs_il/ukb/benchmark_truth.vcf.gz --samplesheet examples/seqc2/hybrid/csv/seqc2_wgs_il_hybrid.csv --fasta /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta --out /path/to/fresh_validation/seqc2_wgs_il/negative_collision_check.json
.venv/bin/python examples/seqc2/scripts/check_negative_truth_collisions.py --vcf /path/to/fresh_validation/hg008_wgs/consensus/candidates.vcf.gz --truth /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/data/giab/data/giab_benchmarks/smvar_v0.3/HG008-T_somatic_smvar_benchmark_v0.3_tumorvariants.vcf.gz --samplesheet examples/seqc2/hybrid/csv/hg008_wgs_hybrid.csv --fasta /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta --out /path/to/fresh_validation/hg008_wgs/negative_collision_check.json
```

HG008 orthogonal normal corroboration (repeat for each stage):

```bash
.venv/bin/python examples/seqc2/scripts/check_three_class_normal_gvcf.py --pilot /path/to/fresh_validation/hg008_wgs/consensus/bam_pilot.json --normal-gvcf /path/to/HG008-N-P.GRCh38.deepvariant.g.vcf.gz --out /path/to/fresh_validation/hg008_wgs/consensus/normal_gvcf.json
.venv/bin/python examples/seqc2/scripts/record_separated_three_class.py --root /path/to/fresh_validation --outdir /path/to/fresh_lightweight_archive
```

The normal gVCF source path and SHA256 are in the archived `normal_gvcf.json` files.
Reproduction requires the retained heavy inputs; they are not bundled into Git.
The validation reports also record exact benchmark and evidence-assessment commands.
The cohort wrapper is a separate, candidate-only preparation step; biological approval
and a reviewed cohort pilot remain necessary before any model-training use.
