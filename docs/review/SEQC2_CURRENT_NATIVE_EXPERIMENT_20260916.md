# Current native policy: SEQC2 reproduction and bounded optimization

The subsequent [indel refinement and rescue audit](SEQC2_INDEL_RESCUE_FOLLOWUP_20260916.md)
records additional TP gains without added FP relative to this fixed candidate.
The values below remain the preserved first-round experiment.

Status: completed standalone development experiment, not promoted to workflow
defaults. HG008 was not read or used for selection. All source VCF checksums
were unchanged after the baseline run. No mapping, variant calling, or Nextflow
process was launched; original workflow outputs and caches were untouched.

## What the reproduction corrected

The historical winning replay is not the current implementation. Its indels
were DeepSomatic-derived and its WGS candidate scope was restricted. Here the
current consensus driver reads all verified DNA caller records; truth and
target intervals apply during scoring. Current native SNP rules plus ordinary
2-of-3 indels do not reproduce historical gains.

In WGS UKB, current SNP consensus adds 9 TP and 12 FP versus DeepSomatic.
All 12 added FP have Mutect2 artifact filters missed by the two exact-string
vetoes, including normal_artifact, strand_bias, haplotype, and other combinations.
WES loses 6 indel TP and removes 3 indel FP; the six lost TP have Mutect2
contamination;weak_evidence calls with two alternate reads. Simply requiring
more alternate reads or higher TLOD would also remove real variants.

## Shared candidate rules

The candidate starts with current native SNP / threshold-indel consensus.
It applies identical rules in both datasets and both target domains:

```text
Current native SNP classified Somatic
  DeepSomatic native PASS/.? -------------------------- retain
  otherwise Mutect2 filters only PASS/./contamination/weak_evidence?
    yes ---------------------------------------------- retain
    no or missing ------------------------------------ remove

Current threshold-consensus indels -------------------- retain
Additional biallelic indel
  DeepSomatic PASS/.; finite QUAL >=20; tumor ALT count >=3
  Mutect2 tumor ALT count >=2
  Mutect2 normal DP >=10 AND normal ALT count =0
  Mutect2 finite GERMQ >=20 AND finite TLOD >0
  Mutect2 filters only PASS/./contamination/weak_evidence
    all satisfied ------------------------------------ add
    any failed/missing ------------------------------- reject addition
```

The original native SNP prerequisites (DeepSomatic positive QUAL and Mutect2
TLOD >=12 / GERMQ >=60 for non-backbone additions) remain in force because
this assay only prunes already admitted SNPs. DeepSomatic PASS/. is not itself
a new override for other rejected consensus records.

QUAL and GERMQ are caller confidence scores; TLOD measures tumor alternate
evidence. The thresholds are explicit development choices, not calibrated
ensemble probabilities. DP/AD refer to the resolved sample in the source VCF.
The normal requirement supplies evidence against germline/normal artifacts;
absence of a normal record is not evidence of absence. Shared-read caller
counts are not summed or treated as independent measurements.

The SNP rule removes all 12 added WGS FP and one added TP (map_qual), retaining
eight added TP versus DeepSomatic. It retains all five added WES TP. Indel
corroboration adds 2 TP/0 FP to current WES UKB and 7 TP/2 FP to current WGS
UKB. Thus the indel change is not FP-free relative to current consensus.

## Aggregate comparison

Each entry is TP/FP/FN and F1. All queries use som.py -N, the same merged
SEQC2 truth, HC -R, target -T and assembly38 FASTA. No -P is used. Production
FILTER labels are preserved; PASS conversion occurs only in benchmark copies.

| Dataset / target | DeepSomatic | Current native code | Shared candidate | Candidate + frozen rescue replay |
|---|---|---|---|---|
| WES / UKB | 1048/38/1252; 0.619019 | 1047/35/1253; 0.619160 | 1049/35/1251; 0.619976 | 1060/37/1240; 0.624080 |
| WES / MedExome | 563/19/266; 0.798016 | 562/18/267; 0.797729 | 563/18/266; 0.798582 | 569/20/260; 0.802539 |
| WGS / UKB | 2168/19/132; 0.966347 | 2165/24/135; 0.964580 | 2171/14/129; 0.968116 | 2171/14/129; 0.968116 |
| WGS / MedExome | 702/6/127; 0.913468 | 700/5/129; 0.912647 | 704/4/125; 0.916070 | 704/4/125; 0.916070 |

| Dataset / target | DeepSomatic precision / recall | Candidate precision / recall | Candidate + rescue precision / recall |
|---|---|---|---|
| WES / UKB | 0.965009 / 0.455652 | 0.967712 / 0.456087 | 0.966272 / 0.460870 |
| WES / MedExome | 0.967354 / 0.679131 | 0.969019 / 0.679131 | 0.966044 / 0.686369 |
| WGS / UKB | 0.991312 / 0.942609 | 0.993593 / 0.943913 | 0.993593 / 0.943913 |
| WGS / MedExome | 0.991525 / 0.846803 | 0.994350 / 0.849216 | 0.994350 / 0.849216 |

The shared consensus candidate improves aggregate F1 and precision in all
four cells. Rescue replay improves WES recall/F1 but WES MedExome precision
falls below DeepSomatic. It therefore fails a strict per-cell precision
non-regression requirement, despite higher F1.

## SNP and indel separation

| Dataset / target | SNP F1: DeepSomatic / candidate / candidate+rescue | Indel F1: DeepSomatic / candidate |
|---|---|---|
| WES / UKB | 0.620456 / 0.622578 / 0.626838 | 0.585714 / 0.556391 |
| WES / MedExome | 0.798817 / 0.801476 / 0.805576 | 0.779661 / 0.727273 |
| WGS / UKB | 0.969739 / 0.971654 / 0.971654 | 0.890052 / 0.883978 |
| WGS / MedExome | 0.916383 / 0.917120 / 0.917120 | 0.848485 / 0.892308 |

Indel-only F1 still trails DeepSomatic in three cells. This is not a claim of
SOTA or superiority for every variant type. WES rescue-replay F1 also remains
below the historical DeepSomatic-indel replay (0.624706 UKB, 0.803949 MedExome).
No old result is deleted or relabeled as the new method.

## Artifacts and reproduction

Root: `examples/seqc2/comparison/current_native_audit_20260916_v2/`.

- `audit.json`: eight baseline comparisons, exact CLI, code/source hashes,
  complete status, unchanged-source verification.
- `{wes_ll,wgs_il}/current_native.vcf.gz`: full annotated current-code output.
- `wes_indel_attribution/` and `wgs_attribution/`: normalized caller evidence,
  benchmark scratch, lost/added allele attribution. Benchmark parity verified.
- `{wes,wgs}_combined_candidate/evaluation.json`: fixed-rule candidate metrics,
  selected alleles, removed SNPs, source/script hashes and commands.
- `{wes,wgs}_candidate_gate_replay/evaluation.json`: frozen gate transfer metrics.
- Under each candidate/replay directory, `{ukb,medexome}/query.vcf.gz` is the
  benchmark-only allele query. These are not annotated training labels.
- `wes_indel_candidate/` retains the earlier indel-only experiment.

The failed first root `current_native_audit_20260916/` is explicitly marked
failed_input_staging, with no metrics. The successful v2 run stages sample-prefixed
VCF/index symlinks because consensus caller recognition expects that naming.

Checksums:

| Artifact under root | SHA-256 |
|---|---|
| audit.json | 7fc99782c808112c462d91bad99cf54738756bb79b6d8fe624141775e8d9587e |
| wes_combined_candidate/evaluation.json | 10a7c49462b8921ac285d58a5da1769af349a59dc9eace964bf1d36ccade818c |
| wgs_combined_candidate/evaluation.json | 8a9d998c1f8015b1eea2794ad33b16c1573bd7ade241e0f52ce427fa16dbb1d9 |
| wes_candidate_gate_replay/evaluation.json | b7f67933bd27a0899554586ebf3586cd1d9dfd50ba62431241df63cb79c904d7 |
| wgs_candidate_gate_replay/evaluation.json | fec8e85cda1a6b1e5ef1553cbd13e1f9087a88a8bd551e9ebb4229216b859ba9 |

The standalone helpers are `audit_current_native_policy.py`,
`attribute_native_indel_errors.py`, `test_corroborated_indel_policy.py`
(`--strict-snv-evidence` for the combined candidate), and
`benchmark_frozen_gate_transfer.py`, under `examples/seqc2/scripts/`.
Each requires a fresh output directory; exact executed commands are preserved
in the JSON records. Eight focused tests and Python compilation passed.

## Decision and next boundary

Keep this candidate as a SEQC2 development result. Do not change workflow
defaults yet. The rescue experiment transfers the preserved *realignment*
gate additions, inheriting historical candidate scope and annotations; it
does not rerun either rescue driver or validate first-round rescue portability.

Next investigate the remaining indel recall gap and the two WES rescue FP
using existing caller/normal evidence. Preserve the aggregate-winning candidate
as the comparator. If a rule is selected, reproduce its allele set with the
actual consensus/rescue drivers before freezing it for HG008 evaluation.
