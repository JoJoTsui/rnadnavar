# Three-class validation checkpoint — 2026-09-18

**Decision: do not start the 66-sample regeneration or approve training labels.**
The native negative rules remain provisional. This validation found both a
Somatic admission regression (fixed separately) and unresolved biological
limitations. No alignment, calling, source output or cache was modified.

## What ran

- Full-input standalone consensus on SEQC2 WES LL, SEQC2 WGS IL and HG008 WGS,
  using the normalized caller inputs identified in the frozen manuscript bundle.
- Somatic, Germline and Reference queries separately compared with somatic
  truth, in both UKB and MedExome target regions plus each dataset's HC BED.
- HG008 queries and both baselines also compared with the release-recommended
  single-sample tumorvariants truth. Historical multi-sample results are retained.
- Truth-blind, deterministic paired-BAM spot checks of 32 Germline and 32
  Reference SNPs per dataset, using MAPQ/BQ >=20 and BAQ. No indel accuracy claim
  is made from these SNP pileups.
- The same HG008 pilot sites checked against the public PacBio normal gVCF.
  That sample is HG008-N-P, whereas the Illumina workflow normal is HG008-N-D.
  This is cross-tissue corroboration, not certified genotype truth.

Initial full screens used policy code at `664aff40` (policy introduced in
`d9c34e13`). The indel fallback correction is `50705ab2`. The initial metrics
below are **pre-fix**, not a completed full benchmark of the corrected code.
Actual affected loci were replayed through the corrected CLI separately.
All initial screens reported unchanged source/code hashes.

## Negative labels overlapping somatic truth

These are matched known-somatic events, **not Germline/Reference precision or
FP counts**. Raw som.py negative-query statistics use a Somatic truth set and
must not be interpreted as negative-class accuracy.

| Dataset | Germline UKB / MedExome | Reference UKB / MedExome |
|---|---:|---:|
| SEQC2 WES LL | 0 / 0 | 4 / 3 |
| SEQC2 WGS IL | 0 / 0 | 10 / 3 |
| HG008 WGS | 1 / 0 | 0 / 0 |

All SEQC2 Reference overlaps are SNPs. The WES four have two DeepSomatic tumor
ALT reads among 41–47 observed allele reads. They pass the provisional native
Reference gate, but that evidence does not exclude a low-VAF somatic allele.
Zero observed collisions elsewhere does not establish class precision or
exhaustiveness. UKB and MedExome are overlapping domains; do not add their counts.

The HG008 overlap is an indel, chr10:17796654 G>GT. DeepSomatic calls GERMLINE;
Strelka has normal TAR/TIR=63/57 and tumor TAR/TIR=1/86. The independent normal
gVCF has GT=0/1, GQ=42, DP=31, AD=16,15,0. The full assembly-context VCF has
normal 1|0 and tumor 1|1, while its somatic-only representation uses normal 0|0
and tumor 0|1. Thus a shared germline allele and an additional somatic event can
map to the same GRCh38 allele. Do not declare this a simple Germline FP or use
it to tune thresholds. It needs haplotype-aware adjudication and an explicit
training-label ambiguity policy (the existing NoConsensus vocabulary permits
abstention without inventing another training class).

The [official HG008 release README](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/NIST_HG008-T_somatic-smvar_DraftBenchmark_V0.3-20260425/README.md)
recommends tumorvariants for ordinary somatic benchmarking, warns about
somatic/germline representation interactions, and explicitly says the full
tumor/normal file's germline genotypes are unvalidated and unsuitable as a
germline benchmark. We preserve it as contextual evidence only.

## Paired-read and orthogonal pilot

| Dataset | Germline normal-read corroborated / sampled | Reference paired-zero-ALT corroborated / sampled |
|---|---:|---:|
| SEQC2 WES LL | 30 / 32 | 0 / 32 |
| SEQC2 WGS IL | 32 / 32 | 1 / 32 |
| HG008 WGS | 32 / 32 | 1 / 32 |

All remaining pilot outcomes were inconclusive, not proven errors. Germline
corroboration required normal depth >=20, ALT >=5 and AF >=0.2. Reference
corroboration required >=60 reads and zero ALT in both DNA samples. The latter
corresponds only to an idealized approximately 5% one-sided detection limit,
not proof of no low-frequency mutation. Counts from shared-read callers are
not independent experiments.

For HG008's orthogonal normal gVCF, 31/32 Germline SNPs had qualifying normal ALT
evidence; 24/32 Reference SNPs had qualifying normal reference evidence. The
rest were inconclusive. Normal-only support does not establish tumor Reference.
The gVCF was downloaded to ignored storage, never over the original inputs:
`https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/PacBio_Revio_20240125/pacbio-wgs-wdl_germline_20240206/HG008-N-P.GRCh38.deepvariant.g.vcf.gz`.

Some existing indexes warn that their timestamps predate their data files.
These warnings are retained in logs; no source index was rewritten. Small
pilot checks are not exhaustive index or BAM validity certification.

## Somatic regression isolated and fixed

The opt-in branch evaluated native Somatic admission but inadvertently skipped
the existing threshold fallback for unadmitted indels. A minimal two-caller
fixture reproduced Somatic -> NoConsensus in 0.13 seconds. The fix asks the
existing classifier for Somatic admission first, then applies the separate
negative nomination/conflict logic. It never imports the legacy classifier's
negative labels as negative evidence. Germline/Reference thresholds were not
tuned against HG008 or changed by this fix.

| Actual-input replay | Initial full-input Somatic label losses | Restored by fix | Remaining explicit conflicts |
|---|---:|---:|---:|
| WES LL | 64 | 64 | 0 |
| WGS IL | 200 | 197 | 3 |
| HG008 | 33 | 30 | 3 |

These are all input loci, not only benchmark-region TP/FN. No replay loci were
missing and caller source hashes stayed unchanged. The residual loci and their
rule traces are in `evidence/*/fallback_fix_replay/replay.json`.
A fresh full post-fix benchmark remains pending; do not use the pre-fix metrics
as a claimed performance result for the corrected implementation.

Verification: 294 focused tests passed (66 existing NumPy deprecation warnings).
The debugging skill's failing-test/replay loop identified and verified the fix.

## Files and reproduction

`summary.json` indexes results, source locations and archived evidence checksums.
`evidence/` contains only small JSON reports, exact conflict loci, metrics and
commands. Heavy VCFs/logs remain under the ignored directory
`examples/seqc2/comparison/three_class_validation_20260918/`; the orthogonal gVCF
is under `.artifacts/three_class_validation_20260918/orthogonal/`.

Example new-output rerun, using current code (not a replay of pre-fix metrics):

```bash
.venv/bin/python examples/seqc2/scripts/validate_three_class_policy.py --frozen-validation docs/manuscript/frozen_native_gate_20260916/evidence/seqc2_wes_ll/validation.json --outdir examples/seqc2/comparison/three_class_postfix_wes_new
```

Repeat with seqc2_wgs_il and hg008_wgs and distinct output directories. For
HG008, subsequently use benchmark_three_class_truth_sensitivity.py with
--validation pointing to the completed new report and --truth pointing to the
release's tumorvariants VCF. Exact initial commands are archived in each report.
The scripts refuse overwriting their result directories/files.

## Next decision and remaining validation

1. Keep candidate nomination separate from training eligibility. Native-only
   Reference is not currently a high-confidence training label; neither a
   missing normal measurement nor a failed somatic call establishes Reference.
2. Define the Reference detection limit and acceptable uncertainty, then
   validate a paired-evidence eligibility gate on SEQC2. Keep HG008 held for
   evaluation, not threshold fitting.
3. Adjudicate indels and somatic-on-germline/LOH contexts. The present SNP pilot
   cannot approve indel labels. Do not use somatic-truth absence as negative truth.
4. Run the full corrected-policy comparison and both rescue-round validations;
   the current empirical screen covers DNA consensus only, not new rescue VCFs.
5. Only after these gates, pilot representative cohort samples and regenerate
   the 66-sample consensus/rescue/annotation/Parquet in a new output namespace.
   No cohort VCF regeneration was started in this validation.
