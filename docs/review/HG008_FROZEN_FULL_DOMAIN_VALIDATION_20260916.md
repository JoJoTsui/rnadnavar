# Frozen-policy HG008 validation

User approved evaluation after completion of the HG008 hybrid workflow on
2026-09-16. The shared workflow log reports successful completion; both rescue
rounds exist, including the final realignment-rescue VEP VCF. No new mapping,
variant calling or Nextflow execution is part of this validation.

## Frozen experiment

The SEQC2-selected consensus implementation is unchanged from
`driver_parity_v1/parity.json`: every recorded consensus-driver/library hash
matches. Neither consensus nor rescue thresholds are tuned on HG008.
HG008 was used in older work, so this is a subsequent frozen evaluation, not
a claim that HG008 is an untouched independent holdout.

The new standalone driver processes the entire supplied DNA caller files,
not truth-membership or HC/target-selected subsets. This removes the earlier
region-scoped DNA-baseline limitation. It does not imply genome-wide caller
coverage: the completed workflow may itself have called only configured
intervals. Normalization checks the supplied reference and fails on mismatch;
exact duplicate removal matches the SEQC2 development preparation. Native
sample evidence is retained; multiallelic records are not artificially split
to manufacture biallelic evidence.

The experimental adapter is applied separately to existing annotated first
rescue and final realignment rescue using the matching RNA caller panels.
These remain post-annotation reconstructions from historical outputs, not a
fresh annotation or full workflow release. Normalized DNA and historical rescue
representations, annotation coverage, negative labels and population conflicts
still require label-level review before training use.

## Reproducible inputs and command

Explicit source paths are in
`examples/seqc2/hybrid/frozen_hg008_validation.json`. The manifest intentionally
records this local evaluation's absolute paths: shared workflow inputs and a
local MedExome BED. For another host, copy the manifest and adjust paths before
freezing; do not edit an active run's manifest.

```bash
.venv/bin/python examples/seqc2/scripts/validate_frozen_hybrid_policy.py \
  --manifest examples/seqc2/hybrid/frozen_hg008_validation.json \
  --outdir examples/seqc2/comparison/current_native_audit_20260916_v2/hg008_frozen_full_domain_v1
```

Every attempt requires a new directory. `validation.json` records commands,
code hashes, all source hashes (including FASTA), status, output hashes, and
before/after integrity checks. Failed attempts retain their error and files.
Original files and indexes are never rewritten; normalization and indexing
write only to the new output directory.

Comparators: refined consensus, DNA DeepSomatic, optimized first rescue,
optimized realignment rescue, workflow first rescue and workflow realignment
rescue. Each is scored against HG008 v0.3 truth with its matching HC BED (`-R`)
and separately UKB and MedExome targets (`-T`). SNP, indel and aggregate
TP/FP/FN, precision, recall and F1 are retained. Benchmark copies are sampleless
PASS VCFs; original biological FILTER labels are unchanged. No `-P` is used.

## Status and acceptance

The new driver passed 28 focused tests with existing NumPy deprecation warnings.
The subsequent combined regression run passed 264 tests (62 existing NumPy
deprecation warnings), including both new validation/audit test modules.
The final pre-commit check also included the indel-assay and follow-up-assay
tests: 271 tests passed with the same 62 deprecation warnings.
All 16 configured input paths passed preflight. All 12 benchmark cells completed
under `hg008_frozen_full_domain_v1`, with input, code and manifest integrity
checks passing. The [complete results](HG008_FROZEN_POLICY_RESULTS_20260916.md)
include all variant types and both rescue rounds. Optimized consensus improves
UKB aggregate F1 but narrowly loses to DeepSomatic on MedExome; rescue does not
add net TP on HG008. Completion is not acceptance or default-policy promotion.
Preserve inferior metrics as evaluation results; do not revise thresholds to
fit HG008. Workflow defaults and cohort training labels remain unchanged.

## Structural label audit

`examples/seqc2/scripts/audit_refined_label_contract.py` independently scans
complete output VCFs without loading all records into memory. It checks the
eight-column biological-label contract, FILTER/UNIFIED_FILTER consistency,
nonempty decision rationale, and (for adapter outputs) rationale class,
decoded source-snapshot allele identity, baseline/promotion flag consistency,
and distinct caller support for SNP promotions. It hashes the input before
and after and refuses to overwrite its report. This is not a substitute for
biological label QC, annotation validation, or representation-equivalence
checking; `structural_pass_not_training_approval` explicitly preserves that
distinction.

The consensus structural audit passed all 724,384 records with zero reported
issues and unchanged source checksum. First rescue (2,076,392 records) and
realignment rescue (974,609 records) also passed with zero structural issues
and unchanged source checksums. Their reports are named in the results document.
All three structural audits are complete; biological training-label approval
remains outstanding.

```bash
.venv/bin/python examples/seqc2/scripts/audit_refined_label_contract.py \
  --vcf /absolute/path/to/refined.rescue.vcf.gz \
  --report /new/path/to/structural_audit.json
```
