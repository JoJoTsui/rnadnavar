# Frozen optimized consensus/rescue benchmark snapshot

This is the manuscript-facing numerical and provenance archive for the policy
checkpoint at `a7f9ba5a05ffd9abe014e30d0c0379023828bef7`. It does not promote
the policy to workflow defaults or approve generated VCFs as training labels.
Future changes require a new versioned snapshot; do not replace this comparison.

## Contents and authority

- [evidence/](evidence/): 90 copied light reports, native metrics and structural
  audits, with [checksums](evidence_checksums.json). No large VCFs are committed.
- [heavy_artifacts.json](heavy_artifacts.json): 75 retained heavy artifact entries
  with absolute locations and sizes. Ignored aliases under
  `.artifacts/benchmark_archive/frozen_native_gate_20260916/` organize the three
  evaluation roots without moving them or breaking provenance paths.
- [reproduce.sh](reproduce.sh): executable commands for all three full-input
  evaluations and nine structural audits. Run from either repo location with
  `bash docs/manuscript/frozen_native_gate_20260916/reproduce.sh /absolute/new-results`.
  Source manifests retain the recorded absolute input paths; moving inputs needs
  explicit manifest remapping and checksum verification. A new result directory
  is mandatory. [verify.py](verify.py) checks evidence and policy code hashes
  before reproduction; run it alone for a lightweight archive check.
- Export utility: `examples/seqc2/scripts/archive_frozen_benchmarks.py --archive
  /new/archive --heavy-links /ignored/aliases` (requires a provenance.json in the
  new archive). Existing archive exports are never overwritten.
- [metrics.csv](metrics.csv): 108 rows, comprising three datasets, two targets,
  six methods and three variant types. Counts and precision/recall/F1 come
  directly from completed full-input som.py evaluations, without recomputation.
- [provenance.json](provenance.json): source report locations and SHA-256,
  exact input manifests and source hashes, code hashes, output VCF hashes and
  the nine completed full-input structural audit references.
- Original large VCFs, SQLite evidence stores and benchmark scratch remain at
  their original locations. No source file or cache was moved or deleted.
  This is an indexed metadata/results archive, not a backup of those large
  artifacts. Retention and storage backup must preserve the referenced files.

The authoritative evaluations are under
`examples/seqc2/comparison/current_native_audit_20260916_v2/`:

| Dataset | Final full-input evaluation directory |
|---|---|
| SEQC2 WES-LL | `seqc2_wes_frozen_full_domain_v1` |
| SEQC2 WGS-IL | `seqc2_wgs_frozen_full_domain_v1` |
| HG008 WGS | `hg008_frozen_full_domain_v1` |

Within each directory, `validation.json` contains all executed commands and
metrics. `refined.vcf.gz`, `first/refined.rescue.vcf.gz` and
`realignment/refined.rescue.vcf.gz` are biological-label outputs. The separate
`*.query.vcf.gz` PASS copies are for benchmarking only, not training.

## Headline results

Aggregate record F1 (per-type results remain in the CSV):

| Dataset/target | DNA DeepSomatic | Optimized consensus | + First rescue | + Realignment rescue |
|---|---:|---:|---:|---:|
| SEQC2 WES / UKB | 0.619019 | 0.621199 | 0.624301 | 0.625294 |
| SEQC2 WES / MedExome | 0.798016 | 0.799433 | 0.802539 | 0.803383 |
| SEQC2 WGS / UKB | 0.966347 | 0.968346 | 0.968346 | 0.968346 |
| SEQC2 WGS / MedExome | 0.913468 | 0.916070 | 0.916070 | 0.916070 |
| HG008 / UKB | 0.810300 | 0.841043 | 0.837521 | 0.838926 |
| HG008 / MedExome | 0.830189 | 0.828338 | 0.826087 | 0.828338 |

Method identifiers in the CSV:

- `deepsomatic`: DNA DeepSomatic native passing calls.
- `refined_consensus`: experimental native-evidence SNP and corroborated-indel consensus.
- `refined_first_rescue`, `refined_realignment_rescue`: that DNA baseline plus
  the experimental nomination/eligible-RNA/biological gate, separately per round.
- `workflow_first_rescue`, `workflow_realignment_rescue`: existing annotated
  workflow results; these are historical comparators, not the optimized adapter.

## Experimental contract and limitations

All supplied DNA caller records entered candidate selection. Truth HC `-R` and
UKB/MedExome `-T` were used at scoring only. This does not imply whole-genome
coverage of the source callers. som.py used the recorded GRCh38 reference and
`-N`, not `-P`. Native aggregate record metrics must not be replaced with a
sum of SNP and indel rows, which can differ due to representation.

SEQC2 WES/WGS were development data. HG008 is a subsequent frozen evaluation,
but older development already used HG008, so it is not a pristine holdout.
UKB and MedExome overlap; WES/WGS SEQC2 share biological material. These are
not six independent biological validations. No confidence interval or
significance claim is supplied by this archive.

Consensus aggregate F1 improves over DeepSomatic in five of six target-specific
comparisons, not every variant class. In particular, WES MedExome indel F1 and
HG008 MedExome aggregate F1 remain lower. Rescue helps SEQC2 WES, has no scored
gain on SEQC2 WGS, and adds no net TP on HG008. No universal superiority or
cross-tool SOTA claim is justified. The six methods here do not constitute a
new comprehensive comparison against every standalone caller or ClairS.

The WES first and realignment source artifacts come from distinct historical
workflow versions. They must not be described as paired stages of a newly
executed workflow. Current experiments reuse existing annotated rescue
candidates, not newly generated annotation/realignment output.

All structural audits passed and all evaluation integrity checks passed.
Structural consistency is not biological truth. The HG008 baseline conflict
review confirms 48 exact gnomAD AF matches; UKB partitions include 34 TP,
6 FP and 8 unmatched alleles. A blanket AF veto is not adopted. This evidence
does not resolve the biology of each conflict or approve all training labels.

## Historical versus final evidence

Earlier region-scoped assays, native-gate replays, rejected rule ablations and
failed source-discovery attempts remain preserved as development history.
They are not substitutes for the final full-input directories above and must
not be silently combined into a same-policy result table. All eight optimized
SEQC2 rescue comparisons reproduce their earlier region-scoped metrics, but
whole-output VCF identity is not implied.

Supporting records:

- [Full-input validation and conflict review](../../review/FROZEN_VALIDATION_REMAINING_WORK_20260916.md)
- [HG008 detailed performance](../../review/HG008_FROZEN_POLICY_RESULTS_20260916.md)
- [Rule selection and rejected hypotheses](../../review/SEQC2_INDEL_RESCUE_FOLLOWUP_20260916.md)
- [Evidence-preserving rescue adapter](../../review/SEQC2_RESCUE_OUTPUT_ADAPTER_20260916.md)

## Cohort rollout boundary

Freeze this version for a versioned cohort pilot rather than retuning against
HG008. Further optimization is a separate SEQC2-only development experiment,
not a prerequisite for generating candidate cohort labels.

Global workflow-default integration is not required before a cohort rerun.
A dedicated, tested opt-in cohort driver is required: the older cohort config
does not invoke this optimized consensus and standalone gate. Preserve both DNA
consensus and rescued outputs, existing source VCFs/indexes, and separate output,
state and checksum namespaces. Reuse matched RNA realignment caller evidence;
do not rerun mapping/calling. Retain or regenerate downstream annotation as
explicitly required by the final artifact contract.

Full biological training-label approval need not precede generation into an
isolated experimental output directory. It must precede using those candidates
as approved training labels or replacing previous labels. Before launching all
66 samples, require input/evidence preflight, a representative pilot, bounded
resource/concurrency settings and explicit per-sample QC/completion checks.
No cohort jobs were launched while preparing this archive.
