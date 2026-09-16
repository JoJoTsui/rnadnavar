# Experimental refined rescue output adapter

## Scope and release boundary

`bin/apply_refined_rescue.py` is a standalone post-annotation adapter using
existing refined DNA consensus, annotated rescue and six caller VCFs. It does
not invoke Nextflow, mapping, calling, annotation downloads or HG008 evaluation.
No workflow config enables it. Outputs are explicitly marked
`experimental_not_training_ready` and must not replace production labels.

The validation currently uses a **benchmark-domain-scoped DNA consensus**.
Even though the adapter retains the union of its inputs, that does not make
the resulting VCF a full-genome refined DNA label set. Whole-domain baseline
generation, annotation consistency checks and label QC remain release gates.

## Decisions and evidence

The reconstructed gate is documented in
[the refinement report](SEQC2_INDEL_RESCUE_FOLLOWUP_20260916.md#rescue-eligibility-reconstruction-and-benchmark-validation).
The adapter adds the following output safeguards:

1. Keep the refined DNA Somatic baseline. Count baseline conflicts with the
   available population/editing annotations for subsequent review.
2. Preserve existing rescue negative labels; protected DNA Artifact, Germline,
   Reference or RNAedit cannot be promoted by the new gate.
3. For an existing Somatic rescue SNP not in the DNA baseline, require the
   positive DNA nomination, two distinct eligible RNA callers and biological
   gates. Rejected/inconclusive DNA verification withholds an otherwise
   qualifying addition. Gate-rejected positives become NoConsensus, not PASS.
4. Preserve baseline indels, but do not admit indels through RNA rescue.
5. Update FILTER, UNIFIED_FILTER, the current classification rationale,
   UNIFIED_FILTER_DNA, PASSES_CONSENSUS_DNA and rescued/promotion flags together.

The chosen source record's annotations and per-caller fields are retained.
`GATE_SOURCE_RECORD` stores a percent-encoded copy of that record before
relabeling. Where both inputs contain the allele, `GATE_DNA_RECORD` preserves
the refined DNA record separately. Decode once to recover the original
eight-column record. These are explicitly historical snapshots, not competing
active classifications. `GATE_DNA_NOMINATORS` and `GATE_RNA_ELIGIBLE` record
the reconstructed caller identities; `GATE_ALIGNMENT_ROUND` identifies the
supplied RNA round without rewriting historical evidence origins.

Source evidence may contain old round tags or annotation conflicts. Retaining
it does not certify or repair it; this is another reason the output is not yet
a production training-label release.

## Safety and scalability

- Caller identity is supplied explicitly as `caller=path`, not guessed from
  filenames. Each modality requires three distinct caller files.
- VCF sample roles use the existing header/name/caller conventions, not SEQC2
  sample-name literals. More than two caller samples is rejected.
- SQLite stores the full union and compressed source records on disk, with a
  32 MiB page-cache limit. Caller files are scanned sequentially; duplicate
  observations never create extra caller votes.
- Original inputs and indexes are read-only. Before/after SHA-256 checks cover
  all eight inputs. Existing output directories are refused.
- The SQLite database is retained as an experimental artifact. It is not a
  Nextflow cache and may be substantially larger than the compressed VCFs.

## Validation entry point

From the repository root, use a new output directory for every attempt:

```bash
.venv/bin/python examples/seqc2/scripts/validate_refined_rescue_output.py \
  --root examples/seqc2/comparison/current_native_audit_20260916_v2 \
  --dataset wes --region ukb \
  --outdir examples/seqc2/comparison/current_native_audit_20260916_v2/rescue_output_wes_ukb_v1
```

The validation driver processes first rescue and realignment rescue separately,
using their matching RNA caller bundles. It writes full biological-label VCFs,
then separate PASS-only benchmark queries. som.py inherits the frozen truth,
reference, HC `-R` and selected target `-T`; it does not use `-P`.

Outputs for each round:

- `refined.rescue.vcf.gz` and index: experimental evidence-preserving VCF.
- `report.json`: input hashes, label/decision counts and integrity status.
- `union.sqlite`: on-disk union and source records.
- `benchmark.query.vcf.gz`, metrics and logs: benchmark-only artifacts.
- Parent `validation.json`: commands, comparison metrics and round completion.

The WES first-round source is the policy-default workflow; its realignment
source is the archived `realign.latest` workflow recorded by the verified
bundle. Their provenance is distinct; this is not a claim that the two rounds
were generated together in a new pipeline execution.

## Validation results and remaining checks

Completed evidence-preserving output benchmarks (aggregate records):

| Dataset/target | Rescue round | TP | FP | FN | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| WES/UKB | First | 1061 | 38 | 1239 | 0.965423 | 0.461304 | 0.624301 |
| WES/UKB | Realignment | 1063 | 37 | 1237 | 0.966364 | 0.462174 | 0.625294 |
| WES/MedExome | First | 569 | 20 | 260 | 0.966044 | 0.686369 | 0.802539 |
| WES/MedExome | Realignment | 570 | 20 | 259 | 0.966102 | 0.687575 | 0.803383 |
| WGS/UKB | First | 2172 | 14 | 128 | 0.993596 | 0.944348 | 0.968346 |
| WGS/UKB | Realignment | 2172 | 14 | 128 | 0.993596 | 0.944348 | 0.968346 |
| WGS/MedExome | First | 704 | 4 | 125 | 0.994350 | 0.849216 | 0.916070 |
| WGS/MedExome | Realignment | 704 | 4 | 125 | 0.994350 | 0.849216 | 0.916070 |

All eight rounds are complete and report unchanged input hashes. All four
realignment assays reproduce the frozen replay's SNP, indel and aggregate
metrics. Completion does not remove the benchmark-domain baseline limitation.

All paths below are relative to
`examples/seqc2/comparison/current_native_audit_20260916_v2/`:

- `rescue_output_wes_ukb_v1/validation.json`: both rounds complete.
- `rescue_output_wes_medexome_v1/validation.json`: both rounds complete.
- `rescue_output_wgs_ukb_v2/validation.json` and
  `rescue_output_wgs_medexome_v2/validation.json`: both rounds complete.
- The two WGS `v1` attempts failed source discovery before adapter execution.
  They are preserved with `failed_source_discovery` reports, not benchmark
  results. The driver now selects WGS first-round `*.filtered.vcf.gz` separately
  from realignment `*.rescue.filtered.stripped.vep.vcf.gz`; a regression test
  covers the four dataset/round combinations.

### Baseline annotation conflicts: investigated, not silently filtered

The WES/UKB output flags 14 retained DNA-baseline alleles with gnomAD AF above
0.001. Read-only audits are retained in
`baseline_conflict_audit_v1/audit.json` and
`baseline_conflict_audit_v2/audit.json`; v2 adds paired caller evidence.
All 14 AF values match the exact REF/ALT allele in the configured local gnomAD
exomes v4.1 database. Saved som.py query partitions classify **10 as FP and
4 as TP**. A blanket AF veto on these baseline alleles would therefore remove
four truth positives as well as ten false positives; it is not adopted.

The paired evidence does not provide a clean normal-alt-count discriminator:
13 sites have Strelka normal evidence, all with tier-one alternate count zero;
normal DP spans 5–103. The remaining FP has no paired-caller record. DeepSomatic
exports only the tumor sample here, so its missing normal FORMAT evidence
must not be interpreted as zero normal support. Only two sites have Mutect2
records (one FP and one TP). Missing records are not negative evidence.
Database record filters are preserved in the audit as well: exact AF matching
does not imply that every gnomAD record itself is PASS.

The four TP sites are `chr14:106268545:C:A`, `chr14:106324560:C:A`,
`chr17:81906207:G:A`, and `chr22:29141977:G:A`. These are benchmark labels,
not a determination of why a population-database match occurs. No new
threshold, baseline veto, workflow default or training-label release follows
from this small, selected subset. The existing population veto still applies
to rescue additions, separately from baseline retention.

### Resource and implementation checkpoint

The first WES/UKB output assay stages 999,080 union records. Observed process
RSS during staging/writing was roughly 70–80 MB. The assay is retained under
`rescue_output_wes_ukb_v1/`; check its `validation.json` and per-round reports
for completed metrics rather than treating this progress note as a pass.

The first-round process was launched before additional contig/union guards
and conflict counters were added; its report records that earlier script hash.
The later realignment process uses the updated adapter. These changes do not
alter the decision rules on valid inputs; fixture tests cover the safeguards.
Both sources were independently checked to have compatible INFO definitions
and shared-contig lengths.

The relevant Python regression run passed 259 tests (62 existing/dependency
warnings). No mapping, variant calling, Nextflow run, HG008 evaluation, source
index rewrite or workflow-default promotion was performed for these checks.

Subsequent structural audit: the complete WES/UKB realignment output passed
all implemented structural checks across 441,954 records with zero reported
issues and unchanged source checksum. Evidence:
`rescue_output_wes_ukb_v1/realignment/structural_audit_v1.json` under the root
above. All remaining seven output scans subsequently completed with zero
reported structural issues and unchanged source checksums. This is not
biological training-label approval.
HG008 evaluation was subsequently authorized and is recorded in
[the frozen full-domain validation runbook](HG008_FROZEN_FULL_DOMAIN_VALIDATION_20260916.md).
