# v2 three-sample model-interface review

Superseded execution scope: the user confirmed training uses another host
and another user's model code, and requested the full cohort instead.
See [COHORT_HANDOFF.md](COHORT_HANDOFF.md). Findings below apply only to
the inspected local copy, not the destination implementation.

Scope: local `/t9k/mnt/hdd/work/Vax/EvoSomatic` copy, without Git metadata.
The model owner's current execution repository/version has not been confirmed.
This is a CPU-only parser probe, not a cache build or training run.

## Somatic quarantine review

Of the 3,464 quarantined Somatic records in the 63-sample proposal:

- 3,462 have `GNOMAD_AF > 0.001`: 3,007 exceed 0.01 and 455 are in
  (0.001, 0.01]. All use `branch:retained_dna_baseline`.
- Two use `branch:dna_nominated_rna_supported` and conflict with native
  negative nominations, as documented in the paired-read pilot.
- There are 3,330 SNPs and 134 indels in total.

This is not an export bug: the frozen rescue gate explicitly applies biological
vetoes to additions, not retained DNA baseline calls. These exclusions are a
conservative **training-subset proposal**, not retrospective proof that every
record is a false positive, nor an update to the benchmark callset. Approval
must name this choice. Population annotations were read from the v2 export;
this review did not independently reannotate gnomAD or establish genotype truth.

## Concrete interface findings

1. The inspected model expects JSON `schema_version=1` with full sample ID,
   `tumor_dna`, `normal_dna`, `tumor_rna`, and a labeled `vcf` per sample.
   It does not consume the new TSV/Parquet directly. Existing local three-sample
   config still points to old labels in Zhan's Dataset directory.
2. Its FILTER mapping is Reference=0, Germline=1, Somatic=2. It ignores approval
   fields. Do not launch it on review files merely because parsing succeeds.
3. Its split helper uses chr1=test, chr21/22=validation, **every other header
   contig=train**, including chrX/Y/M. The existing cohort split documentation
   instead names chr2–20 for training. Resolve this explicitly before release.
4. Cache build identity includes the manifest's path/size/mtime but does not
   hash the VCF/BAM contents referenced by the manifest. Use an immutable
   release-bound cache namespace; do not assume old packed caches match v2.

## Verified bridge

Root:
`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/model_review_bridge_20260921`.

`samples.review.json` references three new indexed `SAMPLE.review.vcf.gz`
files. These minimal VCFs contain exactly the selected Parquet alleles/classes;
annotation-rich evidence remains in the original handoff Parquet. All 976,677
rows passed ordered allele/class hash roundtrip checks. No candidate labels or
original outputs changed. `TRAINING_ELIGIBLE=NO` and manifest
`training_approved=false` remain explicit.

The actual local model loader and candidate parser imported successfully and
parsed all three samples. Report: `local_loader_probe.json`; exact inspected
code hashes are included. It returned 976,675 candidates, omitting two Germline
records outside its allele-length scope. Coordinate conversion uses the model's
zero-based candidate anchor, not the one-based VCF POS.

| Sample | Bridge rows | Parsed | Scope omissions | Non-autosomal training rows |
| --- | ---: | ---: | ---: | ---: |
| PRJNA298376_4007 | 72,958 | 72,957 | 1 | 1,883 |
| PRJNA298376_4060 | 298,540 | 298,540 | 0 | 9,384 |
| PRJNA298376_4072 | 605,179 | 605,178 | 1 | 28,244 |

Thus the local default assigns 39,511 sex-chromosome/mitochondrial records to
training. This is an observed scope mismatch, not a claim those records are
biologically wrong. The inspected model supports only indel length differences
<50 bp and normalized alleles <=50 bp; record omitted candidates in accounting.

## Reproduction

Focused verification: `test_v2_model_review_bridge.py`,
`test_v2_training_handoff.py`, and `test_three_class_export_verification.py`
passed together (10 tests). These checks cover the review bridge and export
contracts, not model training or biological label approval.

Fresh destinations required; this does not initialize Evo2/GPU or build features:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 prlimit --as=17179869184 -- \
  .venv/bin/python examples/seq2neo/scripts/build_v2_model_review_bridge.py \
  --handoff /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_three_class_v2_20260919/handoff_review_20260921 \
  --outdir /path/to/fresh/model_review_bridge

PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 prlimit --as=17179869184 -- \
  .venv/bin/python examples/seq2neo/scripts/probe_v2_model_review_loader.py \
  --model-root /t9k/mnt/hdd/work/Vax/EvoSomatic \
  --bridge /path/to/fresh/model_review_bridge \
  --out /path/to/fresh/model_review_bridge/local_loader_probe.json
```

Next: confirm the actual model execution copy, approve the scoped weak-label
release and resolve the chromosome split; then check the final configuration
and launch the **three samples together**, not a separate 4007-only stage.
Reference-specific 299-read research remains outside these release requirements.
