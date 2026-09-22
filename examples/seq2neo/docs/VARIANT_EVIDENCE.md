# Approved-release variant evidence, schema seq2neo-info-v1

`bin/export_variant_evidence.py` recomputes read evidence from the three BAMs
explicitly registered in the hash-bound approved 63-sample release. It never
uses historical measurements or calls a classification, alignment, variant
calling, quantification, or training workflow. The approved FILTER is copied
unchanged. Immutable stale eligibility flags are deliberately not selection
criteria: RELEASE_APPROVAL.json authorizes these exact source bytes.

Output must be outside the approved release and historical stats directories.
Development files live under `development/<full sample_id>/part-*.parquet`;
reserved files under `reserved/<full sample_id>/part-*.parquet`. Read only those
parquet files, never JSON receipts or temporary files. A bounded pilot is an
explicit partial release and is not the full dataset.

## Data dictionary and feature boundary

The five-column primary key is `sample_id, CHROM, POS, REF, ALT`. POS is the
original 1-based VCF position. BAM queries use 0-based half-open coordinates.
All original keys remain unchanged; no normalization, decomposition, key
collapse, sample-suffix lookup, or additional label filtering occurs. Duplicate
full keys abort. Distinct ALT rows at one locus are counted independently.

`FILTER` is approved label provenance, not a model feature. `pool` is
`train_pool` or `reserved`. `split` is `reserved` for every reserved record;
otherwise chr1=test, chr21/chr22=val, everything else=train, following the
2026-09-22 decision in TRAINING_DATA_GUIDE.md. Sample identity, chromosome split,
labels, and missingness explanations are excluded from the numerical feature
allowlist in `schema.json`. No caller, rescue, tier, rationale, or historical
annotation columns are copied. Labels remain workflow-derived weak labels.

For each prefix `bam_DN_`, `bam_DT_`, `bam_RT_`:

| Suffix | Meaning |
|---|---|
| depth | Usable aligned read observations spanning this allele; an indel also requires a sequenced right flank |
| ref_count, alt_count, other_count | Exact REF, exact ALT, other observed allele; sum equals depth |
| vaf_denominator | depth, including other alleles (not REF+ALT only) |
| vaf | alt_count / vaf_denominator, null if denominator zero or unavailable |
| ref_forward, ref_reverse, alt_forward, alt_reverse | Support by alignment strand |
| ref_f1r2, ref_f2r1, alt_f1r2, alt_f2r1 | R1 forward/R2 reverse vs R1 reverse/R2 forward; paired read flags required |
| ref_orientation_unknown, alt_orientation_unknown | Unpaired or ambiguous read-number flags |
| mean_bq | Mean per-observation minimum BQ across allele bases, inserted bases and required indel right flank |
| mean_mq | Mean effective MQ across usable observations, with STAR convention below |
| mq_count | Number of usable observations contributing effective MQ |
| status, missing_reason | Measurement state and explanation; never silently fill null measurements with zero |

`expression_status=unavailable_no_versioned_quantification` always. RNA locus
depth is not gene expression. No TPM, FPKM, gene assignment or transcript
abundance is inferred.

## Frozen read semantics

Exclude unmapped, secondary, QC-failed, duplicate-flagged and supplementary
alignments (SAM mask 3844). Require effective MQ >=20 and minimum event BQ >=20.
Non-ACGT observed bases and absent qualities are unusable. No BAQ, realignment,
pair properness, insert-size, or strand restriction. Duplicate exclusion trusts
the BAM flag; these mapped BAMs may precede duplicate marking and unmarked PCR
duplicates can remain. Overlapping mates count as two read observations,
explicitly not independent fragments. Orientation counts are read counts.

STAR MAPQ=255 is accepted only when the BAM header declares STAR and read NH=1,
using effective MQ=60. Other MQ=255 observations are excluded. The original
header and program chain are preserved in each sample's `inputs.json`; no
unrecorded SplitNCigarReads/BQSR lineage is assumed. CIGAR N is a splice gap,
never deletion evidence. Soft/hard clips do not advance the reference position.

Supported alleles are uppercase A/C/G/T SNVs and simple VCF-anchored insertions
or deletions. Indel counting reconstructs the observed sequence over the whole
REF interval and requires a sequenced right flank. It counts complete inserted
sequence/deleted interval, never just the shared anchor. This is literal
alignment evidence without haplotype realignment; equivalent shifted indels can
be missed. Complex replacements, MNVs, symbolic alleles, comma-separated ALTs,
ambiguous bases and identical REF/ALT are explicitly unsupported and remain
rows with null measurements. Separate biallelic ALT rows are supported.

Each supported REF must exactly match the supplied hash-recorded FASTA. Contig
names are exact (no chr aliases); reference and BAM contig lengths must agree.
Headers, read-group sample IDs, reference identities and index identities are
recorded. Expected BAM RG SM is full sample_id plus DN/DT/RT; disagreement makes
that modality unavailable. A header reference path/length is provenance rather
than independent certification that every alignment used identical bases.

Depth cap is 10,000 usable observations per allele. More than that produces
`depth_cap_exceeded` and null metrics, not truncated estimates. Accumulators
stop growing after cap detection, although alignment traversal continues.
`zero_usable_depth` has count fields zero and VAF/quality means null. Covered
loci with zero ALT have `status=ok`, ALT=0 and VAF=0. Absent BAM, absent index,
unreadable BAM, sample mismatch, query failure, missing contig, reference
mismatch, invalid coordinate and unsupported allele have distinct statuses;
all measurements are null. A query error invalidates every row in that query
window rather than publishing partially read counts.

## Bounded processing, provenance and resume

The approved parquet is projected to six columns and processed one sample per
worker, then in batches (default 10,000 records). Queries group loci in 100 kb
windows, fetching only the first through last target in a window. Original keys
and FILTER are compared exactly after writing and on every resume; no key-only
join can silently multiply rows. Input sample/class counts are reconciled to
the approval at startup and all approval hashes are verified again afterward.

Every BAM, index, FASTA, FAI and approved input is SHA-256 hashed. The manifest
records extractor revision, relevant current source/document/lockfile hashes,
dirty status and library versions. File identities include realpath, size,
mtime, ctime, inode and device. Input metadata is checked again after extraction
to detect concurrent mutation. Resume rehashes the inputs and each output and
requires matching extractor/settings/allele identity; existence alone is never
enough. Change extraction code/settings in a new version directory. The output
lock prevents concurrent writers. Files and completion receipts use atomic
rename; per-part receipts and sample QC preserve completed progress.

`manifest.json`: global provenance and approval acknowledgement.
`schema.json`: exact Arrow types and numerical model-feature allowlist.
`*/<sample_id>/inputs.json`: full BAM/index hashes and headers.
`*/<sample_id>/part-*.json`: content-bound completion receipts.
`*/<sample_id>/qc.json`: classes, split counts, modality statuses, covered rows,
ALT-positive rows, runtime and process resource use.
`qc.json` and `run-*.json`: cohort reconciliation, completeness and exact command.
RSS is the process high-water mark; worker CPU counters are cumulative if that
worker processes several samples. Full-release completion is true only when all
63 samples and the exact 14,342,597 approved rows/classes reconcile.

## Reproduction

Run from the repository root using the existing `.venv` (do not bare uv sync):

```bash
.venv/bin/python -m pytest tests/test_export_variant_evidence.py -q

.venv/bin/python -u bin/export_variant_evidence.py \
  --release-root examples/seq2neo/output_three_class_v2_20260919 \
  --reference /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
  --output examples/seq2neo/info_v1_20260922_pilot_r2 \
  --sample PRJNA298376_4007 --pilot-rows 2000 --batch-rows 1000
```

Append `--resume` to exactly the same command to verify/reuse completed parts.
For the full cohort use a new output directory
`examples/seq2neo/info_v1_20260922`, omit `--sample` and `--pilot-rows`, and choose
`--workers 2`. Identical frozen measurement settings apply to reserved data.
`--max-runtime-seconds` optionally imposes a resumable scheduling budget checked
between batches and before queued samples begin hashing; it does not interrupt
an active hash or batch. Absent that option there is no time limit. Rehashing inputs is
included in runtime. Concurrency affects scheduling, not counting semantics.

An input-identity mismatch on resume is rejected before the existing sample
provenance can be overwritten. The standalone verifier additionally detects
duplicates across part files and offers independent htslib-based SNV and indel
checks for development samples:

```bash
.venv/bin/python bin/verify_variant_evidence.py \
  --output examples/seq2neo/info_v1_20260922_pilot_r2 \
  --allow-partial --check-bam-snvs 25 --check-bam-indels 12
```

Omit `--allow-partial` when verifying the full release; an incomplete cohort then
fails. Verification of held-out data covers membership and numerical integrity;
the optional BAM spot checks are restricted to development samples.

## Why historical BAM extraction is not reused

The 66 historical parquets and launcher were inspected only as a schema/code
reference. The old launcher selects a superseded manifest. The existing Rust
pileup uses reference offset directly as query offset, truncates REF/ALT to one
byte and names strand-only counts F1R2/F2R1. Those semantics are unsuitable for
spliced RNA or indels; missing index/query paths can also silently default.
The new implementation retains the useful grouping and complete-key validation
approach but uses pysam indexed readers and explicit CIGAR projection.
