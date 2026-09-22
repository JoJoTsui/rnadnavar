# Training Data Guide — Search, Split, and Preparation

The single source of truth for **where the training data is**, **how it is split into
train/val/test**, and **how to prepare it for model training** on the
`separated_three_class_v2` cohort.

This guide targets the **approved 63-sample release of 2026-09-21**. If you are
looking for "the manifest and the parquet to train on," start at
[§1 Artifacts](#1-artifacts--authoritative-paths). If you are about to write a
loader, read [§3 Critical gotchas](#3-critical-gotchas-read-before-writing-a-loader)
first — the obvious columns are booby-trapped.

Supersedes the pre-release framing in [`data_split_strategy.md`](data_split_strategy.md)
for training-data purposes (that doc remains the design rationale and the older
reconsensus-cohort record; where the two disagree on other-contig handling, this
guide governs — the recorded decision routes them to **train**, §5.3).

**Scope.** This guide covers **variant-level features only** — the columns of
`cohort63.review.parquet` (labels, per-caller evidence, DP/VAF aggregates) plus
the train/val/test fold assignment. **Read-level tensor construction** (pileups,
per-read / image-like tensors from the DT/DN/RT BAMs) is deliberately **out of
scope** and left to each specific model implementation. The BAM paths are
provided here (§2.1, §2.3) for a model's loader to consume.

---

## 1. Artifacts — authoritative paths

All paths are relative to the release root:

```
examples/seq2neo/output_three_class_v2_20260919/
```

which exists under **both** roots (identical content; the SHA-256 pins, not the
path, are authoritative):

- repo working tree: `/t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/`
- shared/rsynced root (the `shared_root` recorded in the approval sidecar):
  `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/`

| Role | File | Notes |
|------|------|-------|
| **Approval sidecar (authoritative)** | `cohort63_handoff_tools_20260921/RELEASE_APPROVAL.json` | Hash-bound release approval. **This is what "approved" means.** |
| **Variant manifest (TSV)** | `handoff_review_20260921/cohort63.review.tsv` | 63 rows, one per sample: identity + alignment paths + provenance. |
| **Variant parquet** | `handoff_review_20260921/cohort63.review.parquet` | 14,342,597 rows, one per variant: `CHROM/POS/REF/ALT/FILTER` + rich annotation. **This is the training tensor source.** |
| Train-pool sample list | `cohort63_model_review_20260921/train_pool.review.json` | 58 samples → training. Also carries DN/DT/RT BAM paths. |
| Reserved sample list | `cohort63_model_review_20260921/reserved.review.json` | 5 samples → held-out evaluation **only**. |
| Alignment inventory | `cohort63_model_review_20260921/report.json` | 189 DN/DT/RT paths + indexes verified present. |
| Handoff narrative | `cohort63_handoff_tools_20260921/COHORT_HANDOFF.md` | Human-readable summary of the release. |
| Selection provenance + hashes | `handoff_review_20260921/audit.json` | Per-stage counts and hashes. |

SHA-256 pins (from `RELEASE_APPROVAL.json`):

| Artifact | SHA-256 |
|----------|---------|
| `cohort63.review.parquet` | `8e61c3c73b51cbca095b0da26e34f76b5adfe80c594723acaaedf4d740662a6d` |
| `cohort63.review.tsv` | `05602cf063e96e15d9c4c9485a854446c79d0d8d39c8fc68e9e7bf3664ac4cdc` |
| `train_pool.review.json` | `499f9fbc9e03b7bd4aa5b26fe2a475acb0bf54bcdcda29ab1325afab7071b11c` |
| `reserved.review.json` | `219c9fd8c1e0da700636b7963a84a256914cd7f27243477ac78f509725838843` |
| `report.json` (verification) | `eb6f11ef47e8c7c2e9c32853de2843894c98fa68d640f6dbb8d4b96a3cb554b8` |

**Verify before training.** A directory existing is not a successful release; the
hashes must match:

```bash
cd examples/seq2neo/output_three_class_v2_20260919/handoff_review_20260921
sha256sum -c <<'EOF'
8e61c3c73b51cbca095b0da26e34f76b5adfe80c594723acaaedf4d740662a6d  cohort63.review.parquet
05602cf063e96e15d9c4c9485a854446c79d0d8d39c8fc68e9e7bf3664ac4cdc  cohort63.review.tsv
EOF
```

### Do NOT use these for training

| File | Why not |
|------|---------|
| `output_three_class_v2_20260919/exports/` | `candidate_not_training_approved` — export stage, pre-approval. |
| `data/processed/refined_three_class_20260918/`, `refined_native_v2_20260917/` | `candidate_not_training_approved`, `label_qc_verdict: NOT_APPROVED`. |
| `data/processed/sample_manifest*.tsv`, `sample_manifest.parquet` (2026-06) | Original cohort, superseded. |
| `data/processed/training_manifest_65.tsv` | Explicitly superseded by the rerun manifest; history only. |
| `data/processed/selected_variants.parquet` + `sample_split.tsv` (2026-09-03) | The older reconsensus split build. Useful as a cross-check and for its materialized `split` column (see §5.4), but **not** the approved release. |

---

## 2. Data model

### 2.1 Manifest — `cohort63.review.tsv` (sample level, 63 rows)

One row per sample. Grouped by purpose:

- **Identity**: `sample_id`, `project_id`, `patient_id`, `set_number`, `disease`,
  `disease_normalized`, `status` (`standard` 31 / `extra` 32), `status_reason`.
- **Alignment (BAM) paths**: `bam_dn` (normal DNA), `bam_dt` (tumor DNA),
  `bam_rt` (tumor RNA).
- **Per-caller VCF paths**: `caller_{dna,rna}_{mutect2,deepsomatic,strelka}`.
- **Label provenance**: `original_rescue_vcf_path`, `consensus_vcf_path`,
  `rescue_vcf_path`, `training_label_vcf`, `new_rescue_sha256`.
- **Stage/provenance**: `variant_parquet_path`, `source_variant_parquet`,
  `candidate_status`, `candidate_policy` (`separated_three_class_v2`),
  `label_qc_verdict`, `handoff_status`, `is_complete`, `vcf_prefix`, `dir_name`,
  `base_output_dir`.

> **The manifest has no split/pool column.** It does not say which samples are
> train vs reserved, and it carries no train/val/test assignment. Pool membership
> lives in `train_pool.review.json` / `reserved.review.json` (join on `sample_id`);
> the train/val/test split is a property of each **variant's** `CHROM` (§5).

### 2.2 Parquet — `cohort63.review.parquet` (variant level, 14,342,597 rows)

The training tensor source. Key columns:

- **Identity**: `sample_id`, `CHROM`, `POS`, `REF`, `ALT`, `variant_type`
  (`SNP` 14,239,387 / `INDEL` 103,210).
- **Label**: `FILTER` — **this is the training label** (see §4).
- **Per-caller evidence** — encoded as `caller:value` pairs joined by `|` (e.g.
  `AD_BY_CALLER` = `deepsomatic:43,2|strelka:44,2`; `.` = missing call; no JSON):
  `GT_BY_CALLER`, `DP_BY_CALLER`, `AD_BY_CALLER`, `VAF_BY_CALLER`, and the
  `NORMAL_*` equivalents; plus the pre-counted `N_DNA_CALLERS_SUPPORT`,
  `N_RNA_CALLERS_SUPPORT`, `N_DNA_CALLERS_SOMATIC`.
- **Rescue**: `RESCUED`, `RESCUE_PROMOTED`, `PASSES_CONSENSUS_DNA`.
- **Population / curation**: `GNOMAD_AF`, `REDI_ACCESSION`, `REDI_CANONICAL`.
- **Rationale / review**: `CLASSIFICATION_RATIONALE`, `review_reason`,
  `THREE_CLASS_*`, `NEGATIVE_EVIDENCE_*`, `GATE_*`, `DNA_VERIFICATION`,
  `UNIFIED_FILTER_DNA`.
- **Quantitative features** (ready-made tensor inputs): `DP_DNA_MEAN`,
  `DP_RNA_MEAN`, `VAF_DNA_MEAN`, `VAF_RNA_MEAN`.
- **Stale flags — see §3**: `candidate_status`, `label_confidence`,
  `training_eligible`, `TRAINING_ELIGIBLE`.

The split is **not** a column here. It is derived from `CHROM` (§5).

### 2.3 Sample pool manifests — `*.review.json`

```jsonc
{
  "schema_version": 1,
  "training_approved": false,          // stale pre-approval snapshot — see §3
  "samples": [
    { "sample_id": "...",
      "tumor_dna": "/t9k/.../PRJNA298330_3812DT.sorted.bam",
      "normal_dna": "/t9k/.../PRJNA298330_3812DN.sorted.bam",
      /* ... tumor_rna, and index paths ... */ }
  ]
}
```

`train_pool.review.json` = 58 training samples, `reserved.review.json` = 5
held-out samples. Join to the manifest/parquet on `sample_id`.

---

## 3. Critical gotchas (read before writing a loader)

### 3.1 The in-file approval flags are stale. The sidecar is the approval.

Every row of the approved parquet and manifest still says **"not approved"**:

| Column / field | Value in the approved release | Meaning |
|----------------|------------------------------|---------|
| `training_eligible` (parquet) | `False` for all 14,342,597 rows | **Stale.** |
| `TRAINING_ELIGIBLE` (parquet) | `"NO"` for all rows | **Stale.** |
| `candidate_status` (both) | `candidate_not_training_approved` | **Stale.** |
| `label_qc_verdict` (manifest) | `NOT_APPROVED` for all 63 | **Stale.** |
| `label_confidence` (parquet) | `unvalidated` for all rows | **Stale.** |
| `handoff_status` (manifest) | `review_ready_not_training_approved` | **Stale.** |
| `training_approved` (pool JSONs) | `false` | **Stale.** |

These are deliberate **immutable pre-approval snapshots**. Approval is recorded
*only* in `RELEASE_APPROVAL.json` (`training_approved: true`). Consequences:

- **Do not filter on `training_eligible` / `TRAINING_ELIGIBLE`.** It would drop
  every record. These columns carry no selection signal in this release.
- **Do not treat the stale flags as the quarantine.** The 3,464-record Somatic
  quarantine (§4.2) was applied *before* the parquet was written; the parquet is
  the already-selected release. The flags do not re-express it.
- **Do not "fix" or re-derive the flags.** The release is bound to these exact
  hashes. Any loader must explicitly acknowledge `RELEASE_APPROVAL.json` rather
  than silently bypassing eligibility checks.

### 3.2 `chrX` / `chrY` / `chrM` route to **train** (explicit decision)

The parquet contains 25 contigs: `chr1`–`chr22`, `chrX`, `chrY`, `chrM`. The
chromosome map pins `chr1`→test, `chr21`–`chr22`→val, `chr2`–`chr20`→train. The
remaining **532,281 records (~3.7%)** — mostly `chrX` (≈511K), plus `chrY`
(16,334) and `chrM` (4,452) — and any unrecognized contig are assigned to
**train** by an explicit, recorded decision (§5.3).

This is a deliberate choice, not a silent fallthrough. The approval sidecar had
flagged `"other_contigs": "Not silently admitted to training"`; the recorded
decision (§5.3) resolves it by admitting them to **train**. If you change this
(§5.3), update the split function and this note together.

### 3.3 Reserved samples must never reach training

5 samples are held out entirely for zero-shot / low-VAF / rescue evaluation
(§5.1). Filter them at the **sample** level before any variant-level split.

---

## 4. Label contract

### 4.1 `FILTER` is the training label

The label vocabulary is **frozen** across the pipeline→model contract:

```
{Somatic, Germline, Reference, Artifact, NoConsensus, RNAedit}
```

The approved parquet has already been filtered to the training whitelist — its
`FILTER` column contains exactly three classes:

| Class (`FILTER`) | Records | Role |
|------------------|---------|------|
| `Reference` | 12,791,828 | negative (weak) label |
| `Germline` | 1,522,528 | negative (weak) label |
| `Somatic` | 28,241 | positive label |
| **Total** | **14,342,597** | |

`Artifact`, `NoConsensus`, and `RNAedit` are already excluded — downstream never
sees them. Map `FILTER` directly to your class index. Do not invent a
`PASS` interpretation and do not migrate to `PASS`+`INFO` without a coordinated
multi-repo change.

These are **workflow-derived weak labels**, not independently verified genotype
truth. The negatives are native-evidence Germline/Reference calls; they were
retained **without** a new read-count requirement. Treat class accuracy as
unvalidated (`label_confidence: unvalidated`, §3.1).

### 4.2 Somatic quarantine (already applied)

The approved selection quarantined **3,464 Somatic records** (3,462 for common
population AF, 2 native-class conflicts) as conservative training exclusions.
These are **not proven false positives** and are preserved outside the release for
sensitivity analysis — they are simply absent from `cohort63.review.parquet`
(hence the 28,241 figure, not 31,705). Do not re-add them to training; do not
treat them as negatives.

### 4.3 Cohort exclusions (already applied)

Three samples are excluded entirely and appear in **no** pool:
`PRJNA298330_4032`, `PRJNA298376_4081`, `PRJNA298376_4255`. Re-admission requires
a cleaned label and an explicit manifest change (see `data_split_strategy.md`
Stage 0) — do not add them back ad hoc.

---

## 5. Split strategy and train/val/test method

The split has **two independent levels**. Apply them in order.

```
                    63 approved samples (cohort63.review.tsv)
                                    │
                 Level 1 — SAMPLE pool gate (by sample_id)
                                    │
              ┌─────────────────────┴─────────────────────┐
              ▼                                           ▼
     reserved (5 samples)                       train_pool (58 samples)
     held-out evaluation ONLY                   the only training candidates
     never train/val/test                                │
                                         Level 2 — VARIANT split (by CHROM)
                                                         │
                    ┌───────────────────┬────────────────┴───────────────────┐
                    ▼                   ▼                                      ▼
                  test                val                                    train
                  chr1              chr21–22                chr2–20 + chrX/Y/M + other
             1,288,502            525,937                             11,722,395

        (counts are over the 58 train_pool samples only. chrX/chrY/chrM and any
         unrecognized contig route to train by explicit decision — §5.3.)
```

### 5.1 Level 1 — sample pool (reserved vs train_pool)

| Pool | Samples | Use |
|------|---------|-----|
| `reserved` | 5 | Held-out evaluation **only**: zero-shot (unseen disease), low-VAF/DP, rescue validation. **Never in train/val/test.** |
| `train_pool` | 58 | The only samples eligible for the Level-2 chromosome split. |

Reserved samples (all PASS-era, non-digestive cancers):

```
PRJNA298310_3812   melanoma
PRJNA298376_4166   lung
PRJNA298376_4214   lung cancer
PRJNA298376_4231   ampullary cancer
PRJNA298376_4242   gastric cancer
```

The other 58 are digestive-system cancers and form the training candidates.

**Implementation** — the authoritative membership is the two JSONs (this also
gets you the BAM paths):

```python
import json
from pathlib import Path

ROOT = Path("examples/seq2neo/output_three_class_v2_20260919/cohort63_model_review_20260921")

def load_pools():
    pools = {}
    for pool in ("train_pool", "reserved"):
        d = json.loads((ROOT / f"{pool}.review.json").read_text())
        # NOTE: d["training_approved"] is a stale pre-approval snapshot (§3.1).
        # Approval comes from RELEASE_APPROVAL.json, not this field.
        for s in d["samples"]:
            pools[s["sample_id"]] = {"pool": pool, **s}
    return pools

POOLS = load_pools()               # sample_id -> {pool, tumor_dna, normal_dna, ...}
RESERVED = {sid for sid, v in POOLS.items() if v["pool"] == "reserved"}
```

### 5.2 Level 2 — variant split (chromosome map)

Within `train_pool` samples only, assign each variant by its `CHROM` — the
`deepsomatic` strategy (reference impl `neo_var/src/neo_var/data/split_dataset.py`,
a cross-repo path that **may be obsolete** — see References):

| Split | Chromosomes | Records (train_pool) | Description |
|-------|-------------|---------|-------------|
| **train** | `chr2`–`chr20` **+ `chrX`/`chrY`/`chrM` + any other contig** | 11,722,395 | 19 autosomes plus the non-standard contigs (§5.3) |
| **val** | `chr21`, `chr22` | 525,937 | 2 small autosomes |
| **test** | `chr1` | 1,288,502 | largest autosome; primary metrics |

(Counts over the 58 `train_pool` samples = 13,536,834 variant records. The 5
`reserved` samples add another 805,763 and are excluded from these folds.)

This is a **variant-level, chromosome-blocked** split: the same sample contributes
variants to train, val, and test. That is intentional — it keeps loci
non-overlapping across folds (no position leakage) while using all samples. It is
*not* a sample-level holdout; the sample-level holdout is the reserved pool (§5.1).

Reference implementation:

```python
import polars as pl

TRAIN_CHROMS = {f"chr{i}" for i in range(2, 21)}   # chr2 .. chr20
VAL_CHROMS   = {"chr21", "chr22"}
TEST_CHROMS  = {"chr1"}

def assign_split(chrom: str) -> str:
    """Chromosome-blocked train/val/test. Returns 'train'|'val'|'test'.

    chr1 -> test; chr21/chr22 -> val; **everything else -> train**, which includes
    chr2-20 AND chrX/chrY/chrM AND any unrecognized contig (explicit decision, §5.3).
    """
    if chrom in TEST_CHROMS:  return "test"
    if chrom in VAL_CHROMS:   return "val"
    return "train"            # chr2-20 + chrX/Y/M + other contigs (§5.3)


def split_expr(col: str = "CHROM") -> pl.Expr:
    """Vectorized form of `assign_split` (prefer this over `map_elements`, which
    runs a per-row Python UDF and is slow on 14M rows)."""
    return (
        pl.when(pl.col(col).is_in(list(TEST_CHROMS))).then(pl.lit("test"))
         .when(pl.col(col).is_in(list(VAL_CHROMS))).then(pl.lit("val"))
         .otherwise(pl.lit("train"))   # chr2-20 + chrX/Y/M + other contigs (§5.3)
         .alias("split")
    )
```

### 5.3 Decision: the `other` contigs go to **train** (recorded)

**Decision (recorded 2026-09-22): `chrX`, `chrY`, `chrM`, and any unrecognized
contig are assigned to `train`.** 532,281 records (mostly `chrX`) are affected.

Rationale: sex-chromosome and mitochondrial loci are legitimate training signal,
and folding them into `train` keeps `val`/`test` as clean, small autosomal
holdouts. The approval sidecar had cautioned
`"other_contigs": "Not silently admitted to training"`; this decision is the
**explicit** admission it asked for — it is recorded here rather than left to a
silent `else` fallthrough. The older `data_split_strategy.md` line "Unrecognized
chromosomes … are assigned to train" is consistent with this outcome.

If you need a different policy (e.g. hold `chrX` out, or route it to `test`),
change `assign_split` / `split_expr` (§5.2), the fold-size tables (§6), and this
note together — and record the new decision here.

### 5.4 Cross-check / pre-materialized split (optional)

The older build `data/processed/selected_variants.parquet` **does** carry a
materialized `split` column (`train`/`val`/`test`/`reserved`) plus sub-pool tags
(`is_zero_shot`, `is_low_vaf_a/b`, `is_low_dp`, `is_rescued`, `is_non_rescued`,
`is_indel`) and `label_verdict`. It is built by
`examples/seq2neo/scripts/build_split_manifest.py` over the earlier reconsensus
cohort and is **not** the approved release — but it is a useful reference for the
split semantics and for the reserved-pool sub-pool tags. If you use its tags,
join on `sample_id`+`CHROM`+`POS`+`REF`+`ALT` and re-derive the split from §5.2
rather than trusting its `split` for the approved rows. Design notes:
[`adr/0001-split-manifest-design.md`](adr/0001-split-manifest-design.md).

---

## 6. Training data preparation — recipe

End-to-end, from the approved parquet to train/val/test **variant-level feature
tables**. This stops at the per-variant feature table; **read-level tensor
construction is out of scope** (see the Scope note at the top) and left to each
model implementation.

### Step 0 — verify the release (§1)

Match the SHA-256 pins. Abort if any differ.

### Step 1 — load pools and drop reserved + excluded samples

```python
import polars as pl

parquet = "examples/seq2neo/output_three_class_v2_20260919/handoff_review_20260921/cohort63.review.parquet"

lf = (
    pl.scan_parquet(parquet)
      # Reserved samples are held-out evaluation only — never train/val/test.
      .filter(~pl.col("sample_id").is_in(list(RESERVED)))
      # The 3 excluded samples are already absent; assert it (§4.3).
)
```

### Step 2 — map the label (§4)

```python
# FILTER is the training label. Exactly 3 classes are present.
LABEL = {"Reference": 0, "Germline": 1, "Somatic": 2}
lf = lf.with_columns(pl.col("FILTER").replace_strict(LABEL).alias("label"))
```

> Do **not** filter on `training_eligible`/`TRAINING_ELIGIBLE` (§3.1) — they are
> uniformly stale `False`/`NO` and would empty the dataset.

### Step 3 — assign the variant split (§5.2)

```python
lf = lf.with_columns(split_expr())   # train/val/test; chrX/Y/M + other -> train (§5.3)
# No filtering needed: every variant lands in exactly one of train / val / test.
```

### Step 4 — assemble variant-level features (tensors are the model's job)

Per-variant quantitative features are already aggregated in the parquet:
`DP_DNA_MEAN`, `DP_RNA_MEAN`, `VAF_DNA_MEAN`, `VAF_RNA_MEAN` — use these directly
for a feature vector. For richer per-caller features, parse the `*_BY_CALLER`
columns (encoding: `caller:value` pairs joined by `|`, `.` = missing). This step
yields a flat per-variant feature table; **turning it into read-level tensors is
left to the specific model implementation** (Scope note, top).

```python
def parse_by_caller(s: str) -> dict:
    """'deepsomatic:43,2|strelka:44,2' -> {'deepsomatic': '43,2', 'strelka': '44,2'}"""
    out = {}
    for entry in (s or "").split("|"):
        if entry and entry != "." and ":" in entry:
            k, v = entry.split(":", 1)
            out[k] = None if v == "." else v
    return out
# Apply to GT/DP/AD/VAF_BY_CALLER and the NORMAL_* equivalents. AD is "ref,alt".
```

Collect per fold:

```python
folds = {
    s: lf.filter(pl.col("split") == s).collect()
    for s in ("train", "val", "test")
}
```

### Step 5 — class balance

Class mix is heavily skewed toward `Reference` (≈89%). Class balancing /
class-ratio sampling is a **trainer-side** concern, not a split-layer concern —
keep the split pure and rebalance in the sampler.

### Step 6 — evaluation discipline

- **Primary metrics** on `test` (`chr1`), PASS-quality stratum.
- **Validation** on `val` (`chr21`–`chr22`) for model selection.
- **Reserved pool** (5 samples) reported separately as zero-shot / low-VAF /
  rescue evaluation — never pooled into the primary metrics, and never trained on.

(`chrX`/`chrY`/`chrM` live inside `train` per §5.3, so there is no separate
`other` stratum to report.)

### Expected fold sizes (computed from `cohort63.review.parquet`)

**Training set — `train_pool` (58 samples).** These are the real sizes after the
Step 1 reserved-removal. `train` includes `chrX`/`chrY`/`chrM` + other contigs
per §5.3:

| Fold | Somatic | Germline | Reference | Total |
|------|---------|----------|-----------|-------|
| train (`chr2`–`chr20` + `chrX`/`Y`/`M` + other) | 22,723 | 1,215,023 | 10,484,649 | 11,722,395 |
| val (`chr21`–`chr22`) | 781 | 54,603 | 470,553 | 525,937 |
| test (`chr1`) | 2,383 | 129,156 | 1,156,963 | 1,288,502 |
| **total (58 samples)** | **25,887** | **1,398,782** | **12,112,165** | **13,536,834** |

Class mix in `train` is ≈0.19% Somatic / 10.4% Germline / 89.4% Reference —
heavily skewed; rebalance in the sampler (§6 Step 5).

**Reserved pool (5 samples) — held out, never trained on.** Evaluated as one set
(zero-shot / low-VAF / rescue), **not** split into folds:

| Pool | Somatic | Germline | Reference | Total |
|------|---------|----------|-----------|-------|
| reserved (5 samples) | 2,354 | 123,746 | 679,663 | 805,763 |

**Reconciliation:** 13,536,834 (train_pool) + 805,763 (reserved) = 14,342,597,
matching the release total; class totals 28,241 Somatic / 1,522,528 Germline /
12,791,828 Reference match `RELEASE_APPROVAL.json` exactly.

---

## 7. Provenance & reproducibility

- **Approval**: `RELEASE_APPROVAL.json`, release id
  `seq2neo_separated_three_class_v2_cohort63_20260921`, dated **2026-09-21**,
  authority "User explicit instruction: ok, approve and commit".
- **Policy**: `separated_three_class_v2` (see `candidate_policy` on every row).
- **Verification**: `report.json` (sha `eb6f11ef…`), status
  `review_bridge_verified_not_training_approved` — this is *technical*
  verification of ordered `CHROM/POS/REF/ALT/FILTER` roundtrips and class counts;
  it is deliberately **not** the training approval. The approval is the sidecar.
- **Alignment inventory**: 189 DN/DT/RT BAM paths + indexes found. Availability
  ≠ BAM integrity ≠ sample identity — the inventory checks presence/size only.
- **Historical flags** are immutable pre-approval snapshots (§3.1). Consumers must
  acknowledge this release explicitly; do not globally bypass eligibility checks.

### Conditions attached to the approval (verbatim intent)

1. Never convert uncertain / unselected / `Artifact` / `NoConsensus` / `RNAedit`
   records to `Reference` automatically.
2. Preserve reserved-sample membership and the original results/quarantine for
   sensitivity analysis.
3. The model must use the correct label mapping, split policy, and a
   content-bound (release-specific) feature-cache identity.
4. Approval is **not** independent biological-truth certification, nor evidence
   of destination-model compatibility.

### Feature-cache identity (condition 3) — recipe

Bind any derived feature / tensor cache to this exact release so a stale cache can
never be silently reused against different data. Derive the cache key from the
release id + content hashes, not from a mutable path:

```python
import hashlib, json

RELEASE_ID  = "seq2neo_separated_three_class_v2_cohort63_20260921"
PARQUET_SHA = "8e61c3c73b51cbca095b0da26e34f76b5adfe80c594723acaaedf4d740662a6d"

def feature_cache_id(feature_config: dict) -> str:
    """Content-bound cache key: release + data hash + feature/config hash."""
    h = hashlib.sha256()
    h.update(RELEASE_ID.encode())
    h.update(PARQUET_SHA.encode())
    h.update(json.dumps(feature_config, sort_keys=True).encode())
    return f"{RELEASE_ID}-{h.hexdigest()[:16]}"
```

Include every knob that changes the output in `feature_config` — which columns,
the `parse_by_caller` fields, the label map (§4), and the split policy (§5). Store
the cache under that id so changing the release, the data, or the config forces a
rebuild.

### Known limitations

- Labels are **workflow-derived weak labels**, not verified genotype truth.
  Negatives are native-evidence Germline/Reference, unvalidated.
- `label_confidence` is `unvalidated` across the board.
- Training runs on another host with another user's code (not the local
  EvoSomatic copy); this release is data-only. Feature-cache identity must be
  bound to this release's hashes (recipe above).
- Read-level tensor construction is out of scope here (Scope note, top); it lives
  in each model implementation and must honor this split + label mapping.

---

## 8. Quick reference

| Question | Answer |
|----------|--------|
| Which parquet do I train on? | `handoff_review_20260921/cohort63.review.parquet` |
| Which manifest? | `handoff_review_20260921/cohort63.review.tsv` |
| What makes it "approved"? | `cohort63_handoff_tools_20260921/RELEASE_APPROVAL.json` (hash-bound, `training_approved: true`) |
| What is the label? | the `FILTER` column ∈ {`Somatic`, `Germline`, `Reference`} |
| How do I split train/val/test? | by `CHROM`: test=`chr1`, val=`chr21`–`chr22`, train=`chr2`–`chr20` + `chrX`/`Y`/`M` + other (§5.3) |
| Read-level tensors? | Out of scope — variant-level features only; tensors are each model's job (Scope note) |
| Which samples train? | `train_pool.review.json` (58); `reserved.review.json` (5) are held out |
| Should I filter `training_eligible`? | **No** — stale `False` everywhere (§3.1) |
| Where is the split column? | Not in the approved files — derive from `CHROM` (§5.2); the older `selected_variants.parquet` has one (§5.4) |

---

## References

- Approval sidecar: `…/cohort63_handoff_tools_20260921/RELEASE_APPROVAL.json`
- Handoff narrative: `…/cohort63_handoff_tools_20260921/COHORT_HANDOFF.md`
- Split design & reconsensus-cohort record: [`data_split_strategy.md`](data_split_strategy.md)
- Split-manifest design ADR: [`adr/0001-split-manifest-design.md`](adr/0001-split-manifest-design.md)
- Chromosome-split reference impl *(cross-repo, sibling `neo_var` — **may be
  obsolete**)*: `neo_var/src/neo_var/data/split_dataset.py`
- Older split builder: `examples/seq2neo/scripts/build_split_manifest.py`
- Label/FILTER contract *(cross-repo path, **may be obsolete**)*: `docs/ENSEMBLEVAR_OPTIMIZATION_GUIDELINES.md` §3

> **Cross-repo references** (anything under `neo_var/`, or the
> `ENSEMBLEVAR_OPTIMIZATION_GUIDELINES.md` link) point outside this repository /
> workspace. They are kept as-is for provenance but **may be stale, moved, or
> obsolete** — confirm against the current tree before relying on them. In-repo
> paths (`data_split_strategy.md`, `adr/`, `build_split_manifest.py`) are
> authoritative.
