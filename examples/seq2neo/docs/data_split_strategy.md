# Data Split Strategy: QC Gate + Reserved Downstream Pools + Chromosome Train/Val/Test

Supersedes the pre-reconsensus split plan in `neo_var/docs/data_split_strategy.md`.
Updated for the reconsensus rerun cohort (`output_reconsensus/`), G0 FASTQ MD5
verification, and label_qc verdicts.

## Overview

The tensor extraction dataset is built from the **63-sample QC-passed cohort**
(56 PASS + 7 WARN; 12 diseases) drawn from the 66-sample reconsensus rerun. It is
partitioned into:

1. **Reserved downstream pool** — 5 non-digestive cancer samples held out entirely for
   zero-shot evaluation, low VAF/DP assessment, and rescue validation.
2. **Chromosome-based train/val/test split** — the remaining 58 samples split via the
   `deepsomatic` strategy (chr1=test, chr21–22=val, chr2–20=train).

```
                        66 rerun samples
                               │
                    QC gate (root-cause check)
                               │
              ┌────────────────┴────────────────┐
              ▼                                 ▼
    Excluded (3 samples)                Working cohort (63 samples)
    4032 / 4081 / 4255                  56 PASS + 7 WARN
    RCA concluded — keep excluded               │
                                ┌───────────────┴───────────────┐
                                ▼                               ▼
                      Reserved (5 samples)              Remaining (58 samples)
                      non-digestive, all PASS           digestive-system cancers
                                │                               │
                    ┌───────────┼───────────┐         ┌─────────┼─────────┐
                    ▼           ▼           ▼         ▼         ▼         ▼
                  zero_shot  low_vaf    rescue      Train     Val      Test
                             low_dp   /non_rescued chr2-20  chr21-22   chr1
```

## Stage 0: Sample Selection (QC Gate) — Before Any Split

Three samples are **excluded** from training and evaluation. The root-cause
analysis they were gated on has concluded (`docs/RERUN_LABEL_ANOMALY_ANALYSIS.md`):
4032's diff is a raw-rerun-vs-QC-cleaned artifact, and 4081/4255's clean baselines
predate the 2026-03-18 germline rule, but all three remain TruthQC-suspect samples
whose rerun labels are raw (un-QC'd). Decision: **keep excluded**; re-admission
would require re-running their labels through `label_qc --apply` first.

| sample_id | label_qc verdict | set | disease | status |
|-----------|-----------------|-----|---------|--------|
| PRJNA298330_4032 | FAIL | — | — | excluded (RCA concluded) |
| PRJNA298376_4081 | PASS | 1 | colorectal cancer | excluded (RCA concluded) |
| PRJNA298376_4255 | WARN | 1 | colorectal cancer | excluded (RCA concluded) |

Note: 4081 kept a PASS gate verdict and 4255 a WARN verdict — the exclusion is the
RCA decision layered on top of the gate, recorded as `status=useless` in
`sample_manifest_rerun.tsv`. The gate verdict alone is not the inclusion rule.

- **Training cohort**: the other 63 samples (56 PASS + 7 WARN) — exactly the rows
  with a non-empty `training_label_vcf` in `sample_manifest_rerun.tsv`.
- **Re-admission**: only via a cleaned rerun label (`label_qc --apply`) and an
  explicit manifest change; until then these samples appear in no split and no
  reserved pool.
- The remaining 7 WARN samples stay in the cohort but are **train-only** (see
  `label_verdict` rule below).

## Label Provenance

Labels come from the reconsensus rerun, not the original pipeline VCFs. Source of
truth: `examples/seq2neo/data/processed/sample_manifest_rerun.tsv`
(`training_label_vcf` column — the resolved per-sample label; empty for the 3
excluded samples). It resolves the verdict-dependent choice:

- **PASS** → `output_reconsensus/<sample>/rescue/...rescued...filtered.vcf.stripped.vcf.gz`
- **WARN** → `runs/label_qc/cohort66_apply/cleaned_vcf/<sample>.bgz.vcf.gz`

(`training_manifest_65.tsv` encoded the same rule in a separate 65-row file and is
**superseded** by the rerun manifest; kept in the repo for history only. All paths
above live under the rsynced tree at
`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo`.)

Cohort-wide label caveats (manifest `known_limitations` column): labels are
nuclear-only (chrM records were dropped by the gnomAD scatter-gather; fixed in
`d336c7f` after the cohort — ~500 records/sample incl. ~3 rescue-stage Somatics)
and predate the Rule-2 common-AF/prior-artifact vetoes (zero Rule-2 firings
measured cohort-wide, so no practical impact).

The FILTER vocabulary is frozen (`{Somatic, Germline, Reference, Artifact,
NoConsensus, RNAedit}`); extraction keeps the whitelist `{Somatic, Germline,
Reference}`.

### `label_verdict` column

`split_assignments.parquet` gains a `label_verdict` column (`PASS` / `WARN`).

- **Evaluation pools — PASS only**: reserved sub-pools and the test split are built
  exclusively from PASS-verdict samples.
- **WARN → train**: the 7 remaining WARN samples (PRJNA298330_3948, PRJNA298376_3812,
  PRJNA298376_3948, PRJNA298376_3978, PRJNA298376_3995, PRJNA298376_4112,
  PRJNA298376_4232) flow to train only, where cleaned-label noise does the least harm.
- All 5 reserved samples are PASS, so this rule costs nothing in the reserved pool.

## Data Provenance (MD5 Pinning)

Every input FASTQ is pinned to the G0 checksum contract:

- `examples/seq2neo/data/processed/merged.json` — per-pair `r1_md5` / `r2_md5`
  (from `neo_gate/contracts/G0_FASTQ_MD5_RESULTS.tsv`, `status=match` rows only);
  `summary.md5_source` records provenance.
- `examples/seq2neo/data/raw/md5sums.txt` — `md5sum -c`-compatible manifest over the
  `data/raw/<PRJ>/<patient>/<modality>/` symlink tree (638 files, 0 broken links).

`split_assignments.parquet` (or its sidecar metadata) should reference this manifest
so the exact FASTQ bytes behind every label are reproducible.

## Cross-Project Patient Identity — Option A (Independent)

Six numeric patient IDs appear in more than one project (3812 in all three;
3942, 3948, 3978, 3995, 4007 in PRJNA298330 + PRJNA298376). **Decision: treat them
as independent individuals (Option A).** Evidence:

- MD5 checksums prove the FASTQ files are **distinct sequencing runs**, not
  byte-identical redepositions (e.g. 3812 DN: `SRR2634745` vs `SRR2635083` vs
  `SRR9697858`, all different hashes).
- The G0 cohort manifest adjudicates them `standalone_participant_surrogate` with
  distinct biosamples (`SAMN04158xxx` for PRJNA298330 vs `SAMN15374xxx` for
  PRJNA298376) and `duplicate_policy` of `resequencing_redeposition_candidate` /
  `newer_specimen_observation`.

At worst these are related-but-resequenced individuals — no duplicated tensors.
No group-level partitioning by numeric ID is applied. (The rejected conservative
Option B — grouping by numeric ID — would have pulled PASS `PRJNA298330_3812` and
WARN `PRJNA298376_3812` into the reserved pool alongside `PRJNA298310_3812`,
colliding with the PASS-only evaluation rule.)

## Reserved Downstream Pool

### Sample Selection

5 samples from non-digestive cancers, all present in the rerun cohort with PASS
verdicts:

| sample_id | disease_normalized | set | verdict | rationale |
|-----------|-------------------|-----|---------|-----------|
| PRJNA298310_3812 | melanoma | 3 | PASS | skin cancer |
| PRJNA298376_4166 | tumor dna sample from lung... | 4 | PASS | lung tissue |
| PRJNA298376_4214 | lung cancer | 4 | PASS | lung cancer |
| PRJNA298376_4231 | ampullary cancer | 3 | PASS | rare, non-colorectal |
| PRJNA298376_4242 | gastric cancer | 4 | PASS | stomach, non-colorectal |

All other 58 working-cohort samples are digestive-system cancers (colon, colorectal,
pancreatic, rectal, cholangiocarcinoma, esophageal, bile duct) and flow into the
chromosome split.

### Variant Filtering

Only variants with FILTER in `{Somatic, Germline, Reference}` are included.
`NoConsensus`, `Artifact`, and `RNAedit` are excluded from the split assignments.

### Sub-Pool Tags (Boolean, Overlapping)

Within the reserved pool, each variant carries boolean tags for flexible downstream use:

| Tag | Criteria | Description |
|-----|----------|-------------|
| `is_zero_shot` | all reserved variants | unseen disease evaluation |
| `is_low_vaf_a` | DNA_VAF_mean < 0.15 AND RNA non-zero | low allele-fraction calls with RNA evidence |
| `is_low_vaf_b` | 0.15 ≤ DNA_VAF_mean < 0.30 AND RNA non-zero | moderate-low allele-fraction calls with RNA evidence |
| `is_low_dp` | BAM_DT_DP < 20 AND RNA non-zero | low tumor depth calls with RNA evidence |
| `is_rescued` | RESCUED = "YES" | cross-modality rescued variants |
| `is_non_rescued` | RESCUED = "NO" | non-rescued variants |
| `is_indel` | variant_type IN ("INS", "DEL") | insertion and deletion variants |

**RNA non-zero** is defined as: `(RNA_VAF_mean > 0) AND (RNA_DP_mean > 0)`.

Tags are **not mutually exclusive** — a variant can be both `is_low_vaf_a` and
`is_rescued`. The `is_low_vaf_a` and `is_low_vaf_b` pools are disjoint by
definition (VAF ranges do not overlap).

### Expected Counts

**Stale — re-measure from rerun VCFs.** The previous tables (277k zero-shot,
~60k rescued, 8.6k INDEL, etc.) were computed on the pre-rerun 61-sample label set
and must not be quoted for the new cohort.

## Chromosome Train/Val/Test Split

The remaining 58 samples are split using the **deepsomatic** strategy
(`neo_var/src/neo_var/data/split_dataset.py`), unchanged:

| Split | Chromosomes | Description |
|-------|-------------|-------------|
| Train | chr2–chr20 | 19 autosomes, majority of data |
| Val | chr21–chr22 | 2 small autosomes |
| Test | chr1 | largest autosome, **PASS-verdict samples only** |

Unrecognized chromosomes (e.g., `chrUn_*`, scaffolds) are assigned to train.

### Expected Counts

**Stale — re-measure.** Previous split totals (train ~5.86M / val ~286k /
test ~646k over 56 samples) predate the rerun labels and the 58-sample pool.

## Implementation Notes

- Reserved-sample identification and tag logic are unchanged from
  `neo_var` `reserve_downstream.py`; inputs change to:
  - manifest: `examples/seq2neo/data/processed/sample_manifest_rerun.tsv`
    (rows with non-empty `training_label_vcf`) instead of `sample_manifest.parquet`
  - labels: the manifest's `training_label_vcf` paths per the provenance rules above
- Output schema `split_assignments.parquet` gains `label_verdict` (str:
  `PASS` / `WARN`; from the manifest's `label_qc_verdict`) alongside the
  existing columns.
- Excluded samples (4032/4081/4255) re-join only via the Stage 0 re-admission
  rule; if ever re-admitted, they join at their manifest `set_number` with the
  same assignment rules (evaluation pools only if PASS).

## Design Decisions

1. **QC gate before split, with RCA-concluded exclusions.** 4032/4081/4255 are
   excluded per the concluded root-cause analysis; re-admission stays explicit and
   cheap because the gate sits upstream of all split assignments.

2. **5 samples reserved, not by set.** The 4 disease-based sets don't cleanly
   separate non-digestive from digestive cancers; per-sample reservation gives
   precise control. The reserved 5 span sets 3 and 4 and are all PASS.

3. **Evaluation on PASS labels only.** WARN labels went through label_qc cleaning;
   fine for training signal, not for benchmark ground truth.

4. **Option A identity handling.** MD5 evidence shows cross-project same-ID samples
   are distinct data; treating them as independent avoids shrinking the cohort on an
   unproven identity hypothesis.

5. **Overlapping boolean tags, not mutually exclusive categories.** Boolean columns
   let downstream tasks compose subsets flexibly without losing information.

6. **Deepsomatic chromosome split reused.** The existing `split_dataset.py` code and
   its validation logic are unchanged; only the reserving/QC layer above it is new.

## References

- Rerun labels: `examples/seq2neo/output_reconsensus/`, `data/processed/sample_manifest_rerun.tsv`
  (`training_manifest_65.tsv` superseded)
- FASTQ integrity: `neo_gate/contracts/G0_FASTQ_MD5_RESULTS.tsv`,
  `examples/seq2neo/data/raw/md5sums.txt`, `merged.json` `r1_md5`/`r2_md5`
- Identity adjudication: `neo_gate/contracts/G0_COHORT_MANIFEST.draft.csv`
- Chromosome split: `neo_var/src/neo_var/data/split_dataset.py`
- Prior version of this strategy: `neo_var/docs/data_split_strategy.md`
- Re-consensus runbook: `docs/RECONSENSUS_RERUN.md`
