# Rerun vs Clean Label Anomaly Analysis (4032 / 4081 / 4255)

Date: 2026-09-01. Branch: `stats`. Related: `docs/RECONSENSUS_RERUN.md`, `docs/ENSEMBLEVAR_OPTIMIZATION_GUIDELINES.md`.

## Scope

Per-sample comparison of the 2026-08/09 re-consensus rerun labels
(`examples/seq2neo/output_reconsensus/<SAMPLE>/rescue/*_rescued_*/*.filtered.vcf.stripped.vcf.gz`)
against the previous QC-cleaned truth VCFs
(`truth_qc_out/cleaned_vcf/`, TruthQC run `20260826_204134`), 66 samples, using
`examples/seq2neo/scripts/compare_vcf_counts.py` (totals + per-FILTER).

Result: 63/66 samples share one consistent relabeling signature
(Germline +25k…+50k, offset by −NoConsensus; Somatic/Reference/RNAedit nearly
flat; totals change ≤0.01%). Three samples deviate and were investigated here:

| sample | total Δ | signature |
|---|---|---|
| PRJNA298330_4032 | +4,288 | Somatic 2,826→7,115 (+4,289) |
| PRJNA298376_4081 | +10,205 | **inverted**: Germline −30,261, NoConsensus +26,239, Reference +9,569, Somatic −546 |
| PRJNA298376_4255 | +19,211 | **inverted**: Germline −43,044, NoConsensus +37,756, Reference +20,820, Somatic −576 |

## Findings

Both anomalies trace to the **clean baselines**, not to rerun inputs or config:

1. **4032 — raw-vs-QC'd comparison artifact.** The clean VCFs went through
   `label_qc --apply`; the rerun outputs are raw. All 4,288 sites TruthQC
   dropped for 4032 reappear in the rerun (4,285 as Somatic, 3 as Germline).
   These records are RNA-only (`rna_consensus_only`, `RESCUE_PROMOTED=NO`,
   zero DNA support, `PASSES_CONSENSUS_RNA=YES`) at common gnomAD AF
   (93.6% with AF≥1%) — exactly the R1+R2+R3 pattern the QC gate drops.
   The +4,289 Somatic delta is not a rescue-path bug.

2. **4081 / 4255 — code vintage of the clean baselines.** The germline rule in
   `bin/vcf_utils/variant_classifier.py` changed 2026-03-18 (`92ada18`):
   pre-Mar-18, any site with gnomAD AF above threshold got FILTER=Germline
   with **no caller-evidence requirement**. Of the 15 short-named clean files,
   exactly 4081 (Mar 6) and 4255 (Mar 9) predate that change; all others
   postdate it — precisely the inverted pair. At every flipped site the
   aggregated caller evidence (`FILTERS_ORIGINAL`) is byte-identical between
   clean and rerun; effective rerun params differ from normal samples only in
   `input`/`outdir`. The rerun (M6-era rule, `variant_classifier.py:490`)
   correctly refuses Germline where the only evidence is a DeepSomatic
   `RefCall` or RNA-only `GERMLINE`/`RefCall` (`insufficient_modality_support`
   → NoConsensus). Their Somatic losses are justified demotions of
   self-contradictory labels (COSMIC>0 but no Somatic-labeled caller;
   `UNIFIED_FILTER=NoConsensus` even in the old files).

3. **These samples are intrinsically abnormal anyway.** RefCall-heavy
   DeepSomatic DNA inputs (4081: 4.88M records, 4.84M RefCall, vs 4060's
   0.89M/0.86M) and TruthQC already flagged exactly these three as
   疑似异常样本 (Somatic drop rates 60–93% vs ~3% typical).

## Verdicts

- **4081, 4255**: rerun labels are correct and strictly preferable to the old
  baselines; no re-run needed (inputs and code confirmed identical to the
  cohort). Keep excluded from training on biological grounds, as TruthQC
  advised. If re-admitted, use the rerun VCF gated by a fresh
  `label_qc --apply`, never the old cleaned files.
- **4032**: exclusion from the training manifest (commit `50452ee`) is the
  right call. Its rerun VCF needs `label_qc --apply` before any use, like any
  other sample.

## Follow-up gap

The consensus-stage `rna_consensus_only` rule emits Somatic at common-AF sites
with zero DNA support (the 4032 pattern, ~QC-drop-count records cohort-wide).
The M6/`d336c7f` common-AF vetoes only guard annotation-stage reclassification;
today only `label_qc` catches these downstream. A common-AF veto at consensus
stage would close the gap.

## Evidence artifacts

- Count tables: `tests/results/rerun_vs_clean_counts.tsv`,
  `tests/results/rerun_vs_clean_counts_by_filter.tsv`,
  `tests/results/rerun_vs_clean_per_filter_report.md`
- Analysis scripts + per-sample JSONs: `tests/debug_tools/label_anomaly/`
  (`transitions.py` — FILTER transition matrix + rerun-side INFO breakdown;
  `vc_compare.py` — clean-vs-rerun caller-level classification comparison;
  `cosmic_rescue_check.py`, `qc_drop_check_4032.py` — one-off checks).

## Caveats

- Caller-input identity was verified at flipped sites via `FILTERS_ORIGINAL`
  and via checksum manifests (`runs/rerun_checksums/`), not by whole-file
  comparison against the March originals (not retained).
- Clean VCFs lack `CLASSIFICATION_RATIONALE`/`RESCUE_PROMOTED`/
  `ALT_COUNT_BY_CALLER` (newer output-contract fields), so clean-side rule
  attribution relied on `FILTERS_ORIGINAL`/`UNIFIED_FILTER*`/rescue flags.
