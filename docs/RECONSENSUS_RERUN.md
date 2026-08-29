# Consensus+Rescue-Only Rerun Path

## Overview

The rerun path regenerates **consensus and rescue VCFs** for the seq2neo cohort from the **existing, read-only per-caller VCFs** (Mutect2 / Strelka2 / DeepSomatic, DNA branch and realigned RNA branch), writing results to a **new, separate output location**. It exists to re-label the cohort with fixed consensus/rescue code without re-paying alignment and variant-calling compute (see `dev_docs/audit/2026-08-28_adversarial_review.md`, finding M5, and ticket `.scratch/consensus-label-quality/issues/08-rerun-path.md`).

## Guarantees

1. **Consensus+rescue only.** The rerun enters at `--step consensus` with `--tools consensus,rescue,filtering,vep` — no caller names in `--tools`, and the generated samplesheet contains VCF rows only (no FASTQ/BAM/CRAM columns). FASTQ→BAM alignment (`BAM_ALIGN`, gated on `params.step == 'mapping'`) and per-caller variant calling (modules gated on their caller name in `--tools`) are structurally unreachable. The driver additionally refuses any config that names a caller in `tools` or uses a step other than `consensus`.
2. **Normalization guaranteed (M5 fix).** Entering at a post-calling step that feeds consensus now routes caller VCFs through `VCF_NORMALIZE` (VT decompose + `bcftools norm --multiallelics -any --check-ref w`) before consensus. Previously, steps `consensus`/`annotate`/`filtering`/`rna_filtering` triggered `VCF_CONSENSUS_WORKFLOW` but skipped `VCF_NORMALIZE`, silently breaking indel consensus via exact-key (`chrom:pos:ref:alt`) matching.
3. **Read-only inputs + new outdir.** Input VCFs are opened read-only (md5 hashing only). Before each run the driver records a checksum manifest (md5 + size for every input VCF and its tabix index); after the run it re-verifies — a mismatch fails the sample and blocks retry. Outputs land under `output_reconsensus/<sample_id>/`, and the driver refuses an outdir inside any source output root.

## The M5 wiring fix

`subworkflows/local/vcf_normalize/main.nf` previously gated normalization on steps `mapping..norm` only. The fix:

- `subworkflows/local/vcf_normalize/main.nf` — the normalization gate now also covers steps `consensus`, `annotate`, `filtering`, `rna_filtering` when `consensus` is in `--tools`; `input_sample` is routed in for `norm`/`consensus` entries.
- `subworkflows/local/bam_variant_calling/main.nf` — `--step consensus` routes the samplesheet VCFs into `vcf_to_normalise` (previously only `annotate`/`norm`).
- `subworkflows/local/samplesheet_to_channel/main.nf` — VCF samplesheet rows are accepted at `--step consensus` (same handling as `norm`; sample IDs already in `TUMOR_vs_NORMAL` form are used as-is).
- `subworkflows/local/utils_nfcore_rnadnavar_pipeline/main.nf` — `retrieveInput` resolves `--step consensus` without `--input` to `<outdir>/csv/variantcalled.csv`.
- `conf/modules/prepare_resources/prepare_genome.config` — the sequence dictionary build stays skipped at the `consensus` entry (as for `annotate`): nothing in the normalize/consensus/rescue chain consumes it — `VCF_NORMALIZE` takes only `fasta`/`fasta_fai`, and `VCF_CONSENSUS`/rescue are Python scripts; all `dict` consumers are BAM-based GATK subworkflows unreachable from a VCF-only consensus entry. The fasta index (`SAMTOOLS_FAIDX`) is **no longer** skipped at `consensus`: it is built when `fasta_fai` is absent, because `VCF_NORMALIZE`'s `bcftools norm -f` needs it.

**Deliberate behavior change (M5).** The normalization gate now covers all post-calling, consensus-feeding entry steps — `consensus`, `annotate`, `filtering`, and `rna_filtering` — not just the mapping-through-`norm` path. These entry points previously fed raw caller VCFs into exact-key consensus matching; they now deliberately decompose/left-align first, so their outputs differ from before the fix.

The normal `--step mapping` path is unchanged: the gate only adds steps, and normalization still runs exactly once per caller VCF.

## Driver usage

Driver: `examples/seq2neo/scripts/run_reconsensus_rerun.py`, config: `examples/seq2neo/config/rerun.yaml` (same style as `run_batch_from_json.py` / `runner.yaml`).

```bash
# Dry-run: print samplesheets, commands, and input md5s without writing anything
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --dry-run

# Run one sample / one set / the whole cohort (sequential, with retry + resume)
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --sample PRJNA298330_4032
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --set 2
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml

# Re-verify recorded input checksums at any time (exit 1 on any mismatch)
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --verify-only
```

Per sample the driver:

1. Locates the 6 raw caller VCFs under `<base_output_dir>/<dir_name>/` (`variant_calling/<caller>/<DT pair>/` for DNA, `vcf_realignment/variant_calling/<caller>/<RT_realign pair>/` for RNA) using the sample list from `data/processed/sample_manifest.tsv` (`build_sample_manifest.py`).
2. Writes a VCF-only samplesheet to `runs/rerun_csv/<sample>.csv` (columns `patient,sample,status,variantcaller,vcf`; status 1 = DNA pair, 2 = realigned RNA pair).
3. Records input checksums to `runs/rerun_checksums/<sample>.input_checksums.json`.
4. Runs nextflow with `--step consensus --tools consensus,rescue,filtering,vep --outdir output_reconsensus/<sample>` (plus `-resume -offline -with-conda`).
5. Checks completion (execution trace + `consensus/**/*.vcf.gz` + `rescue/**/*.filtered.vcf.gz|*.rescued.vcf.gz`) and verifies input checksums are unchanged.

State is tracked in `runs/rerun_state.json`; re-running the same command skips completed samples and retries failures up to `max_retries`. A tripped checksum guard is never auto-retried.

## Verification status

- `nextflow config .` loads cleanly (25.10.2), and a `-stub-run` at `--step consensus` on a synthetic 6-VCF samplesheet submits exactly `VT_DECOMPOSE`×6 → `BCFTOOLS_NORM`×6 → `VCF_CONSENSUS`×2 (DNA + RNA groups) and no alignment or variant-calling processes. Rescue/post-processing processes sit behind the same `tools`-contains-`rescue` gate as the production mapping path; they were not reached in stub mode because the local modules define no `stub:` blocks.
- Driver behavior (samplesheet layout, dry-run, checksum guard, config guards) is covered by `tests/rerun_driver/` (`.venv/bin/python -m pytest tests/rerun_driver/`).
- Smoke reruns on real cohort data (ticket 12) done on **both** smoke samples: `PRJNA298330_4032` (abnormal class — correctly FAILs) and `PRJNA298376_4278` (clean class — WARN on S7 only) — see the two smoke sections below.

## Smoke run: PRJNA298330_4032 (2026-08-29)

Run via `python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --sample PRJNA298330_4032`; state `succeeded` in `runs/rerun_state.json` (finished 04:40, wall ~1h45m), outputs in `output_reconsensus/PRJNA298330_4032/`.

**Structural checks (all pass):**

- Execution trace `output_reconsensus/PRJNA298330_4032/pipeline_info/execution_trace_2026-08-29_02-55-24.txt`: 49 tasks, **all COMPLETED** — `VT_DECOMPOSE`×6, `BCFTOOLS_NORM`×6 (the M5 fix visibly active), `VCF_CONSENSUS`×2, `VCF_FILTER`×2, `VCF_RESCUE`, COSMIC/gnomAD + RNA-editing annotation, `VCF_RESCUE_FILTER`, VEP, MultiQC, plus QC/interval helpers. **Zero alignment or variant-calling tasks** (no BWA/STAR/MarkDuplicates/BQSR/Mutect2/Strelka/DeepSomatic — grep-verified).
- Input checksums: `--verify-only` reports `[OK] PRJNA298330_4032 inputs untouched` — the six read-only caller VCFs and their indexes are byte-identical after the run.

**label_qc Tier A verdicts** (dry-run, artifacts under `runs/label_qc/`):

| output | somatic | HIGH | actioned | RNA-only common-AF (S3) | self-contradiction (S2) | verdict |
|---|---|---|---|---|---|---|
| old (Apr 10, pre-fix pipeline) | 7,114 | 4,288 | 60.28% | 60.01% | 0.014% (R6×1) | **FAIL** (S0,S1,S3) |
| rerun smoke (this path) | 7,115 | 4,288 | 60.27% | 60.00% | 0.0% | **FAIL** (S0,S1,S3) + S7 WARN |

Ti/Tv (2.053), indel fraction (2.25%), and DNA participation (0.567) are identical to four significant figures. **The 4032 abnormality is reproduced almost exactly by the fixed pipeline, i.e. it is intrinsic to this sample's input caller VCFs (RNA-only germline leakage), not a consensus-labeling artifact** — the gate correctly keeps it out of the training set.

Notable quantitative shifts old → rerun, consistent with the post-audit fixes now active in the rerun code:

- `median_vaf` 0.239 → 0.500 and het-like-VAF fraction 20.6% → 47.2%: the C1 tumor-sample-aware genotype-extraction fix changes which sample's AD/DP feed the VAF; a 0.5 median is exactly what common-AF germline-het leakage predicts.
- `median_dp` 32.7 → 13.0, low-DP fraction 2.2% → 44.5% (new S7 WARN): same root cause (DP now from the tumor genotype fields).
- R9 (low caller agreement, N_SUPPORT_CALLERS ≤ 1) 0 → 2,263 hits, LOW severity only (report, not actioned): consistent with the corrected support counting; does not change the verdict.

**Tier B (BAM verification, `runs/label_qc/smoke_PRJNA298330_4032_tierB/`, dry-run, 364 s):** verdict unchanged — **FAIL** (S0,S1,S3; S7 WARN). Tier B adds B1 (tumor strand bias) 86 LOW-tier sites; S4 normal contamination is clean: only 9/7,112 evaluable Somatic sites (0.1%) have normal alt-VAF ≥ 0.05. So unlike 4081/4255 (78–83% normal contamination), 4032's failure is purely the RNA-only germline-leakage signature, not contamination.

**Driver fix found by the second smoke sample.** `PRJNA298376_4278` failed at samplesheet validation: its `vcf_prefix` is numeric-only (`4278`), and nf-schema coerces the CSV `patient` field to an integer, failing the patient-is-string check. 15 of 66 cohort samples have numeric-only prefixes (all in PRJNA298376, incl. 4081 and 4278). Fixed in `run_reconsensus_rerun.py` — `patient_column()` falls back to `sample_id` when `vcf_prefix` is all digits (a no-op for the other 51, where the two are identical); regression test `test_samplesheet_numeric_vcf_prefix_uses_sample_id_as_patient` in `tests/rerun_driver/` (10/10 pass).

## Smoke run: PRJNA298376_4278 (2026-08-29, clean-sample control)

Run via the same driver command with `--sample PRJNA298376_4278` after the patient-column fix; state `succeeded` (finished 08:09, wall ~3h — VEP dominated; ~1h50m to rescue filtering, then ~1h20m VEP + MultiQC). Trace `output_reconsensus/PRJNA298376_4278/pipeline_info/execution_trace_2026-08-29_05-10-29.txt`: 49 tasks (46 COMPLETED + 3 CACHED from the shared work dir), **zero alignment/calling tasks** (same grep as §4032). Input checksums verified unchanged (`[OK] ... input checksums verified unchanged`).

**label_qc Tier A, old vs rerun** (dry-run; artifacts `runs/label_qc/smoke_PRJNA298376_4278/` and `runs/label_qc/old_PRJNA298376_4278/`):

| output | somatic | HIGH | MID | actioned | RNA-only common-AF (S3) | self-contradiction (S2) | Ti/Tv | low-DP frac (S7) | verdict |
|---|---|---|---|---|---|---|---|---|---|
| old (Apr 28, pre-fix) | 2,384 | 53 | 2 | 2.31% | 1.76% | 0.0% | 2.240 | 3.6% | **PASS** |
| rerun smoke | 2,380 | 50 | 4 | 2.27% | 1.76% | 0.0% | 2.239 | 31.6% | **WARN** (S7 only) |

Every label-quality signature is reproduced on the clean sample — the fixed pipeline neither inflates nor deflates labels. The only delta is again the low-DP fraction (3.6% → 31.6%; for 4032: 2.2% → 44.5%), a systematic consequence of the C1 fix (DP now read from the tumor genotype fields rather than the corrupted pre-fix extraction), not label corruption. **Follow-up: S7's 30% WARN threshold was calibrated on pre-fix DP values and should be recalibrated against the fixed DP fields; until then expect cohort-wide S7 WARNs.** Under the inclusion contract below, WARN means include-after-cleaning — for 4278 that drops 50 HIGH + relabels 4 MID of 2,380 records, which is the gate working as intended.

## Cohort runbook — all 66 samples + training gate

All commands assume the repo root as CWD unless noted. The cohort is 66 samples (`data/processed/sample_manifest.tsv`, all `is_complete=True`), partitioned into 4 disease-exclusive sets (15/20/16/15).

### 1. Rerun the cohort (consensus+rescue only)

**Parallel cohort run (the way the 64-sample remainder is being run, launched 2026-08-29).**
Launcher: `examples/seq2neo/scripts/run_reconsensus_cohort.py`, config: `examples/seq2neo/config/rerun_cohort.yaml`.

```bash
cd examples/seq2neo

# show the partition + exact per-group commands; launch nothing
python3 scripts/run_reconsensus_cohort.py --dry-run

# launch 6 detached group drivers (nohup-style; safe to close the shell)
python3 scripts/run_reconsensus_cohort.py

# re-run the same command later: groups with a live driver are skipped,
# completed samples are auto-skipped inside each group
```

The launcher reads the manifest (66 samples), excludes samples marked `succeeded` in `runs/rerun_state.json` (the 2 smoke samples), and round-robin-partitions the remaining 64 into 6 disjoint groups (11/11/11/11/10/10). Each group gets a detached `run_reconsensus_rerun.py` process running its samples sequentially with a comma-separated `--sample` list.

**No-VEP decision (cohort config).** `rerun_cohort.yaml` is `rerun.yaml` with `tools: consensus,rescue,filtering` — VEP dropped. Rationale: the training-label artifact (the rescue filtered stripped VCF, `rescue/**/<pair>.filtered.vcf.stripped.vcf.gz`) is produced *before* VEP in the process chain — smoke-run execution trace: `... VCF_RESCUE → COSMIC_GNOMAD_ANNOTATION → RNA_EDITING_ANNOTATION → VCF_RESCUE_FILTER → ENSEMBLVEP_VEP → MULTIQC`. VEP cost ~1h20m of the ~3h smoke wall time on 4278 and feeds nothing the label gate consumes. The mandatory `vep_cache` param stays set in `seq2neo.shared.config`, so pipeline validation still passes. Expected duration: ~1.7h/sample without VEP; 6-way parallel → **~18–20h total** for the 64-sample remainder.

**State-file race handling.** The driver's state file is loaded once at startup and written back wholesale — concurrent drivers sharing it would lose each other's updates (last writer wins). The launcher therefore gives each group its own state file via the driver's `--state-file` flag (`runs/cohort_state/group<N>.json`), and each group driver runs with its own CWD (`runs/cohort_work/group<N>/`) so nextflow `work/` and `.nextflow/` session dirs are per-group. Per-sample artifacts (`runs/rerun_csv/<sid>.csv`, `runs/rerun_checksums/<sid>.json`, `output_reconsensus/<sid>/`) are disjoint by the partition.

**Monitoring.**

```bash
cd examples/seq2neo

# per-group progress + per-sample OK/FAIL lines
tail -f runs/cohort_logs/group1.log

# per-sample state (running/succeeded/failed + retry counts)
cat runs/cohort_state/group1.json

# which nextflow pipelines are alive
ps -o pid,etime,cmd -p $(cat runs/cohort_logs/group*.pid)

# task-level detail for a running sample
grep -E 'VCF_CONSENSUS|VCF_RESCUE' \
  output_reconsensus/<sid>/pipeline_info/execution_trace_*.txt | tail
```

A group is done when its log ends with `=== Done succeeded=N skipped=M failed=K ===` and its driver PID exits. Completion of the whole cohort: all 6 group state files plus `runs/rerun_state.json` show 66 `succeeded` (failures retry up to `max_retries=2` on launcher re-run; a tripped checksum guard is never auto-retried).

**Sequential fallback (original driver, includes VEP via `rerun.yaml`):**

```bash
cd examples/seq2neo

# preflight: inspect samplesheets/commands/input md5s without writing anything
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --dry-run | less

# full cohort (sequential; per-sample retry up to max_retries=2; -resume)
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml

# ... or one cross-validation set at a time
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --set 1
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --set 2
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --set 3
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --set 4

# smoke first if desired (4032 done; 4278 pending — see §6)
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --sample PRJNA298376_4278
```

Re-running the same command skips `succeeded` samples (`runs/rerun_state.json`) and retries failures. A tripped checksum guard is never auto-retried — investigate before clearing state.

### 2. Verify every input stayed read-only

```bash
cd examples/seq2neo
python3 scripts/run_reconsensus_rerun.py --config config/rerun.yaml --verify-only
# must end: failed=0  no-manifest=0
```

### 3. Structural check per sample (no alignment/calling)

```bash
cd examples/seq2neo
for t in output_reconsensus/*/pipeline_info/execution_trace_*.txt; do
  grep -Ei 'BWA|STAR|MARKDUPLICATES|BQSR|SPLITNCIGAR|GATK4_MUTECT2|STRELKA_SOMATIC|DEEPSOMATIC' "$t" \
    && echo "CALLER TASKS FOUND: $t"
done
# expected: no output
```

### 4. Build the label_qc cohort manifest (points at the NEW outputs)

```bash
cd examples/seq2neo
{
  printf 'sample_id\ttruth_vcf\tnormal_bam\ttumor_bam\n'
  tail -n +2 data/processed/sample_manifest.tsv | \
  while IFS=$'\t' read -r sid proj pid setn dis disn st reason base dir prefix rvcf complete bdn bdt brt rest; do
    vcf=$(ls output_reconsensus/"$dir"/rescue/*/*.filtered.vcf.gz 2>/dev/null | head -1)
    if [ -n "$vcf" ]; then
      printf '%s\t%s\t%s\t%s\n' "$sid" "$PWD/$vcf" "$bdn" "$bdt"
    else
      echo "MISSING rescue VCF: $sid" >&2
    fi
  done
} > runs/label_qc/cohort_samples.tsv
# expect 66 data rows; any MISSING line means that sample's rerun did not complete
```

### 5. Gate the cohort with label_qc

```bash
# back at repo root; .venv python, pure stdlib + samtools for Tier B
cd /t9k/mnt/hdd/work/Vax/pipeline/rnadnavar

# Tier A (VCF-only; cohort mode n=66 activates the z-score outlier gates S1/S5)
.venv/bin/python bin/label_qc.py \
  --samples examples/seq2neo/runs/label_qc/cohort_samples.tsv \
  --out examples/seq2neo/runs/label_qc/cohort_tierA

# Tier B (optional, slower: batched samtools mpileup on DN/DT BAMs;
# enables S4 normal-contamination and B1/S8 strand-bias gates)
.venv/bin/python bin/label_qc.py \
  --samples examples/seq2neo/runs/label_qc/cohort_samples.tsv \
  --out examples/seq2neo/runs/label_qc/cohort_tierB \
  --verify-bam --samtools /t9k/mnt/joey/bio_gizmo/samtools

# apply cleaning (HIGH dropped, MID relabelled) — writes cleaned VCFs,
# never modifies inputs
.venv/bin/python bin/label_qc.py \
  --samples examples/seq2neo/runs/label_qc/cohort_samples.tsv \
  --out examples/seq2neo/runs/label_qc/cohort_apply --apply
```

### 6. Training-inclusion contract (PASS/WARN/FAIL)

Verdicts live in `<out>/samples_qc.tsv` (per-sample) and `<out>/report.md`:

| verdict | training action |
|---|---|
| **PASS** | include the rescue filtered VCF as-is |
| **WARN** | include only the `label_qc --apply` **cleaned** VCF (HIGH-confidence flags dropped, MID relabelled to Germline) |
| **FAIL** | **exclude from training entirely** — this is the abnormal-sample class (old-cohort calibration: 4081, 4255 = self-contradiction/S2 + normal contamination; 4032 = RNA-only common-AF germline leakage/S3) |

The rerun does not rehabilitate FAIL samples whose signature is intrinsic to their caller VCFs — the smoke comparison above shows 4032's signature surviving the fixed pipeline essentially unchanged. Expect (and accept) FAILs; feeding them to training is the failure mode this gate exists to prevent (old-cohort evidence: dropping 4081+4255 alone lifted F1 0.633 → 0.898).

### 7. Known pending items

- **S7 recalibration** (see 4278 smoke section): the low-DP WARN threshold was calibrated on pre-fix DP fields; expect cohort-wide S7 WARNs until recalibrated against the fixed DP extraction.
- Both smoke samples are done; the 64-sample remainder is running 6-way parallel (§1, launched 2026-08-29, `config/rerun_cohort.yaml` without VEP, ~18–20h expected). Remaining work after it completes: the cohort-mode label_qc gate (§5).
- The cohort outputs contain no VEP annotations (deliberate, §1). If VEP-annotated MAFs are needed later for specific samples, re-run those samples with `config/rerun.yaml` — completed samples are auto-skipped, so clear their state entry or use a fresh state file.
