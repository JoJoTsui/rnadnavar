# examples/seqc2 — SEQC2 DNA-only benchmark

SEQC-II (SEQC2) somatic-mutation reference run: WES tumor-normal pair
`WES_LL_T_1` (tumor, status=1) vs `WES_LL_N_1` (normal, status=0), GRCh38,
DNA-only (`rna=false`), callers `deepsomatic,mutect2,strelka` + `norm,consensus`.

Modeled on `examples/seq2neo/` (config copied and edited).

## Layout

- `seqc2.shared.config` — Nextflow `-c` config. `input`/`outdir` are passed per run.
- `config/runner.yaml` — execution config for the driver (env, paths, completion globs).
  Defaults auto-resolve from the repo location; no edits needed on an rsynced copy.
- `csv/seqc2_wes_ll.csv` — BAM samplesheet (`patient,status,sample,lane,bam,bai`).
- `run_wes_ll.sh` — thin wrapper over the driver (same pattern as seq2neo's
  `run_single.sh`).
- `scripts/run_pipeline.py` — driver: builds and runs the micromamba/nextflow
  command, tracks run state in `runs/run_state.json`, checks completion artifacts.
- `scripts/run_benchmark.sh` + `aggregate_benchmark.py` — benchmark tooling
  (compares caller/consensus VCFs against the SEQC2 truth).

## Input modes

- **BAM (default)**: config sets `step = "variant_calling"`; samplesheet carries
  `bam,bai` columns. Input BAMs are expected to be aligned + deduplicated
  (`*.bwa.dedup.bam`); BQSR is skipped in this entry mode.
- **FASTQ**: use `fastq_1,fastq_2` columns and override with `--step mapping`
  (alignment via `bwa-mem`, configured in the shared config).

## Data

- BAMs: `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/data/giab/data/seqc2/WES/`
- Capture intervals: `/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/intervals/ukb.pad50.broad.pad50.union.bed`
- Truth set (v1.2.1) + High-Confidence BED: `.../seqc2/truth/`

## Run

```bash
bash run_wes_ll.sh                      # run (resumes automatically)
bash run_wes_ll.sh --dry-run            # print the nextflow command only
bash run_wes_ll.sh --outdir /path/to/outdir
```

The driver skips the run if the completion artifacts in `config/runner.yaml`
(execution trace + consensus VCF) are already present in the outdir.

## Benchmark

After the pipeline completes, compare consensus + per-caller VCFs against the
SEQC2 truth set (som.py, metrics split by SNV/indel):

```bash
bash scripts/run_benchmark.sh                      # defaults: output/seqc2.wes.ll, WES_LL pair
bash scripts/run_benchmark.sh /path/to/outdir PAIR /path/to/comparison_dir
```

Produces `comparison/<PAIR>/<PAIR>.benchmark_comparison.csv` with
TP/FP/FN/Precision/Recall/F1 per query (consensus/mutect2/deepsomatic/strelka)
× variant type (snv/indel/all). Requires the `happy` micromamba env (som.py)
and bcftools/tabix on PATH.
