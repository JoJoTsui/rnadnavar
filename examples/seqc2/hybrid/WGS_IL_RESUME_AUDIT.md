# WGS-IL resume and downstream argument audit

Checked 2026-09-12 against the failed `mad_lamarck` run, Nextflow 25.10.0,
session `c0128e76-65d4-46cd-be43-0eda266c5d1f`.

## Fix and cache scope

`conf/modules/prepare_realignment/vcf_realignment.config` now overrides
mosdepth arguments only inside `RNA_REALIGNMENT_WORKFLOW`. When the module
receives a BED, WGS extra arguments are `-n --fast-mode`; the module supplies
the single `--by <BED>`. First-pass WGS still uses `-n --fast-mode --by 500`.
WES retains its existing empty extra arguments.

The original trace contains 126 completed tasks. All 126 work directories,
`.command.sh` files, and exit-code-0 markers are present. The original session
cache is present. None of these completed tasks matches the changed selector:
realignment mosdepth failed during command construction before submission.
This fix leaves completed task scripts, inputs, resources, and process names
unchanged. Keep the original work directory and session cache; do not clean
them or change input files before resuming. Actual cache hits are confirmed
only when Nextflow resumes; this audit did not launch production tasks.

Regression check (four WGS/WES × first-pass/realignment cases):

```bash
java -cp /t9k/mnt/joey/micromamba/envs/nextflow/share/nextflow/dist/25.10.0/nextflow-25.10.0-one.jar \
  groovy.ui.GroovyMain tests/realignment_mosdepth_config.groovy
```

The test reproduced the exact duplicate-`--by` exception before the fix and
passes afterward. It uses Nextflow's config parser and selector resolution,
then evaluates the actual mosdepth module script without executing it.

## Downstream audit

The effective root + WGS-IL configuration was resolved with Nextflow's
`ConfigBuilder`. All 29 registered process scripts under RNA realignment,
second rescue, and MultiQC rendered and passed `bash -n`, using representative
input bindings. This includes registered conditional merge processes that
receive no tasks in the single-interval run. Python CLI flags were additionally
checked against the installed parsers' `--help` output. No additional argument
errors were found for this run's configuration.

| Stage | Checked arguments and routing |
| --- | --- |
| Preprocessing | MarkDuplicates settings unchanged; SplitNCigarReads receives the merged realignment BED and `--create-output-bam-index false`; CRAM merge/index settings resolve. Second-pass BQSR is skipped by the workflow. |
| Strelka | Original DNA normal + realigned RNA tumor; `--callRegions` receives indexed realignment intervals; `--exome` is enabled for RNA; 46 CPUs. |
| Mutect2 | Same matched normal, realignment intervals, germline resource and panel of normals; `--normal-sample WGS_IL_N_1`, soft-clipped bases disabled, F1R2 output specified. |
| FilterMutectCalls | Realignment selector supplies `--max-events-in-region 5`; orientation, segmentation, contamination and estimate inputs are empty, so no nonexistent table arguments are emitted. |
| DeepSomatic | Correct normal/tumor ordering; realignment `--regions`; configured WGS model; 4 shards, 64 GB. Model choice is existing policy, not a validation of model accuracy on RNA. |
| VCF QC | bcftools stats and vcftools Ts/Tv-count, Ts/Tv-quality and FILTER-summary options render. |
| Normalization | VT decomposition; bcftools `--output-type z --multiallelics -any --check-ref w`; both TBI and CSI indexing. |
| RNA consensus | Expected panel `deepsomatic,mutect2,strelka`; SNV/indel thresholds 2; minimum alt support 3; native-evidence SNV policy retained. |
| RNA filtering | `--gnomad_thr 0.0001 --min_alt_reads 2`, multiallelic filtering and reference input; whitelist/blacklist are unset. |
| Second rescue | First-pass DNA + realigned RNA consensus and all three individual caller VCFs per modality; `--alignment-round realignment`; promotion enabled; DNA/RNA minimums 1/2; DNA Artifact veto. |
| Annotation | COSMIC/gnomAD thresholds 0.001/2/5 and four workers; REDIportal minimum RNA support 2. Required database paths exist. |
| Rescue filtering | Filtered and stripped outputs; `--filter_multiallelic`; CLI options accepted. |
| VEP | Runs for second rescue even without `vep` in `--tools`; offline, merged cache, GRCh38, human, release 115, VCF/bgzip/tabix, 24 forks. Cache `homo_sapiens_merged/115_GRCh38/info.txt` exists; Conda specification requests VEP 115.2. |
| MultiQC | Standard config/report collection and empty optional title arguments render. |

Resources resolve within the configured per-task ceiling of 46 CPUs / 72 GB.
Rendering is not execution: tool runtime behavior, data-dependent failures,
new Conda environment creation, and scientific benchmark performance remain
to be established by the resumed run.

The post-Nextflow benchmark wrapper was also checked. It requires the final
realignment-rescue VEP VCF and runs `som.py` for consensus, three DNA callers,
and rescue, followed by metrics aggregation. Installed `som.py` accepts its
arguments; truth SNV/indel VCFs, HC BED, reference, bcftools and tabix exist.
The actual comparison domain is **high-confidence regions (`-R`) intersected
with the configured UKB/broad union target BED (`-T`)**. WGS describes the
input assay; these scores are not unrestricted whole-genome scores. `-N`
normalizes both truth and query; non-PASS caller records remain excluded.

## Resume

From the repository root, resume the exact failed session and preserve its
Conda cache location:

```bash
NXF_CONDA_CACHEDIR=/t9k/mnt/joey/nf_conda_envs NXF_CONDA_USEMAMBA=true \
micromamba run -n nextflow nextflow run main.nf \
  -c examples/seqc2/hybrid/seqc2.hybrid.wgs.config \
  --input examples/seqc2/hybrid/csv/seqc2_wgs_il_hybrid.csv \
  --outdir examples/seqc2/hybrid/output/seqc2.wgs.il.hybrid \
  -work-dir /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/nf_work \
  -offline -with-conda -resume c0128e76-65d4-46cd-be43-0eda266c5d1f
```

Alternatively, the existing shared benchmark wrapper resumes Nextflow using
`config_wgs_il.yaml` (`resume: true`, original work/Conda directories), then
runs comparisons automatically. Pass a new fourth-argument log path to retain
the original failure log. Its bare `-resume` selects the most recent session,
so use the explicit command above if another pipeline run has since intervened.
