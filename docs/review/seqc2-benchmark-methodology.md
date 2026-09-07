# SEQC2 benchmark methodology audit

Review date: 2026-09-07. Code baseline: `931c476a081162e6398130a5fd28b1fd8aced11c`.
This appendix separates scoring limitations from demonstrated classification failures.

## Training endpoint clarification

The user clarified that downstream training uses rescue labels derived from realigned RNA, including subsequent re-consensus/QC derivatives, because they have observed fewer false positives there. The local manifest audit supports the realigned-RNA lineage. Accordingly, the primary endpoint is the final training-label artifact versus DeepSomatic DNA; the completed first-round SEQC2 rescue comparison is an intermediate diagnostic. Evaluate realignment's effect by comparing matched annotation/selection stages, with gained/lost TP and FP. A failed or missing second pass cannot be represented by the first-round score.

## What the current scores measure

`examples/seqc2/scripts/run_benchmark.sh:31–37,70–85` compares SEQC2 v1.2.1 merged SNV/indel truth using `som.py -N -R <High-Confidence_Regions_v1.2.bed> -T <ukb.pad50.broad.pad50.union.bed> -r <GRCh38 FASTA>`. It selects consensus/rescue `FILTER=Somatic` into benchmark-only PASS copies; individual callers use PASS selection. Production biological FILTER vocabulary is preserved.

The installed evaluator is `hap.py 0.3.15`, build `py27hcb73b3d_0`, from `/t9k/mnt/joey/micromamba/envs/happy/conda-meta/hap.py-0.3.15-py27hcb73b3d_0.json`. Its `bin/som.py:452–477` applies normalization, regions, targets, and PASS preprocessing to truth and query; line 495 intersects them with `bcftools isec`. Consequently, these are normalized position/allele comparisons, not a haplotype-aware `hap.py`/vcfeval comparison. Complex indels and MNV representations warrant an additional representation-aware sensitivity analysis. This does not negate the observed SNV false-positive burden. [Upstream som.py documentation](https://github.com/Illumina/hap.py/blob/master/doc/sompy.md).

The main tables share a denominator of 2,300 truth records: 2,205 SNVs and 95 indels. Even the WGS-input run uses this exome-target evaluation domain. A WGS input does not make its score genome-wide. WES scores against a padded union BED also conflate caller sensitivity with the sample's actual capture and coverage. Keep the historical domain for regression, and report sample capture and prespecified coverage strata alongside it. Do not choose callable regions from variants that happened to pass a caller.

SEQC2's truth calls and high-confidence regions are intended to be used together. The current high-confidence restriction is appropriate; unexpressed DNA truth sites remain in the denominator for an RT-only comparison unless a separately defined RNA-callable stratum is used. A DNA-somatic truth benchmark does not validate RNA-editing discovery. [SEQC2 consortium protocol](https://sites.google.com/view/seqc2/home/data-analysis/high-confidence-somatic-snv-and-indel-v1-2).

## DeepSomatic training overlap

DeepSomatic's **r1.9** WGS and WES HCC1395 examples explicitly evaluate chromosome 1 held out during training. The repository's whole-target HCC1395 comparison is therefore useful as a pipeline regression benchmark, but it does not establish performance on a held-out cell line or on entirely held-out genomic loci. This is a limitation of the evidence, not proof that the model memorized any particular error. It also does not explain away losses introduced downstream on the same DeepSomatic calls.

Report chr1 separately without tuning against it; use a genuinely independent sample or verified model-specific held-out set for broader superiority claims. Exact checkpoint training provenance should be attached to that evaluation. [DeepSomatic r1.9 WGS example](https://raw.githubusercontent.com/google/deepsomatic/r1.9/docs/deepsomatic-case-study-wgs.md), [r1.9 WES example](https://raw.githubusercontent.com/google/deepsomatic/r1.9/docs/deepsomatic-case-study-wes.md).

The official model list includes DNA WGS/WES and other sequencing technologies, with no listed RNA-seq model. Applying WES to RNA tumor versus DNA normal is a distribution-shift hypothesis requiring direct measurement; it is not evidence of validated RNA performance. The existing RT run must be interpreted on its own measured performance. [DeepSomatic model interface](https://github.com/google/deepsomatic).

## Reproducibility and stage identity

The benchmark script reuses merged truth and Somatic-to-PASS query derivatives whenever an index exists (`run_benchmark.sh:52–68,94–105`). It does not invalidate these caches when the source VCF, truth, code, or policy changes. This is a demonstrated implementation risk, not a claim that the inspected hybrid result is stale: the artifact audit confirmed its rescue derivative matches the filtered first-round output.

For each next comparison, use a fresh directory or content-keyed derivatives and record source checksums, variant-content checksums, source stage, code revision, executed arguments, model/container identity, database releases, BED checksums, and evaluator versions. Capture both raw and annotated rescue, because annotation can change representation and biological labels. The current aggregator resolves SNV/indel rows by semantic labels (`aggregate_benchmark.py:63–95`), avoiding the historical positional row-swap problem; check older CSVs against their native metrics rather than trusting their filenames.

## Required reporting for proposed changes

1. Preserve existing domain metrics for regression, separately by SNV/indel and input pair. Add capture/coverage, chr1 holdout, and RNA-expression strata with explicit denominators.
2. Record paired transitions relative to DeepSomatic and to DNA consensus: gained TP, gained FP, lost TP, removed FP. Net changes alone hide harmful exchanges.
3. Report incremental rescue yield among newly admitted variants and break it down by decision rule, caller combination, DNA/normal support, alt count, RNA-editing overlap, splice proximity, and variant type.
4. Compare initial and realigned RNA under the same annotation settings and on the same evaluation domain. Candidate-only recall is conditional on having entered the realignment candidate set; retain a separate end-to-end denominator.
5. Report uncertainty for the small indel denominator. A few changed indels are material; a single pooled F1 score hides this.

These are review recommendations. No benchmark policy or workflow implementation has been changed.
