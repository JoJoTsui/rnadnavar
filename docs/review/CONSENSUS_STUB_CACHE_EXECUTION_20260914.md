# Synthetic consensus/cache execution evidence (2026-09-14)

This bounded check closes the executable portion of the validation work without
touching SEQC2/HG008 outputs, production inputs, or shared Nextflow caches.
The fixture is three distinct copies of a small caller VCF (Mutect2, Strelka,
DeepSomatic) under `/tmp/rnadnavar-fake-vcfs/`.

## End-to-end stub

Command shape:

```text
nextflow run . -profile test -stub-run \
  -work-dir /tmp/rnadnavar-consensus-stub-clean5 \
  --step consensus --tools consensus,rescue,filtering \
  --skip_tools vcftools \
  --consensus_expected_callers deepsomatic,mutect2,strelka \
  --input /tmp/rnadnavar-fake-vcfs/callers_distinct.csv
```

The run completed successfully through `VCF_NORMALIZE`, `VCF_CONSENSUS`,
`VCF_FILTER`, and `MULTIQC`. It produced valid compressed/indexed artifacts,
including:

- `consensus/FAKE_T_vs_FAKE_N/FAKE_T_vs_FAKE_N.consensus.vcf.gz[.tbi]`
- `filtered/FAKE_T_vs_FAKE_N/FAKE_T_vs_FAKE_N.filtered.vcf.gz[.tbi]`
- normalized caller VCFs and indexes for all three callers

The stub modules now preserve a valid input VCF and create a real tabix index;
empty compressed placeholders are rejected by `cyvcf2` and were the reason the
original synthetic graph stopped at consensus.

## Resume/cache check

The identical command was rerun with `-resume` and the same work directory.
The log showed cache hits for all three `VT_DECOMPOSE` tasks, all three
`BCFTOOLS_NORM` tasks, all three `BCFTOOLS_STATS` tasks, and the reference
`SAMTOOLS_FAIDX` task. Consensus/filter/MultiQC were submitted under the new
run identity, so this check does **not** claim downstream label-task cache hits.
It does prove that a policy-only consensus-stage invocation can reuse the
upstream normalization/reference/statistics cache boundary in an isolated
work directory.

No mapping, variant caller, realignment, source-input, production-output, or
shared-cache path was used or modified.
