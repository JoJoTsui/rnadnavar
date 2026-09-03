# Split manifest is built in seq2neo from truth-label VCFs, not by neo_var's Lance pipeline

The sample/variant split for the seq2neo cohort is produced by a repo-local
script (`examples/seq2neo/scripts/build_split_manifest.py`) that streams the
63 truth-label VCFs from `sample_manifest_rerun.tsv` and emits
`sample_split.tsv` + `selected_variants.parquet` + a provenance sidecar.
We deliberately did not port neo_var's `prepare-datasets` staging: its
`reserve_downstream.py` / `split_dataset.py` are ~90% Lance tensor I/O, and
the only portable core was the chromosome map, the reserved-ID list, and the
tag predicates. Both PASS and WARN labels use the chromosome map; verdict is an
evaluation stratum, not a reason to override chromosome assignment.

## Considered Options

- **Drive neo_var's pipeline on a Lance export** — rejected: would require a
  tensor-extraction step just to produce split assignments, and neo_var's
  code has no label-verdict concept (PASS-primary evaluation with WARN sensitivity).
- **TSV for the variant-level manifest** — rejected: ~10M rows; parquet
  (polars, zstd) is the performance format. The sample-level table stays TSV
  for diffability.
- **Rename `FILTER` to `class_label`** (as neo_var's actual writer did) —
  rejected: both strategy docs say `FILTER`, and FILTER is the frozen
  label contract with the model repos; `class_label` was unilateral drift.

## Consequences

- **PASS-primary evaluation with WARN sensitivity**: both verdicts use the chromosome
  map, so WARN chr1/21/22 records go to test/val rather than train. PASS metrics
  remain primary; WARN metrics are reported separately. The reserved pool is PASS-only.
- **Tags are recomputed from label-VCF INFO fields**, not from neo_var's
  `variant_details/*.parquet`. `DP_DNA_MEAN` (caller-reported mean depth)
  stands in for neo_var's `BAM_DT_DP` (tumor-BAM pileup depth) in the
  `is_low_dp` predicate. This is not a like-for-like swap: BAM pileup depth
  exists at every locus, while `DP_DNA_MEAN` exists only where a DNA caller
  emitted a record — so `is_low_dp` can only fire on DNA-detected variants,
  and the pool shrinks from thousands (neo_var) to 17 Somatic cohort-wide.
  Missing values follow neo_var null rules (missing `DP_DNA_MEAN` ⇒
  `is_low_dp` False; missing RNA fields ⇒ RNA non-zero False).
- **FILTER whitelist (`Somatic/Germline/Reference`) is applied at
  manifest-build time**, absorbing a step that used to live in tensor
  extraction; downstream validation consumes the manifest directly.
- **Reserved sub-pool thinness is accepted, not fixed**: per-sample Somatic
  counts for `low_vaf_a` / `low_vaf_b` / `low_dp` / `rescued` are
  single/double digits for some reserved samples (4231 has zero rescued and
  zero low_vaf_b; pooled `low_dp` is 17 Somatic cohort-wide). These
  sub-pools are evaluated pooled across the 5 reserved samples only;
  the reserved set is chosen for disease-fold diversity, and swapping
  members for tag richness would break fold coverage.
- Row order is deterministic and `--check` regenerates-and-diffs, so the
  manifest is a reproducible build artifact, not a cached result.
