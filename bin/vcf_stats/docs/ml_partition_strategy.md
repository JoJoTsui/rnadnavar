# ML Partition Strategy

## Overview

The pipeline assigns each variant to a train/val/test partition based on
chromosome, providing a deterministic and reproducible data split for
downstream machine learning experiments.

## Partition Rules

| Partition | Chromosomes      | Rationale                              |
|-----------|-----------------|----------------------------------------|
| test      | chr1            | Largest autosome, ~8% of variants      |
| val       | chr21, chr22    | Smallest autosomes, ~2% of variants    |
| train     | All others      | Remaining ~90% of variants             |

## Design Rationale

- **Chromosome-based splitting** avoids data leakage from nearby variants
  on the same chromosome appearing in both train and test sets. Positional
  linkage disequilibrium means variants within the same genomic region are
  not independent.
- **Deterministic** — the same variant always maps to the same partition
  regardless of run order or sample composition.
- **No randomness** — unlike random splits, this is fully reproducible
  across runs without saving split indices.

## Output Files

- `stats/threshold/partition_summary.tsv` — Variant and sample counts
  per partition.
- `stats/threshold/disease_partition_summary.tsv` — Disease-by-partition
  variant counts (checks for disease imbalance across partitions).
- `stats/threshold/filter_vaf_dp_cross_tab.tsv` — FILTER x VAF x DP
  cross-tabulation, optionally partitioned.

## Usage in ML Workflows

```python
import polars as pl

df = pl.read_parquet("variant_details/*.parquet")

# The partition column is already present
train = df.filter(pl.col("partition") == "train")
val   = df.filter(pl.col("partition") == "val")
test  = df.filter(pl.col("partition") == "test")
```

## Implementation

The `partition` column is added as a lazy expression in `cli.py` after
the combined parquet scan, before any aggregation. This means it is
available to all downstream statistics and visualizations.
