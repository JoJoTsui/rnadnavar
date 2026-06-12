## Context

`--resume` was added as a convenience flag to skip variant processing and BAM stats when re-running the pipeline. However, `all_stats` and `bam_stats_df` remain empty (initialized to `[]` and `pl.DataFrame()`), causing cross-sample CSVs and per-sample/BAM charts to be silently skipped.

## Decisions

### Reload from existing CSVs

When `--resume` is active and the CSV files exist, reload `all_stats` and `bam_stats_df`:

```python
if args.resume:
    # ... count variants (existing) ...
    args.no_bam = True

    stats_csv = output_dir / "sample_summary.csv"
    if stats_csv.exists():
        all_stats = pl.read_csv(str(stats_csv)).to_dicts()

    bam_csv = output_dir / "bam_stats.csv"
    if bam_csv.exists():
        bam_stats_df = pl.read_csv(str(bam_csv))

    # Count only filtered samples
    if args.max_samples or args.set or args.sample_ids:
        sample_ids = {r["sample_id"] for r in rows}
        total_variants = 0
        for pq in _glob.glob(os.path.join(variant_dir_str, "*_variants.parquet")):
            sid = os.path.basename(pq).replace("_variants.parquet", "")
            if sid in sample_ids:
                total_variants += pl.scan_parquet(pq).select(pl.len()).collect().item()
```

If CSVs don't exist (partial previous run), `all_stats`/`bam_stats_df` stay empty — graceful degradation.

## Risks

- CSV reload types may differ from original dict types (int↔float) — minor, acceptable for convenience flag
- Stale CSVs from different code version — user should re-run without `--resume`
