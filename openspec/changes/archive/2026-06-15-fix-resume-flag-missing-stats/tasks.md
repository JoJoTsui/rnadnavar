## 1. Fix --resume block in cli.py

- [x] 1.1 Reload `all_stats` from `sample_summary.tsv` when the file exists (line 566-573, with `.csv` fallback)
- [x] 1.2 Reload `bam_stats_df` from `bam_stats.tsv` when the file exists (line 576-583, with `.csv` fallback)
- [x] 1.3 Filter variant count to respect `--max-samples`/`--set`/`--sample-ids` (line 544-552)

## 2. Verification

- [x] 2.1 Run full pipeline (12 samples), then re-run with `--resume` — verify all CSVs and charts regenerate
- [x] 2.2 Run with `--resume --max-samples 6` — verify total_variants matches filtered count
- [x] 2.3 Run with `--resume` when CSVs don't exist — verify graceful degradation (no crash)
