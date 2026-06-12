## 1. Fix --resume block in cli.py

- [ ] 1.1 Reload `all_stats` from `sample_summary.csv` when the file exists
- [ ] 1.2 Reload `bam_stats_df` from `bam_stats.csv` when the file exists
- [ ] 1.3 Filter variant count to respect `--max-samples`/`--set`/`--sample-ids`

## 2. Verification

- [ ] 2.1 Run full pipeline (12 samples), then re-run with `--resume` — verify all CSVs and charts regenerate
- [ ] 2.2 Run with `--resume --max-samples 6` — verify total_variants matches filtered count
- [ ] 2.3 Run with `--resume` when CSVs don't exist — verify graceful degradation (no crash)
