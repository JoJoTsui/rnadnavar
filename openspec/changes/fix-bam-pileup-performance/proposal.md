## Why

BAM pileup (variant-wise depth/strand/quality validation) and BAM stats (whole-genome coverage/insert/mapping metrics) share the same input data but are implemented independently with no shared infrastructure. Three interconnected problems exist:

1. **Performance**: Pileup processes 3 BAM types sequentially with per-window queries. BAM stats does full-genome scans. Both could benefit from shared BED-guided region access and parallel execution.

2. **Coverage accuracy**: BAM stats uses a post-hoc coverage recalculation (`total_bases / bed_total`) that includes off-target reads in the numerator, overestimating WES coverage by 20-40%. The BED regions exist to define targets but aren't used to filter bases.

3. **No shared infrastructure**: The BED file is read for stats (sum total only), then would be re-read for pileup (need merged regions). Regions are not shared between subsystems.

## What Changes

### Part 1 — Shared BED Processing
- `read_and_merge_bed(path, gap=100_000)` — read BED, merge adjacent intervals, return both `bed_total` and merged `[(chrom, start, end)]`
- Merged regions feed both BAM stats (on-target coverage) and BAM pileup (region-guided queries)
- BED read once, used everywhere

### Part 2 — BAM Stats Coverage Accuracy Fix
- Pass BED regions to Rust `whole_genome_stats()`
- Track `on_target_bases` separately from `total_bases`
- `mean_coverage = on_target_bases / bed_total` (accurate for WES)
- Maintain full-genome scan for other metrics (reads, mapping, insert, mapq)

### Part 3 — Rust Combined Multi-BAM Pileup
- New `pileup_variants_multi(bam_paths, positions, bed_regions?)` function
- Opens all BAMs once, shares position grouping and window iteration across DN/DT/RT
- BED-guided regions replace 1Mb windows (~300 queries vs ~2,765)
- Binary search inner loop replaces HashMap iteration
- Returns all 3 result sets in one FFI call

### Part 4 — Parallelize BAM Stats and Variant Processing
- BAM stats (manifest-only, no variant data) runs concurrently with variant processing
- ThreadPoolExecutor launches BAM stats in background while process_single_sample proceeds
- Hides BAM stats latency (~9 min) behind variant processing (~2 min)

### Part 5 — Tests
- 5 new test methods: BED merge, on-target coverage, combined multi-BAM parity, binary search parity, BED-guided region count

## Impact

- **Rust**: `bam.rs` (+BED regions param, on-target tracking), `pileup.rs` (+multi-BAM function, BED regions, binary search)
- **Python**: `bam_stats.py` (+read_and_merge_bed), `cli.py` (parallel BAM stats, pass BED regions to pileup), `rust_bam.py` (multi-BAM wrapper)
- **Tests**: 5 new test methods
- **Performance**: 22.5 → ~75 samples/hour (4 workers, WES with BED)
