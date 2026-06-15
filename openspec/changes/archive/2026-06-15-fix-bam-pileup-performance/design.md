## Context

BAM pileup validates variant calls against alignment data. BAM stats computes whole-genome quality metrics. Both operate on the same BAM files but share no infrastructure. The `--bed` flag defines WES target regions but is only used for coverage denominator post-hoc recalculation — the actual region boundaries are discarded. This design unifies BED processing and implements a Rust combined multi-BAM pileup function.

## Goals / Non-Goals

**Goals:** Shared BED processing, accurate on-target coverage in BAM stats, Rust combined multi-BAM pileup, binary search inner loop, parallel BAM stats + variant processing.

**Non-Goals:** Changing BAM stats metrics beyond coverage, modifying the Nextflow pipeline, changing rescue VCF semantics.

---

## Part 1 — Shared BED Processing

### Problem

`sum_bed_regions()` computes total BED length but discards the merged region coordinates. Pileup would need to re-read and re-merge the same BED file.

### Design

```python
def read_and_merge_bed(bed_path: str, gap: int = 100_000) -> tuple[int, list[tuple[str, int, int]]]:
    """Read BED file, merge adjacent intervals, return total + merged regions.

    Merging within 'gap' bp reduces ~50K raw exons to ~300 contiguous regions.
    These regions feed both BAM stats (on-target coverage) and BAM pileup
    (region-guided queries).
    """
    intervals = []
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3: continue
            try:
                intervals.append((parts[0], int(parts[1]), int(parts[2])))
            except ValueError: continue

    intervals.sort(key=lambda x: (x[0], x[1]))
    merged = []
    for chrom, start, end in intervals:
        if (merged and merged[-1][0] == chrom
            and start - merged[-1][2] <= gap):
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))

    bed_total = sum(end - start for _, start, end in merged)
    return bed_total, merged
```

Calling code in cli.py:
```python
bed_total = 0
bed_regions = None
if args.bed:
    bed_total, bed_regions = read_and_merge_bed(args.bed)
    print(f"BED: {bed_total:,} bp, {len(bed_regions)} merged regions")
```

---

## Part 2 — BAM Stats On-Target Coverage

### Problem

Current Rust `whole_genome_stats()` uses `total_query_length` (all mapped reads) as coverage numerator. Post-hoc recalculation multiplies by `bam_ref/bed_total` but numerator still includes off-target reads.

### Design

Add `bed_regions` parameter to `whole_genome_stats()`:

```rust
pub fn whole_genome_stats(
    bam_path: &Path,
    max_reads: u64,
    bed_regions: Option<&Vec<(String, i64, i64)>>,
) -> Result<BamStats, Box<dyn std::error::Error>> {
    let mut on_target_bases: u64 = 0;

    // Build per-chromosome BED interval list for efficient lookup
    let bed_index = build_bed_index(bed_regions);

    for record in reader.records() {
        total += 1;
        if !flags.is_unmapped() {
            mapped += 1;
            total_query_length += seq_len;

            // On-target check: only for coverage denominator
            if let Some(ref bed) = bed_index {
                if read_overlaps_bed(&record, bed) {
                    on_target_bases += seq_len;
                }
            }
            // ... other metrics unchanged
        }
    }

    let denominator = bed_regions.as_ref()
        .map(|regions| sum_region_lengths(regions))
        .unwrap_or(ref_lengths);

    let coverage_bases = if bed_regions.is_some() { on_target_bases } else { total_query_length };
    mean_coverage = coverage_bases as f64 / denominator as f64;
}
```

Performance: the on-target check uses a two-pointer walk since BAM records and BED regions are both sorted by coordinate. O(reads + regions) per chromosome.

---

## Part 3 — Rust Combined Multi-BAM Pileup

### Function Signature

```rust
pub fn pileup_variants_multi(
    bam_paths: &[String],
    bam_labels: &[String],
    chroms: &[String],
    positions: &[i64],
    ref_bases: &[String],
    alt_bases: &[String],
    bed_regions: Option<&[(String, i64, i64)]>,
) -> Result<Vec<(String, Vec<PileupResult>)>, String>
```

### Algorithm

```
Phase 1: Region grouping (shared across all BAMs)
  1. Determine query regions:
     - BED provided → merged BED intervals
     - No BED → 1Mb sliding windows (current)
  2. Group positions into regions:
     region_map: HashMap<(chrom, start), Vec<(pos, orig_idx, ref_byte, alt_byte)>>

Phase 2: Per-region processing
  3. Open all BAMs + BAI indices
  4. Build chrom_lengths from first BAM header
  5. FOR each region:
     a. Build sorted position array for binary search
     b. FOR each BAM:
        - Clamp region to chromosome length
        - Create Region(chrom, start..=end)
        - reader.query(&header, &index, &region)
        - FOR each read:
          * Skip unmapped/duplicate
          * Binary search positions overlapping [align_start, align_end)
          * FOR each overlapping position:
            → Update PileupResult accumulators
     c. Reads from all BAMs dropped after region

Phase 3: Finalize
  6. Divide BQ/MQ sums by counts → means
  7. Return Vec<(bam_label, Vec<PileupResult>)>
```

### Binary Search Inner Loop

```rust
// Per window: build sorted position array (once)
let sorted_pos: Vec<i64> = positions.iter().copied().sorted().collect();

// Per read: binary search for overlap range
for each read:
    let align_end = align_start + seq_len as i64;
    let start_idx = sorted_pos.binary_search_by(|p| p.cmp(&align_start))
        .unwrap_or_else(|i| i);

    for pos in &sorted_pos[start_idx..] {
        if *pos >= align_end { break; }  // past read end, stop
        // Process this position
        for (orig_idx, ref_byte, alt_byte) in pos_map[pos] {
            // update PileupResult
        }
    }
```

Complexity: O(n_reads × log(n_positions) + n_overlaps) per window.

### BED Region Mode

When `bed_regions` is provided:
- Regions are the merged BED intervals (not 1Mb windows)
- ~300 queries for WES (vs ~2,765 for genome-wide windows)
- Skip regions with no target positions
- Clamp regions to chromosome length (existing safety)

---

## Part 4 — Parallelize BAM Stats and Variant Processing

BAM stats (manifest-only) is independent of variant processing (needs positions from VCF). They can run concurrently:

```python
# In main(), launch BAM stats in background before variant processing
if not args.no_bam:
    with ThreadPoolExecutor(max_workers=1) as bam_bg:
        bam_future = bam_bg.submit(
            compute_all_bam_stats, rows, max_workers=args.bam_workers,
            bed_total=bed_total
        )

# ... variant processing runs in parallel ...

# Collect BAM stats result
if not args.no_bam:
    bam_stats_df = bam_future.result()
```

This hides BAM stats latency (~9 min for large samples) behind variant processing (~2-3 min), making total wall-clock time max(variant, BAM stats) instead of sum.

---

## Part 5 — Tests

### New tests in TestRustPileup:
- `test_pileup_multi_bam_parity` — combined function matches individual calls
- `test_bed_merge_intervals` — merge within 100Kb gap produces correct regions
- `test_pileup_bed_guided_region_count` — BED mode creates fewer regions than 1Mb mode
- `test_pileup_binary_search_matches_hashmap` — binary search produces identical results

### New tests in TestBamStats:
- `test_bam_stats_on_target_coverage` — BED-filtered coverage lower than total coverage
- `test_read_and_merge_bed_total` — merged total matches input intervals

## Risks

- **Multi-BAM memory**: 3 result sets = ~320 MB peak. Acceptable for 200 GB machines.
- **BED interval tree**: Per-chromosome interval lookup needs correct two-pointer implementation. Fall back to full scan if intervals exceed 10K.
- **Binary search correctness**: Must handle edge cases (position at exact alignment start/end). Validated by parity test against HashMap approach.
