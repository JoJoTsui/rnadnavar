## Context

BAM pileup validates variant calls against raw alignment data. For each variant position, it computes total depth, REF vs ALT allele support, strand bias (F1R2/F2R1), base quality, and mapping quality from the alignment BAM files.

The current implementation uses noodles-bam in Rust (with pysam fallback), correctly computing all metrics, but with a per-position BAI query loop that creates unacceptable performance and a 2-column join that causes data misattribution at multiallelic sites.

## Design

### Part 1 — Windowed BAM Queries (pileup.rs)

**Current algorithm:**
```
for each of N positions:
    region = Region(chrom, pos..=pos)       // single-base region
    query  = reader.query(header, index, region)  // BAI lookup + BAM seek
    for each read in query:
        classify base (REF/ALT), count DP, strand, BQ, MQ
```

**New algorithm:**
```
1. Group positions by (chromosome, window)
   — Sort positions by (chrom, pos)
   — Divide into windows of ~1,000,000 bp per chromosome
   — Build HashMap<(chrom, window_start), Vec<(pos, orig_idx, ref_byte, alt_byte)>>

2. For each window:
   region = Region(chrom, window_start..=window_end)  // 1Mb range
   query  = reader.query(header, index, region)        // ONE query per window
   
   // Build position lookup within this window
   pos_map: HashMap<i64, Vec<(usize, u8, u8)>>  // pos → [(orig_idx, ref_byte, alt_byte)]

   for each read in query:
       if read is unmapped or duplicate: skip
       for each position that this read overlaps:
           if position is in pos_map:
               for each (orig_idx, ref_byte, alt_byte) for this position:
                   if base == ref_byte: ref_dp += 1
                   elif base == alt_byte: alt_dp += 1
                   update strand counts, BQ, MQ

3. Return results in original input order (via orig_idx)
```

**Performance estimate:**
- 24 chromosomes × ~3 windows each = ~72 windows per BAM type
- ~72 queries × 3 BAM types = ~216 queries total
- Each query: BAI lookup + sequential BAM read within window
- Estimated: 30-60 seconds per BAM type, ~2-3 minutes total

**Key implementation detail:** A single read can overlap multiple target positions. For a 150bp read, it could overlap up to 150 positions. The algorithm must check each position the read covers against `pos_map`. This is O(read_length × positions_per_window) per read, but read_length ≤ 150 and positions_per_window is bounded.

Alternative approach if per-read position check is slow: use a sliding window over positions. Since positions are sorted and reads are sorted by alignment start, we can advance through both lists in parallel.

### Part 2 — 4-Column Join Fix

**rust_bam.py change:**
```python
# Before:
result["CHROM"] = chroms
result["POS"] = poss
return pl.DataFrame(result)

# After:
result["CHROM"] = chroms
result["POS"] = poss
result["REF"] = refs   # ADDED
result["ALT"] = alts   # ADDED
return pl.DataFrame(result)
```

**cli.py change:**
```python
# Before:
df = df.join(
    pileup_df.select(["CHROM", "POS"] + list(rename.values())),
    on=["CHROM", "POS"], how="left",
)

# After:
df = df.join(
    pileup_df.select(["CHROM", "POS", "REF", "ALT"] + list(rename.values())),
    on=["CHROM", "POS", "REF", "ALT"], how="left",
)
```

This ensures pileup data for (chr1, 100, A, G) joins only to the row with (chr1, 100, A, G), not to (chr1, 100, A, T) at multiallelic sites.

### Part 3 — Wire --pileup-mode Flag

```python
# In process_single_sample(), after building positions:
if not no_pileup:
    # Filter positions by pileup mode
    if pileup_mode == "filtered" and "FILTER" in df.columns:
        mask = df["FILTER"] != "NoConsensus"
        positions = [
            (row[0], row[1], row[2], row[3])
            for row, keep in zip(
                df.select(["CHROM", "POS", "REF", "ALT"]).iter_rows(),
                mask.to_list()
            ) if keep
        ]
    else:
        positions = [(row[0], row[1], row[2], row[3])
                     for row in df.select(["CHROM", "POS", "REF", "ALT"]).iter_rows()]
    
    n_pos = len(positions)
    print(f"  [{sample_id}] BAM pileup: {n_pos} positions ({pileup_mode} mode)")

    for bt in ["DN", "DT", "RT"]:
        ...
        print(f"  [{sample_id}] BAM pileup {bt}: starting ({n_pos} positions)...")
        pileup_df = do_pileup(bam_path, positions)
        _mem(f"after pileup {bt}")
```

### Part 4 — Tests

**New tests in `TestRustPileup`:**

1. `test_pileup_output_includes_ref_alt` — Result dict has "REF" and "ALT" keys
2. `test_pileup_multiallelic_join_correct` — 4-column DataFrame join matches correct alleles
3. `test_pileup_windowed_matches_per_position` — Windowed query produces identical results to per-position for 100 positions
4. `test_pileup_10k_positions_completes_quickly` — 10,000 positions complete in < 60 seconds

**New test class `TestPileupIntegration`:**

5. `test_pileup_mode_filtered_excludes_noconsensus` — Filtered mode positions < all mode positions
6. `test_process_single_sample_writes_pileup_columns` — Parquet output has BAM_DN_DP etc.
7. `test_pileup_join_uses_4_columns` — Join key is (CHROM, POS, REF, ALT)

## Risks

- **Windowed query memory**: A 1Mb window at 100× coverage contains ~30K reads. With decompressed records, memory is < 100 MB per window — acceptable.
- **Rust rewrite complexity**: The windowed algorithm is more complex than per-position queries. Fall back to per-position if windowed mode fails.
- **BAI index requirement**: Windowed queries still need BAI index. If missing, fall back to per-position with pysam.
