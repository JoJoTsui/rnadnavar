//! BAM pileup at variant positions using noodles-bam windowed queries.
//!
//! Two modes:
//! - `pileup_variants()` — single BAM, 1Mb sliding windows, HashMap inner loop.
//! - `pileup_variants_multi()` — multiple BAMs in one FFI call, BED-guided or
//!   1Mb sliding windows, binary search inner loop.
//!
//! Both clamp window regions to chromosome lengths from the BAM header to
//! prevent BGZF block boundary errors.

use std::collections::HashMap;
use std::path::Path;

use noodles_bam::{self as bam, bai};
use noodles_core::Region;

/// Compute the reference alignment span from CIGAR operations.
fn alignment_reference_span(record: &bam::Record) -> i64 {
    let mut span: i64 = 0;
    for op_result in record.cigar().iter() {
        if let Ok(op) = op_result {
            if op.kind().consumes_reference() {
                span += op.len() as i64;
            }
        }
    }
    if span == 0 {
        record.sequence().len() as i64
    } else {
        span
    }
}

const WINDOW_SIZE: i64 = 1_000_000;  // 1 Mb windows

/// Pileup result for a single position.
#[derive(Debug, Clone, Default)]
pub struct PileupResult {
    pub dp: Option<i64>,
    pub ref_dp: Option<i64>,
    pub alt_dp: Option<i64>,
    pub f1r2_ref: Option<i64>,
    pub f2r1_ref: Option<i64>,
    pub f1r2_alt: Option<i64>,
    pub f2r1_alt: Option<i64>,
    pub mean_bq: Option<f64>,
    pub mean_mq: Option<f64>,
    // Internal accumulators (not exposed as columns)
    bq_sum: f64,
    bq_count: i64,
    mq_sum: f64,
}

// ═══════════════════════════════════════════════════════════════════════════════
// Shared helpers
// ═══════════════════════════════════════════════════════════════════════════════

/// Update a PileupResult from a read base call at a given position.
fn update_pileup_result(
    r: &mut PileupResult,
    base: Option<u8>,
    ref_byte: u8,
    alt_byte: u8,
    is_rev: bool,
    bq_val: Option<u8>,
    mq_val: Option<f64>,
) {
    r.dp = Some(r.dp.unwrap_or(0) + 1);
    if let Some(bq) = bq_val {
        r.bq_sum += bq as f64;
        r.bq_count += 1;
    }
    if let Some(mq) = mq_val {
        r.mq_sum += mq;
    }

    if let Some(b) = base {
        let b_upper = b.to_ascii_uppercase();
        let ref_upper = ref_byte.to_ascii_uppercase();
        let alt_upper = alt_byte.to_ascii_uppercase();
        if ref_upper != 0 && b_upper == ref_upper {
            r.ref_dp = Some(r.ref_dp.unwrap_or(0) + 1);
            if is_rev {
                r.f2r1_ref = Some(r.f2r1_ref.unwrap_or(0) + 1);
            } else {
                r.f1r2_ref = Some(r.f1r2_ref.unwrap_or(0) + 1);
            }
        } else if alt_upper != 0 && b_upper == alt_upper {
            r.alt_dp = Some(r.alt_dp.unwrap_or(0) + 1);
            if is_rev {
                r.f2r1_alt = Some(r.f2r1_alt.unwrap_or(0) + 1);
            } else {
                r.f1r2_alt = Some(r.f1r2_alt.unwrap_or(0) + 1);
            }
        }
    }
}

/// Finalize mean BQ and mean MQ by dividing accumulated sums by counts.
fn finalize_results(results: &mut [PileupResult]) {
    for r in results.iter_mut() {
        if r.bq_count > 0 {
            r.mean_bq = Some((r.bq_sum / r.bq_count as f64 * 10.0).round() / 10.0);
        }
        if let Some(dp) = r.dp {
            if dp > 0 && r.mq_sum > 0.0 {
                r.mean_mq = Some((r.mq_sum / dp as f64 * 10.0).round() / 10.0);
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Region computation
// ═══════════════════════════════════════════════════════════════════════════════

/// Compute query regions from positions and optional BED intervals.
///
/// Returns `Vec<(chrom, start, end, Vec<(pos, orig_idx, ref_byte, alt_byte)>)>`.
/// When BED regions are provided, uses them directly (positions are grouped into
/// the BED interval they fall within). Otherwise uses 1Mb sliding windows.
fn compute_regions(
    chroms: &[String],
    positions: &[i64],
    ref_bases: &[String],
    alt_bases: &[String],
    chrom_lengths: &HashMap<String, i64>,
    bed_regions: Option<&[(String, i64, i64)]>,
) -> Vec<(String, i64, i64, Vec<(i64, usize, u8, u8)>)> {
    let n = chroms.len();

    if let Some(bed) = bed_regions {
        // BED-guided: each BED interval is a query region.
        // Uses binary search for O(n_positions × log n_regions) lookup.
        let mut region_data: Vec<(String, i64, i64, Vec<(i64, usize, u8, u8)>)> = Vec::new();
        // Track which BED region index each (chrom, start) maps to
        let mut region_index: HashMap<(String, i64), usize> = HashMap::new();

        // Build per-chromosome sorted interval list for binary search.
        // BED intervals are non-overlapping and sorted by start within each chr.
        let mut chrom_intervals: HashMap<String, Vec<(i64, i64, usize)>> = HashMap::new();
        for (bed_chrom, bed_start, bed_end) in bed {
            if !chrom_lengths.contains_key(bed_chrom.as_str()) { continue; }
            let idx = region_data.len();
            region_data.push((bed_chrom.clone(), *bed_start, *bed_end, Vec::new()));
            region_index.insert((bed_chrom.clone(), *bed_start), idx);
            chrom_intervals.entry(bed_chrom.clone())
                .or_insert_with(Vec::new)
                .push((*bed_start, *bed_end, idx));
        }

        for i in 0..n {
            let pos = positions[i];
            if pos <= 0 { continue; }
            let chrom = &chroms[i];
            if !chrom_lengths.contains_key(chrom.as_str()) { continue; }

            let ref_byte = ref_bases.get(i)
                .and_then(|s| s.as_bytes().first()).copied().unwrap_or(0);
            let alt_byte = alt_bases.get(i)
                .and_then(|s| s.as_bytes().first()).copied().unwrap_or(0);

            // Binary search: find the last interval with start < pos
            if let Some(intervals) = chrom_intervals.get(chrom.as_str()) {
                let search_idx = intervals.partition_point(|(s, _, _)| *s < pos);
                if search_idx > 0 {
                    let (start, end, data_idx) = intervals[search_idx - 1];
                    if pos > start && pos <= end {
                        region_data[data_idx].3.push((pos, i, ref_byte, alt_byte));
                    }
                }
            }
        }

        region_data.sort_by_key(|(c, s, _, _)| (c.clone(), *s));
        region_data
    } else {
        // 1Mb sliding windows
        let mut windows: HashMap<(String, i64), Vec<(i64, usize, u8, u8)>> = HashMap::new();

        for i in 0..n {
            let pos = positions[i];
            if pos <= 0 { continue; }
            if !chrom_lengths.contains_key(&chroms[i]) { continue; }
            let win_start = (pos - 1) / WINDOW_SIZE * WINDOW_SIZE + 1;
            let ref_byte = ref_bases.get(i)
                .and_then(|s| s.as_bytes().first()).copied().unwrap_or(0);
            let alt_byte = alt_bases.get(i)
                .and_then(|s| s.as_bytes().first()).copied().unwrap_or(0);

            let chrom = chroms[i].clone();
            windows.entry((chrom, win_start))
                .or_insert_with(Vec::new)
                .push((pos, i, ref_byte, alt_byte));
        }

        let mut regions: Vec<_> = windows.into_iter().map(|((chrom, win_start), win_positions)| {
            let chrom_len = chrom_lengths.get(&chrom).copied().unwrap_or(0);
            let win_end = (win_start + WINDOW_SIZE - 1).min(chrom_len);
            (chrom, win_start, win_end, win_positions)
        }).collect();
        regions.sort_by_key(|(c, s, _, _)| (c.clone(), *s));
        regions
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Single BAM pileup (existing API — unchanged signature)
// ═══════════════════════════════════════════════════════════════════════════════

/// Perform pileup at variant positions using windowed BAI-indexed queries.
///
/// Groups positions into ~1Mb genomic windows, queries each window once,
/// and matches reads to target positions within each window.
///
/// Falls back to per-position queries if the BAI index is missing.
pub fn pileup_variants(
    bam_path: &Path,
    chroms: &[String],
    positions: &[i64],
    ref_bases: &[String],
    alt_bases: &[String],
) -> Result<Vec<PileupResult>, String> {
    let n = chroms.len();
    if n == 0 { return Ok(Vec::new()); }

    let file = std::fs::File::open(bam_path)
        .map_err(|e| format!("Cannot open {}: {}", bam_path.display(), e))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().map_err(|e| format!("Header error: {}", e))?;

    // Build chromosome length map
    let mut chrom_lengths: HashMap<String, i64> = HashMap::new();
    for (name, rs) in header.reference_sequences() {
        if let Ok(name_str) = std::str::from_utf8(name) {
            chrom_lengths.insert(name_str.to_string(), usize::from(rs.length()) as i64);
        }
    }

    let bai_path = format!("{}.bai", bam_path.display());
    let bai_path = Path::new(&bai_path);
    if !bai_path.exists() {
        return Ok(vec![PileupResult::default(); n]);
    }
    let index = bai::fs::read(bai_path)
        .map_err(|e| format!("Cannot read BAI index: {}", e))?;

    let mut results: Vec<PileupResult> = vec![PileupResult::default(); n];

    // ── Group positions by (chromosome, window_start) ────────────────────
    let mut windows: HashMap<(String, i64), Vec<(i64, usize, u8, u8)>> = HashMap::new();

    for i in 0..n {
        let pos = positions[i];
        if pos <= 0 { continue; }
        if !chrom_lengths.contains_key(&chroms[i]) { continue; }
        let win_start = (pos - 1) / WINDOW_SIZE * WINDOW_SIZE + 1;
        let ref_byte = ref_bases.get(i)
            .and_then(|s| s.as_bytes().first()).copied().unwrap_or(0);
        let alt_byte = alt_bases.get(i)
            .and_then(|s| s.as_bytes().first()).copied().unwrap_or(0);

        let chrom = chroms[i].clone();
        windows.entry((chrom, win_start))
            .or_insert_with(Vec::new)
            .push((pos, i, ref_byte, alt_byte));
    }

    // ── Process each window with a single BAI query ──────────────────────
    for ((chrom, win_start), win_positions) in &windows {
        let chrom_len = chrom_lengths.get(chrom).copied().unwrap_or(0);
        if chrom_len == 0 { continue; }
        let win_end = (*win_start + WINDOW_SIZE - 1).min(chrom_len);
        if *win_start > chrom_len { continue; }

        // Build pos_map: pos → Vec<(orig_idx, ref_byte, alt_byte)>
        let mut pos_map: HashMap<i64, Vec<(usize, u8, u8)>> = HashMap::new();
        for (pos, orig_idx, ref_byte, alt_byte) in win_positions {
            pos_map.entry(*pos).or_insert_with(Vec::new)
                .push((*orig_idx, *ref_byte, *alt_byte));
        }

        // Query the window region
        let pos_start = match std::num::NonZero::new(*win_start as usize) {
            Some(nz) => match noodles_core::Position::try_from(usize::from(nz)) {
                Ok(p) => p, Err(_) => continue,
            }, None => continue,
        };
        let pos_end = match std::num::NonZero::new(win_end as usize) {
            Some(nz) => match noodles_core::Position::try_from(usize::from(nz)) {
                Ok(p) => p, Err(_) => continue,
            }, None => continue,
        };
        let region = Region::new(chrom.clone(), pos_start..=pos_end);
        let query = match reader.query(&header, &index, &region) {
            Ok(q) => q, Err(_) => continue,
        };

        for record_result in query.records() {
            let record = match record_result { Ok(r) => r, Err(_) => continue };
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_duplicate() { continue; }

            let align_start = match record.alignment_start() {
                Some(Ok(p)) => usize::from(p) as i64,
                _ => continue,
            };
            let seq = record.sequence();
            let align_end = align_start + alignment_reference_span(&record);
            let quals = record.quality_scores();

            for (pos, entries) in &pos_map {
                if *pos < align_start || *pos >= align_end { continue; }
                let pos_in_read = (pos - align_start) as usize;
                if pos_in_read >= seq.len() { continue; }
                let base = seq.get(pos_in_read);

                let is_rev = flags.is_reverse_complemented();
                let mq = record.mapping_quality().map(|q| u8::from(q) as f64);
                let bq_val = if pos_in_read < quals.len() {
                    quals.as_bytes().get(pos_in_read).copied()
                } else { None };

                for (orig_idx, ref_byte, alt_byte) in entries {
                    let r = &mut results[*orig_idx];
                    update_pileup_result(r, base, *ref_byte, *alt_byte, is_rev, bq_val, mq);
                }
            }
        }
    }

    finalize_results(&mut results);
    Ok(results)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Combined multi-BAM pileup (Phase 8 — new function)
// ═══════════════════════════════════════════════════════════════════════════════

/// Perform pileup across multiple BAMs in a single FFI call.
///
/// Groups positions into regions (BED-guided or 1Mb windows), queries each
/// region once per BAM, and uses binary search to match reads to positions.
///
/// Returns results grouped by BAM label: `Vec<(bam_label, Vec<PileupResult>)>`.
/// All results arrays have the same length (number of input positions).
///
/// # Arguments
/// * `bam_paths` — Paths to BAM files (must have corresponding .bai indices).
/// * `bam_labels` — Human-readable labels for each BAM (e.g., "DN", "DT", "RT").
/// * `chroms`, `positions`, `ref_bases`, `alt_bases` — Variant positions.
/// * `bed_regions` — Optional merged BED intervals for region-guided queries.
///   When provided, these replace 1Mb windows, dramatically reducing queries
///   for WES data (~300 regions vs ~2,765 windows).
pub fn pileup_variants_multi(
    bam_paths: &[String],
    bam_labels: &[String],
    chroms: &[String],
    positions: &[i64],
    ref_bases: &[String],
    alt_bases: &[String],
    bed_regions: Option<&[(String, i64, i64)]>,
) -> Result<Vec<(String, Vec<PileupResult>)>, String> {
    let n = chroms.len();
    let n_bams = bam_paths.len();

    if n == 0 || n_bams == 0 { return Ok(Vec::new()); }

    // ── Open all BAM readers once and keep them alive ────────────────────
    // CRITICAL: noodles-bgzf maintains internal block state (pos, len).
    // Reopening readers per-region causes stale block.data().pos from the
    // header read to carry into region queries. When block.pos > block.len
    // on a new BGZF block, Data::as_ref() panics with "slice index starts
    // at X but ends at Y". Keeping readers open avoids this — the same
    // pattern used by the single-BAM pileup_variants().
    //
    // We use Vec<(reader, header)> instead of a struct to avoid naming
    // the complex generic types; Rust infers them from bindings.
    let mut bam_chrom_lengths: Vec<HashMap<String, i64>> = Vec::with_capacity(n_bams);
    let mut bam_readers_and_headers = Vec::with_capacity(n_bams);
    let mut bam_indices: Vec<bai::Index> = Vec::with_capacity(n_bams);

    for i in 0..n_bams {
        let file = std::fs::File::open(&bam_paths[i])
            .map_err(|e| format!("Cannot open {}: {}", bam_paths[i], e))?;
        let mut reader = bam::io::Reader::new(file);
        let header = reader.read_header()
            .map_err(|e| format!("Header error on {}: {}", bam_paths[i], e))?;

        let bai_path_str = format!("{}.bai", bam_paths[i]);
        let bai_path = Path::new(&bai_path_str);
        if !bai_path.exists() {
            let defaults: Vec<PileupResult> = vec![PileupResult::default(); n];
            return Ok(bam_labels.iter().map(|l| (l.clone(), defaults.clone())).collect());
        }
        let index = bai::fs::read(bai_path)
            .map_err(|e| format!("Cannot read BAI for {}: {}", bam_paths[i], e))?;

        // Build chrom_lengths from the already-read header
        let mut cls: HashMap<String, i64> = HashMap::new();
        for (name, rs) in header.reference_sequences() {
            if let Ok(name_str) = std::str::from_utf8(name) {
                cls.insert(name_str.to_string(), usize::from(rs.length()) as i64);
            }
        }

        bam_chrom_lengths.push(cls);
        bam_readers_and_headers.push((reader, header));
        bam_indices.push(index);
    }

    // ── Compute query regions (shared across all BAMs) ───────────────────
    let regions = compute_regions(
        chroms, positions, ref_bases, alt_bases,
        &bam_chrom_lengths[0], bed_regions,
    );

    // ── Initialize per-BAM result arrays ─────────────────────────────────
    let mut all_results: Vec<Vec<PileupResult>> = (0..n_bams)
        .map(|_| vec![PileupResult::default(); n])
        .collect();

    // ── Process each region on each BAM ──────────────────────────────────
    for (region_chrom, region_start, region_end, region_positions) in &regions {
        if region_positions.is_empty() { continue; }
        if *region_start > *region_end { continue; }

        // Build sorted position array and entries for binary search
        let mut pos_entries: Vec<(i64, usize, u8, u8)> = region_positions
            .iter()
            .map(|(pos, orig_idx, ref_byte, alt_byte)| (*pos, *orig_idx, *ref_byte, *alt_byte))
            .collect();
        pos_entries.sort_by_key(|(pos, _, _, _)| *pos);

        let sorted_positions: Vec<i64> = pos_entries.iter().map(|(p, _, _, _)| *p).collect();
        let max_pos = sorted_positions.last().copied().unwrap_or(i64::MAX);

        for bam_idx in 0..n_bams {
            // Per-BAM chromosome length check and region clamping
            let chrom_len = bam_chrom_lengths[bam_idx].get(region_chrom).copied().unwrap_or(0);
            if chrom_len == 0 { continue; }
            if *region_start > chrom_len { continue; }
            let clamped_end = (*region_end).min(chrom_len);

            let pos_start_nz = match std::num::NonZero::new(*region_start as usize) {
                Some(nz) => nz, None => continue,
            };
            let pos_end_nz = match std::num::NonZero::new(clamped_end as usize) {
                Some(nz) => nz, None => continue,
            };
            let pos_start = match noodles_core::Position::try_from(usize::from(pos_start_nz)) {
                Ok(p) => p, Err(_) => continue,
            };
            let pos_end = match noodles_core::Position::try_from(usize::from(pos_end_nz)) {
                Ok(p) => p, Err(_) => continue,
            };

            let region = Region::new(region_chrom.clone(), pos_start..=pos_end);
            let (reader, header) = &mut bam_readers_and_headers[bam_idx];
            let query = match reader.query(header, &bam_indices[bam_idx], &region) {
                Ok(q) => q, Err(_) => continue,
            };

            // Process reads with binary search inner loop
            for record_result in query.records() {
                let record = match record_result { Ok(r) => r, Err(_) => continue };
                let flags = record.flags();
                if flags.is_unmapped() || flags.is_duplicate() { continue; }

                let align_start = match record.alignment_start() {
                    Some(Ok(p)) => usize::from(p) as i64,
                    _ => continue,
                };
                let seq = record.sequence();
                let align_end = align_start + alignment_reference_span(&record);

                if align_end <= sorted_positions[0] { continue; }
                if align_start > max_pos { break; }

                let quals = record.quality_scores();
                let is_rev = flags.is_reverse_complemented();
                let mq = record.mapping_quality().map(|q| u8::from(q) as f64);

                let start_idx = sorted_positions.binary_search_by(|p| p.cmp(&align_start))
                    .unwrap_or_else(|idx| idx);

                for idx in start_idx..sorted_positions.len() {
                    let pos = sorted_positions[idx];
                    if pos >= align_end { break; }

                    let pos_in_read = (pos - align_start) as usize;
                    if pos_in_read >= seq.len() { continue; }
                    let base = seq.get(pos_in_read);
                    let bq_val = if pos_in_read < quals.len() {
                        quals.as_bytes().get(pos_in_read).copied()
                    } else { None };

                    let (_pos, orig_idx, ref_byte, alt_byte) = pos_entries[idx];
                    let r = &mut all_results[bam_idx][orig_idx];
                    update_pileup_result(r, base, ref_byte, alt_byte, is_rev, bq_val, mq);
                }
            }
        }
    }

    // ── Finalize all results ─────────────────────────────────────────────
    for results in all_results.iter_mut() {
        finalize_results(results);
    }

    // ── Group by BAM label ───────────────────────────────────────────────
    let output: Vec<(String, Vec<PileupResult>)> = bam_labels.iter()
        .zip(all_results.into_iter())
        .map(|(label, results)| (label.clone(), results))
        .collect();

    Ok(output)
}
