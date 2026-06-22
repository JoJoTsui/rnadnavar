//! BAM parsing using noodles-bam 0.90.
//!
//! Provides whole-genome BAM statistics with optional BED-guided on-target
//! coverage calculation for WES data.

use std::collections::HashMap;
use std::path::Path;

use noodles_bam as bam;

/// Compute the reference alignment span from CIGAR operations.
///
/// Sums lengths of reference-consuming ops (M, D, N, =, X).
/// This is the correct alignment_end offset — record.sequence().len()
/// includes insertions and excludes deletions/skipped regions.
fn alignment_reference_span(record: &bam::Record) -> i64 {
    let mut span: i64 = 0;
    for op_result in record.cigar().iter() {
        if let Ok(op) = op_result {
            if op.kind().consumes_reference() {
                span += op.len() as i64;
            }
        }
    }
    // Fallback to sequence length if CIGAR is empty (unmapped reads, etc.)
    if span == 0 {
        record.sequence().len() as i64
    } else {
        span
    }
}

/// Whole-genome BAM statistics.
#[derive(Debug, Clone, Default)]
pub struct BamStats {
    pub total_reads: u64,
    pub mapped_reads: u64,
    pub mapping_rate: f64,
    pub mean_coverage: f64,
    pub mean_insert_size: f64,
    pub mean_mapq: f64,
    pub duplication_rate_pct: f64,
    pub properly_paired_pct: f64,
    pub insert_size_stddev: f64,
}

/// Compute whole-genome BAM statistics (sample up to max_reads, 0 = no limit).
pub fn whole_genome_stats(bam_path: &Path, max_reads: u64) -> Result<BamStats, Box<dyn std::error::Error>> {
    whole_genome_stats_impl(bam_path, max_reads, None, 0)
}

/// Compute BAM statistics with BED-guided on-target coverage.
///
/// When `bed_regions` is provided, only bases overlapping BED intervals are
/// counted toward mean_coverage. Other metrics (total reads, mapping rate,
/// insert size, MAPQ) are unaffected — they remain whole-genome.
pub fn whole_genome_stats_bed(
    bam_path: &Path,
    max_reads: u64,
    bed_chroms: &[String],
    bed_starts: &[i64],
    bed_ends: &[i64],
) -> Result<BamStats, Box<dyn std::error::Error>> {
    let n = bed_chroms.len();
    let regions: Vec<(String, i64, i64)> = (0..n)
        .map(|i| (bed_chroms[i].clone(), bed_starts[i], bed_ends[i]))
        .collect();
    let bed_total: u64 = regions.iter().map(|(_, s, e)| (e - s) as u64).sum();
    whole_genome_stats_impl(bam_path, max_reads, Some(&regions), bed_total)
}

/// Shared implementation: whole-genome scan with optional BED on-target tracking.
fn whole_genome_stats_impl(
    bam_path: &Path,
    max_reads: u64,
    bed_regions: Option<&Vec<(String, i64, i64)>>,
    bed_total: u64,
) -> Result<BamStats, Box<dyn std::error::Error>> {
    let mut reader = bam::io::Reader::new(std::fs::File::open(bam_path)?);
    let header = reader.read_header()?;

    // Sum of all reference sequence lengths for coverage estimation
    let ref_lengths: u64 = header.reference_sequences().values()
        .map(|rs| usize::from(rs.length()) as u64)
        .sum();

    // Build per-chromosome BED interval index for O(reads + regions) lookup.
    let bed_index: Option<HashMap<String, Vec<(i64, i64)>>> = bed_regions.map(|regions| {
        let mut map: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        for (chrom, start, end) in regions {
            map.entry(chrom.clone()).or_insert_with(Vec::new).push((*start, *end));
        }
        map
    });

    // Build reference ID → chromosome name mapping for on-target check
    let ref_id_to_name: Option<HashMap<usize, String>> = bed_index.as_ref().map(|_| {
        let mut map = HashMap::new();
        let rs = header.reference_sequences();
        for (idx, (name, _)) in rs.iter().enumerate() {
            if let Ok(name_str) = std::str::from_utf8(name) {
                map.insert(idx, name_str.to_string());
            }
        }
        map
    });

    let mut total: u64 = 0;
    let mut mapped: u64 = 0;
    let mut mq_sum: f64 = 0.0;
    let mut insert_sum: f64 = 0.0;
    let mut insert_sum_sq: f64 = 0.0;
    let mut insert_count: u64 = 0;
    let mut total_query_length: u64 = 0;
    let mut on_target_bases: u64 = 0;
    let mut dup_count: u64 = 0;
    let mut proper_pair_count: u64 = 0;

    // Per-chromosome BED cursor for two-pointer walk
    let mut bed_cursors: HashMap<String, usize> = HashMap::new();

    for result in reader.records() {
        let record = match result { Ok(r) => r, Err(_) => continue };
        total += 1;

        let flags = record.flags();
        if flags.is_duplicate() { dup_count += 1; }
        if !flags.is_unmapped() {
            mapped += 1;
            if flags.is_properly_segmented() { proper_pair_count += 1; }
            if let Some(mq) = record.mapping_quality() {
                mq_sum += u8::from(mq) as f64;
            }
            let seq_len = record.sequence().len() as u64;
            total_query_length += seq_len;

            // On-target check: count bases overlapping BED regions
            if let Some(ref index) = bed_index {
                let chrom = ref_id_to_name.as_ref().and_then(|id_map| {
                    record.reference_sequence_id()
                        .and_then(|r| r.ok())
                        .and_then(|id| id_map.get(&id).cloned())
                });
                if let Some(ref chrom_name) = chrom {
                    if let Some(intervals) = index.get(chrom_name) {
                        let align_start = match record.alignment_start() {
                            Some(Ok(p)) => usize::from(p) as i64,
                            _ => 0,
                        };
                        let align_end = align_start + alignment_reference_span(&record);
                        // Two-pointer walk: advance cursor past intervals before this read
                        let cursor = bed_cursors.entry(chrom_name.clone()).or_insert(0);
                        while *cursor < intervals.len() && intervals[*cursor].1 <= align_start {
                            *cursor += 1;
                        }
                        // Check current and subsequent intervals for overlap
                        let mut i = if *cursor > 0 { *cursor - 1 } else { 0 };
                        while i < intervals.len() && intervals[i].0 < align_end {
                            let (int_start, int_end) = intervals[i];
                            if align_start < int_end && align_end > int_start {
                                let overlap_start = align_start.max(int_start);
                                let overlap_end = align_end.min(int_end);
                                if overlap_end > overlap_start {
                                    on_target_bases += (overlap_end - overlap_start) as u64;
                                }
                            }
                            i += 1;
                        }
                    }
                }
            }

            // Insert size: only count properly paired reads with positive TLEN.
            // Supplementary/improper pairs can have arbitrarily large TLEN values
            // that would bias the mean.
            if flags.is_properly_segmented()
                && !flags.is_supplementary()
                && !flags.is_secondary()
            {
                let tlen = record.template_length();
                if tlen > 0 {
                    insert_sum += tlen as f64;
                    insert_sum_sq += (tlen as f64) * (tlen as f64);
                    insert_count += 1;
                }
            }
        }

        if max_reads > 0 && total >= max_reads { break; }
    }

    // Coverage: use on-target bases when BED is provided, else total bases
    let coverage_bases = if bed_regions.is_some() && bed_total > 0 {
        on_target_bases
    } else {
        total_query_length
    };
    let denominator = if bed_total > 0 { bed_total } else { ref_lengths };

    let duplication_rate_pct = if total > 0 { 100.0 * dup_count as f64 / total as f64 } else { 0.0 };
    let properly_paired_pct = if mapped > 0 { 100.0 * proper_pair_count as f64 / mapped as f64 } else { 0.0 };
    let insert_size_stddev = if insert_count > 0 {
        let mean = insert_sum / insert_count as f64;
        let variance = (insert_sum_sq / insert_count as f64) - (mean * mean);
        if variance > 0.0 { variance.sqrt() } else { 0.0 }
    } else { 0.0 };

    Ok(BamStats {
        total_reads: total,
        mapped_reads: mapped,
        mapping_rate: if total > 0 { mapped as f64 / total as f64 * 100.0 } else { 0.0 },
        mean_coverage: if denominator > 0 { coverage_bases as f64 / denominator as f64 } else { 0.0 },
        mean_insert_size: if insert_count > 0 { insert_sum / insert_count as f64 } else { 0.0 },
        mean_mapq: if mapped > 0 { mq_sum / mapped as f64 } else { 0.0 },
        duplication_rate_pct,
        properly_paired_pct,
        insert_size_stddev,
    })
}

/// Compute per-base coverage depth bins across BED regions.
///
/// Returns a map of "{threshold}x_pct" → percentage of on-target bases
/// covered at that depth or higher. Thresholds: 1x, 10x, 20x, 50x, 100x.
pub fn coverage_bins(bam_path: &Path, bed_regions: &[(String, u32, u32)]) -> Result<Option<HashMap<String, f64>>, String> {
    use bam::bai;
    use noodles_core::Region;

    if bed_regions.is_empty() { return Ok(None); }

    let mut reader = bam::io::Reader::new(std::fs::File::open(bam_path).map_err(|e| e.to_string())?);
    let header = reader.read_header().map_err(|e| e.to_string())?;

    // Load BAI index
    let bai_path_str = format!("{}.bai", bam_path.display());
    let bai_path = Path::new(&bai_path_str);
    if !bai_path.exists() {
        return Err(format!("BAI index not found: {}", bai_path.display()));
    }
    let index = bai::fs::read(bai_path)
        .map_err(|e| format!("Cannot read BAI index: {}", e))?;

    let mut total_bases: u64 = 0;
    let mut bases_1x: u64 = 0; let mut bases_10x: u64 = 0;
    let mut bases_20x: u64 = 0; let mut bases_50x: u64 = 0; let mut bases_100x: u64 = 0;

    for (chrom, start, end) in bed_regions {
        let region_len = (end - start) as u64;
        total_bases += region_len;
        let mut depths = vec![0u32; region_len as usize];

        let pos_start = std::num::NonZero::new(*start as usize)
            .and_then(|nz| noodles_core::Position::try_from(usize::from(nz)).ok());
        let pos_end = std::num::NonZero::new(*end as usize)
            .and_then(|nz| noodles_core::Position::try_from(usize::from(nz)).ok());
        let (Some(pos_start), Some(pos_end)) = (pos_start, pos_end) else { continue; };
        let region = Region::new(chrom.clone(), pos_start..=pos_end);

        let query = reader.query(&header, &index, &region)
            .map_err(|e| e.to_string())?;

        for record_result in query.records() {
            let record = match record_result { Ok(r) => r, Err(_) => continue };
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_duplicate() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            let rec_start = match record.alignment_start() {
                Some(Ok(p)) => usize::from(p) as i64,
                _ => continue,
            };
            let cigar = record.cigar();
            let mut pos = rec_start;
            for op_result in cigar.iter() {
                let op = match op_result { Ok(o) => o, Err(_) => continue };
                if op.kind().consumes_reference() {
                    let op_len = op.len() as i64;
                    // Only accumulate depth for alignment matches (M, =, X),
                    // not for deletions or skipped reference regions (D, N).
                    use noodles_sam::alignment::record::cigar::op::Kind;
                    let is_alignment = matches!(op.kind(), Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch);
                    if is_alignment {
                        for offset in 0..op_len {
                            let base_pos = pos + offset - (*start as i64);
                            if base_pos >= 0 && (base_pos as u64) < region_len {
                                depths[base_pos as usize] += 1;
                            }
                        }
                    }
                    pos += op_len;
                }
            }
        }
        for d in &depths {
            if *d >= 1 { bases_1x += 1; }
            if *d >= 10 { bases_10x += 1; }
            if *d >= 20 { bases_20x += 1; }
            if *d >= 50 { bases_50x += 1; }
            if *d >= 100 { bases_100x += 1; }
        }
    }
    let mut map = HashMap::new();
    if total_bases > 0 {
        map.insert("cov_1x_pct".to_string(), 100.0 * bases_1x as f64 / total_bases as f64);
        map.insert("cov_10x_pct".to_string(), 100.0 * bases_10x as f64 / total_bases as f64);
        map.insert("cov_20x_pct".to_string(), 100.0 * bases_20x as f64 / total_bases as f64);
        map.insert("cov_50x_pct".to_string(), 100.0 * bases_50x as f64 / total_bases as f64);
        map.insert("cov_100x_pct".to_string(), 100.0 * bases_100x as f64 / total_bases as f64);
    }
    Ok(Some(map))
}
