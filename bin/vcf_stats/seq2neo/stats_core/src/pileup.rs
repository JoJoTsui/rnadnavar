//! BAM pileup at variant positions using noodles-bam indexed queries.

use std::path::Path;

use noodles_bam::{self as bam, bai};
use noodles_core::Region;

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
}

/// Perform pileup at variant positions using BAI-indexed queries.
pub fn pileup_variants(
    bam_path: &Path,
    chroms: &[String],
    positions: &[i64],
    ref_bases: &[String],
    alt_bases: &[String],
) -> Result<Vec<PileupResult>, String> {
    let n = chroms.len();
    if n == 0 { return Ok(Vec::new()); }

    // Open BAM
    let file = std::fs::File::open(bam_path)
        .map_err(|e| format!("Cannot open {}: {}", bam_path.display(), e))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().map_err(|e| format!("Header error: {}", e))?;

    // Read BAI index
    let bai_path = format!("{}.bai", bam_path.display());
    let bai_path = Path::new(&bai_path);
    if !bai_path.exists() {
        // No index — all positions return None
        return Ok(vec![PileupResult::default(); n]);
    }
    let index = bai::fs::read(bai_path)
        .map_err(|e| format!("Cannot read BAI index: {}", e))?;

    // Sort positions by (chrom, pos) for efficient chunked access
    let mut indexed: Vec<(usize, &str, i64, u8, u8)> = (0..n)
        .map(|i| {
            let ref_byte = ref_bases.get(i).and_then(|s| s.as_bytes().first()).copied().unwrap_or(0);
            let alt_byte = alt_bases.get(i).and_then(|s| s.as_bytes().first()).copied().unwrap_or(0);
            (i, chroms[i].as_str(), positions[i], ref_byte, alt_byte)
        })
        .collect();
    indexed.sort_by(|a, b| a.1.cmp(b.1).then(a.2.cmp(&b.2)));

    let mut results: Vec<PileupResult> = vec![PileupResult::default(); n];

    for (orig_idx, chrom, pos, ref_byte, alt_byte) in &indexed {
        if *pos <= 0 || chrom.is_empty() { continue; }

        let pos_usize = *pos as usize;
        let pos_nz = match std::num::NonZero::new(pos_usize) {
            Some(nz) => nz,
            None => continue,
        };
        let pos_val = match noodles_core::Position::try_from(usize::from(pos_nz)) {
            Ok(p) => p,
            Err(_) => continue,
        };
        let region = Region::new(chrom.to_string(), pos_val..=pos_val);

        let query = match reader.query(&header, &index, &region) {
            Ok(q) => q,
            Err(_) => continue,
        };

        let mut dp: i64 = 0;
        let mut ref_dp: i64 = 0;
        let mut alt_dp: i64 = 0;
        let mut f1r2_ref: i64 = 0;
        let mut f2r1_ref: i64 = 0;
        let mut f1r2_alt: i64 = 0;
        let mut f2r1_alt: i64 = 0;
        let mut bq_sum: f64 = 0.0;
        let mut bq_count: i64 = 0;
        let mut mq_sum: f64 = 0.0;

        for record_result in query.records() {
            let record = match record_result { Ok(r) => r, Err(_) => continue };
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_duplicate() { continue; }

            let alignment_start = match record.alignment_start() {
                Some(Ok(p)) => usize::from(p) as i64,
                _ => continue,
            };

            let pos_in_read = *pos - alignment_start;
            if pos_in_read < 0 { continue; }

            let seq = record.sequence();
            if pos_in_read as usize >= seq.len() { continue; }

            let base = seq.as_bytes().get(pos_in_read as usize).copied();

            dp += 1;
            if let Some(mq) = record.mapping_quality() { mq_sum += u8::from(mq) as f64; }
            let quals = record.quality_scores();
            if pos_in_read < quals.len() as i64 {
                let qual = quals.as_bytes().get(pos_in_read as usize).copied().unwrap_or(0);
                bq_sum += qual as f64;
                bq_count += 1;
            }

            let is_rev = flags.is_reverse_complemented();
            if let Some(b) = base {
                if b == *ref_byte && *ref_byte != 0 {
                    ref_dp += 1;
                    if is_rev { f2r1_ref += 1 } else { f1r2_ref += 1 }
                } else if b == *alt_byte && *alt_byte != 0 {
                    alt_dp += 1;
                    if is_rev { f2r1_alt += 1 } else { f1r2_alt += 1 }
                }
            }
        }

        results[*orig_idx] = PileupResult {
            dp: Some(dp), ref_dp: Some(ref_dp), alt_dp: Some(alt_dp),
            f1r2_ref: Some(f1r2_ref), f2r1_ref: Some(f2r1_ref),
            f1r2_alt: Some(f1r2_alt), f2r1_alt: Some(f2r1_alt),
            mean_bq: if bq_count > 0 { Some((bq_sum / bq_count as f64 * 10.0).round() / 10.0) } else { None },
            mean_mq: if dp > 0 { Some((mq_sum / dp as f64 * 10.0).round() / 10.0) } else { None },
        };
    }

    Ok(results)
}
