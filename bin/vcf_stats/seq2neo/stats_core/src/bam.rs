//! BAM parsing using noodles-bam 0.90.
//!
//! NOTE: Full BAM pileup with per-position strand/quality analysis is deferred
//! due to noodles-sam version conflicts in the dependency chain. The Python
//! pysam fallback (rust_bam.py) handles pileup. This module provides only
//! whole-genome BAM statistics.

use std::path::Path;

use noodles_bam as bam;

/// Whole-genome BAM statistics.
#[derive(Debug, Clone, Default)]
pub struct BamStats {
    pub total_reads: u64,
    pub mapped_reads: u64,
    pub mapping_rate: f64,
    pub mean_insert_size: f64,
    pub mean_mapq: f64,
}

/// Compute whole-genome BAM statistics (sample up to max_reads).
pub fn whole_genome_stats(bam_path: &Path, max_reads: u64) -> Result<BamStats, Box<dyn std::error::Error>> {
    let mut reader = bam::io::Reader::new(std::fs::File::open(bam_path)?);
    reader.read_header()?;

    let mut total: u64 = 0;
    let mut mapped: u64 = 0;
    let mut mq_sum: f64 = 0.0;
    let mut insert_sum: f64 = 0.0;
    let mut insert_count: u64 = 0;

    for result in reader.records() {
        let record = match result { Ok(r) => r, Err(_) => continue };
        total += 1;
        if !record.flags().is_unmapped() {
            mapped += 1;
            if let Some(mq) = record.mapping_quality() { mq_sum += u8::from(mq) as f64; }
            let tlen = record.template_length();
            if tlen != 0 {
                let t = tlen.unsigned_abs() as u64;
                insert_sum += t as f64;
                insert_count += 1;
            }
        }
        if max_reads > 0 && total >= max_reads { break; }
    }

    Ok(BamStats {
        total_reads: total,
        mapped_reads: mapped,
        mapping_rate: if total > 0 { mapped as f64 / total as f64 * 100.0 } else { 0.0 },
        mean_insert_size: if insert_count > 0 { insert_sum / insert_count as f64 } else { 0.0 },
        mean_mapq: if mapped > 0 { mq_sum / mapped as f64 } else { 0.0 },
    })
}
