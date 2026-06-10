//! Caller VCF parsing using noodles-vcf 0.88.
//!
//! Parses normalized caller VCFs (*.dec.norm.vcf.gz) to extract FORMAT fields
//! at target positions matching (CHROM, POS, REF, ALT) against the rescue VCF.

use std::collections::HashSet;
use std::io::BufReader;
use std::path::Path;

use noodles_bgzf as bgzf;
use noodles_vcf as vcf;
use noodles_vcf::variant::record::samples::series::Value;
use noodles_vcf::variant::record::AlternateBases as _;

use vcf::variant::record::samples::Sample as SampleTrait;

/// Describes which FORMAT fields to extract for a given caller type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CallerKind {
    Mutect2,
    DeepSomatic,
    Strelka,
}

impl CallerKind {
    pub fn from_name(name: &str) -> Self {
        let lower = name.to_lowercase();
        if lower.contains("strelka") { CallerKind::Strelka }
        else if lower.contains("mutect2") { CallerKind::Mutect2 }
        else { CallerKind::DeepSomatic }
    }
    pub fn has_gt_ad(&self) -> bool { matches!(self, CallerKind::Mutect2 | CallerKind::DeepSomatic) }
    pub fn has_strand_fields(&self) -> bool { matches!(self, CallerKind::Mutect2) }
}

/// Holds extracted FORMAT field values for matched variants as column vectors.
#[derive(Debug, Clone, Default)]
pub struct CallerResults {
    pub chroms: Vec<String>,
    pub positions: Vec<i64>,
    pub refs: Vec<String>,
    pub alts: Vec<String>,
    pub dp: Vec<Option<i64>>,
    pub ad_ref: Vec<Option<i64>>,
    pub ad_alt: Vec<Option<i64>>,
    pub gt: Vec<Option<String>>,
    pub vaf_caller: Vec<Option<f64>>,
    pub tar: Vec<Option<i64>>,
    pub tir: Vec<Option<i64>>,
    pub tor: Vec<Option<i64>>,
    pub au: Vec<Option<i64>>,
    pub cu: Vec<Option<i64>>,
    pub gu: Vec<Option<i64>>,
    pub tu: Vec<Option<i64>>,
    pub sb: Vec<Option<String>>,
    pub fad: Vec<Option<String>>,
}

// ── FORMAT field extraction helpers ──────────────────────────────────────

fn get_int(s: &impl SampleTrait, h: &vcf::Header, key: &str) -> Option<i64> {
    match s.get(h, key) {
        Some(Ok(Some(Value::Integer(n)))) => Some(n as i64),
        Some(Ok(Some(Value::Float(f)))) => Some(f as i64),
        Some(Ok(Some(Value::String(v)))) => v.parse::<i64>().ok(),
        Some(Ok(Some(Value::Array(arr)))) => match arr {
            vcf::variant::record::samples::series::value::Array::Integer(vals) =>
                vals.iter().next().and_then(|r| r.ok().flatten()).map(|v| v as i64),
            vcf::variant::record::samples::series::value::Array::Float(vals) =>
                vals.iter().next().and_then(|r| r.ok().flatten()).map(|v| v as i64),
            _ => None,
        },
        _ => None,
    }
}

fn get_float(s: &impl SampleTrait, h: &vcf::Header, key: &str) -> Option<f64> {
    match s.get(h, key) {
        Some(Ok(Some(Value::Float(f)))) => Some(f as f64),
        Some(Ok(Some(Value::Integer(n)))) => Some(n as f64),
        Some(Ok(Some(Value::String(v)))) => v.parse::<f64>().ok(),
        // AF/VAF may be parsed as single-element Float array
        Some(Ok(Some(Value::Array(arr)))) => {
            match arr {
                vcf::variant::record::samples::series::value::Array::Float(vals) =>
                    vals.iter().next().and_then(|r| r.ok().flatten()).map(|v| v as f64),
                vcf::variant::record::samples::series::value::Array::Integer(vals) =>
                    vals.iter().next().and_then(|r| r.ok().flatten()).map(|v| v as f64),
                _ => None,
            }
        }
        _ => None,
    }
}

fn get_str(s: &impl SampleTrait, h: &vcf::Header, key: &str) -> Option<String> {
    match s.get(h, key) {
        Some(Ok(Some(Value::String(v)))) => {
            if v.is_empty() { None } else { Some(v.to_string()) }
        }
        Some(Ok(Some(Value::Integer(n)))) => Some(n.to_string()),
        Some(Ok(Some(Value::Float(f)))) => Some(f.to_string()),
        Some(Ok(Some(Value::Character(c)))) => Some(c.to_string()),
        Some(Ok(Some(Value::Genotype(g)))) => {
            let s = format!("{:?}", g);
            if s.is_empty() { None } else { Some(s) }
        }
        Some(Ok(Some(Value::Array(arr)))) => {
            let s = format!("{:?}", arr);
            if s.is_empty() { None } else { Some(s) }
        }
        _ => None,
    }
}

fn get_int_first(s: &impl SampleTrait, h: &vcf::Header, key: &str) -> Option<i64> {
    match s.get(h, key) {
        Some(Ok(Some(Value::String(v)))) => v.split(',').next().and_then(|x| x.parse::<i64>().ok()),
        Some(Ok(Some(Value::Integer(n)))) => Some(n as i64),
        Some(Ok(Some(Value::Float(f)))) => Some(f as i64),
        Some(Ok(Some(Value::Array(arr)))) => match arr {
            vcf::variant::record::samples::series::value::Array::Integer(vals) =>
                vals.iter().next().and_then(|r| r.ok().flatten()).map(|v| v as i64),
            vcf::variant::record::samples::series::value::Array::Float(vals) =>
                vals.iter().next().and_then(|r| r.ok().flatten()).map(|v| v as i64),
            _ => None,
        },
        _ => None,
    }
}

/// Parse a normalized caller VCF and extract FORMAT fields at target positions.
pub fn parse_caller_vcf(
    path: &Path,
    target_positions: &HashSet<(String, i64, String, String)>,
    sample_suffix: &str,
    caller_kind: CallerKind,
) -> Result<CallerResults, String> {
    if target_positions.is_empty() {
        return Ok(CallerResults::default());
    }

    let file = std::fs::File::open(path)
        .map_err(|e| format!("Cannot open {}: {}", path.display(), e))?;
    let buf_reader = BufReader::new(bgzf::Reader::new(file));
    let mut reader = vcf::io::Reader::new(buf_reader);
    let header = reader.read_header().map_err(|e| format!("Header error: {}", e))?;

    let sample_names: Vec<String> = header.sample_names().iter().map(|s| s.to_string()).collect();
    let sample_idx = sample_names.iter()
        .position(|name| name.ends_with(sample_suffix) || name == sample_suffix)
        .ok_or_else(|| format!("Sample suffix '{}' not found", sample_suffix))?;

    let mut results = CallerResults::default();
    let mut remaining = target_positions.clone();
    let has_gt_ad = caller_kind.has_gt_ad();
    let has_strand = caller_kind.has_strand_fields();
    let is_strelka = matches!(caller_kind, CallerKind::Strelka);
    let vaf_field = match caller_kind {
        CallerKind::Mutect2 => Some("AF"),
        CallerKind::DeepSomatic => Some("VAF"),
        CallerKind::Strelka => None,
    };

    for result in reader.records() {
        let record = match result { Ok(r) => r, Err(_) => continue };

        let chrom = record.reference_sequence_name().to_string();
        let pos = match record.variant_start() {
            Some(Ok(p)) => usize::from(p) as i64,
            _ => continue,
        };
        let ref_base = record.reference_bases().to_string();
        let alt_bases = record.alternate_bases();
        let alt = {
            let mut iter = alt_bases.iter();
            match iter.next() {
                Some(Ok(s)) => s.to_string(),
                _ => continue,
            }
        };

        let key = (chrom.clone(), pos, ref_base.clone(), alt.clone());
        if !remaining.contains(&key) { continue; }

        results.chroms.push(chrom);
        results.positions.push(pos);
        results.refs.push(ref_base);
        results.alts.push(alt);

        let samples = record.samples();
        let sample = samples.get_index(sample_idx);

        // DP
        results.dp.push(sample.as_ref().and_then(|s| get_int(s, &header, "DP")));

        // AD / GT
        if has_gt_ad {
            let (ad_ref, ad_alt) = match sample.as_ref().and_then(|s| s.get(&header, "AD")) {
                Some(Ok(Some(Value::Array(arr)))) => {
                    let vals: Vec<Option<i64>> = match arr {
                        vcf::variant::record::samples::series::value::Array::Integer(vals) =>
                            vals.iter().map(|r| r.ok().flatten().map(|v| v as i64)).collect(),
                        vcf::variant::record::samples::series::value::Array::Float(vals) =>
                            vals.iter().map(|r| r.ok().flatten().map(|v| v as i64)).collect(),
                        _ => Vec::new(),
                    };
                    let r = vals.first().and_then(|x| *x);
                    let a = vals.get(1).and_then(|x| *x);
                    (r, a)
                }
                _ => {
                    // Fallback: try as string
                    let ad_str = sample.as_ref().and_then(|s| get_str(s, &header, "AD"));
                    match ad_str {
                        Some(ref s) => {
                            let p: Vec<&str> = s.split(',').collect();
                            (p.first().and_then(|x| x.parse::<i64>().ok()),
                             p.get(1).and_then(|x| x.parse::<i64>().ok()))
                        }
                        None => (None, None),
                    }
                }
            };
            results.ad_ref.push(ad_ref);
            results.ad_alt.push(ad_alt);
            results.gt.push(sample.as_ref().and_then(|s| get_str(s, &header, "GT")));
        } else {
            results.ad_ref.push(None); results.ad_alt.push(None); results.gt.push(None);
        }

        // VAF caller
        if let Some(f) = vaf_field {
            results.vaf_caller.push(sample.as_ref().and_then(|s| get_float(s, &header, f)));
        } else { results.vaf_caller.push(None); }

        // Strelka
        if is_strelka {
            results.tar.push(sample.as_ref().and_then(|s| get_int_first(s, &header, "TAR")));
            results.tir.push(sample.as_ref().and_then(|s| get_int_first(s, &header, "TIR")));
            results.tor.push(sample.as_ref().and_then(|s| get_int_first(s, &header, "TOR")));
            results.au.push(sample.as_ref().and_then(|s| get_int_first(s, &header, "AU")));
            results.cu.push(sample.as_ref().and_then(|s| get_int_first(s, &header, "CU")));
            results.gu.push(sample.as_ref().and_then(|s| get_int_first(s, &header, "GU")));
            results.tu.push(sample.as_ref().and_then(|s| get_int_first(s, &header, "TU")));
        } else {
            for _ in 0..7 { results.tar.push(None); results.tir.push(None); results.tor.push(None);
                results.au.push(None); results.cu.push(None); results.gu.push(None); results.tu.push(None); }
        }

        // Mutect2 strand
        if has_strand {
            results.sb.push(sample.as_ref().and_then(|s| get_str(s, &header, "SB")));
            results.fad.push(sample.as_ref().and_then(|s| get_str(s, &header, "FAD")));
        } else { results.sb.push(None); results.fad.push(None); }

        remaining.remove(&key);
        if remaining.is_empty() { break; }
    }

    Ok(results)
}
