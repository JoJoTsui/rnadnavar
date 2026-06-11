//! VCF parsing using noodles-vcf 0.88.

use std::collections::HashMap;
use std::path::Path;

use std::fs::File;
use std::io::{BufRead, BufReader};

use noodles_bgzf as bgzf;
use noodles_vcf as vcf;
use noodles_vcf::variant::record::{AlternateBases as _, Filters as _};

/// A single variant record extracted from a rescue VCF (legacy row-oriented).
#[derive(Debug, Clone)]
pub struct RescueRecord {
    pub chrom: String,
    pub pos: i64,
    pub ref_base: String,
    pub alt: String,
    pub filter: String,
    pub info: Vec<(String, String)>,
}

/// Column-oriented rescue VCF parse results.
///
/// Each column is a separate Vec — no per-record Python dicts are created.
/// INFO columns use HashMap keyed by INFO field name. `info_keys` preserves
/// the VCF header declaration order for deterministic column output.
#[derive(Debug, Clone)]
pub struct RescueColumns {
    pub chrom: Vec<String>,
    pub pos: Vec<i64>,
    pub ref_base: Vec<String>,
    pub alt: Vec<String>,
    pub filter: Vec<String>,
    pub info: HashMap<String, Vec<Option<String>>>,
    pub info_keys: Vec<String>,
    pub variant_type: Vec<String>,
    pub ti_tv: Vec<Option<bool>>,
}

/// Classify a variant as SNV, INS, DEL, or MNV based on REF/ALT lengths.
fn derive_variant_type(ref_len: usize, alt_len: usize) -> &'static str {
    if ref_len == 1 && alt_len == 1 {
        "SNV"
    } else if alt_len > ref_len {
        "INS"
    } else if ref_len > alt_len {
        "DEL"
    } else {
        "MNV"
    }
}

/// Determine if a single-nucleotide substitution is a transition.
/// Returns None for non-SNV variants.
fn is_transition(ref_base: &str, alt_base: &str) -> Option<bool> {
    if ref_base.len() != 1 || alt_base.len() != 1 {
        return None;
    }
    let r = ref_base.to_uppercase();
    let a = alt_base.to_uppercase();
    match (r.as_str(), a.as_str()) {
        ("A", "G") | ("G", "A") | ("C", "T") | ("T", "C") => Some(true),
        _ => Some(false),
    }
}

/// Build a VCF header from raw text, deduplicating FILTER/INFO lines.
/// Rescue VCFs have duplicate FILTER IDs (e.g., "Somatic" appears twice)
/// that noodles-vcf rejects. This function reads the header text manually
/// and builds a clean header.
fn build_header(path: &Path) -> Result<(vcf::Header, Vec<String>), Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let mut reader = bgzf::Reader::new(file);
    let mut header_lines = String::new();
    let mut buf = String::new();

    // Read lines until the "#CHROM" header line
    loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 { break; }
        header_lines.push_str(&buf);
        if buf.starts_with("#CHROM") { break; }
    }

    // Deduplicate FILTER and INFO lines: keep first occurrence
    let mut seen_filters = std::collections::HashSet::new();
    let mut seen_infos = std::collections::HashSet::new();
    let mut cleaned = String::new();
    let mut info_keys = Vec::new();

    for line in header_lines.lines() {
        if line.starts_with("##FILTER=<ID=") {
            if let Some(id) = extract_meta_id(line, "FILTER") {
                if !seen_filters.insert(id.clone()) {
                    continue; // skip duplicate
                }
            }
        } else if line.starts_with("##INFO=<ID=") {
            if let Some(id) = extract_meta_id(line, "INFO") {
                if !seen_infos.insert(id.clone()) {
                    continue; // skip duplicate
                }
                info_keys.push(id);
            }
        }
        cleaned.push_str(line);
        cleaned.push('\n');
    }

    // Read back the clean header into a BufReader
    let cursor = std::io::Cursor::new(cleaned);
    let mut header_reader = vcf::io::Reader::new(cursor);
    let header = header_reader.read_header()?;

    Ok((header, info_keys))
}

/// Extract the ID from a ##FILTER=<ID=X,...> or ##INFO=<ID=X,...> line.
fn extract_meta_id(line: &str, prefix: &str) -> Option<String> {
    let marker = format!("##{}=", prefix);
    line.strip_prefix(&marker)
        .and_then(|rest| rest.strip_prefix('<'))
        .and_then(|rest| {
            rest.split(',')
                .find(|part| part.starts_with("ID="))
                .and_then(|id_part| id_part.strip_prefix("ID="))
                .map(|id| id.trim().to_string())
        })
}

/// Parse a rescue VCF file into column-oriented data.
///
/// Returns a `RescueColumns` struct where each field is a separate column vector.
/// INFO fields are stored in a `HashMap<String, Vec<Option<String>>>` with keys
/// ordered by `info_keys` for deterministic column output.
/// `variant_type` and `ti_tv` are computed during the parse pass.
pub fn parse_rescue_columns(path: &Path) -> Result<RescueColumns, Box<dyn std::error::Error>> {
    let (header, info_keys) = build_header(path)?;

    // Initialize column vectors
    let mut chrom: Vec<String> = Vec::new();
    let mut pos: Vec<i64> = Vec::new();
    let mut ref_base: Vec<String> = Vec::new();
    let mut alt: Vec<String> = Vec::new();
    let mut filter: Vec<String> = Vec::new();
    let mut variant_type: Vec<String> = Vec::new();
    let mut ti_tv: Vec<Option<bool>> = Vec::new();

    // Initialize INFO columns: one Vec per declared INFO key
    let mut info_columns: HashMap<String, Vec<Option<String>>> = HashMap::with_capacity(info_keys.len());
    for key in &info_keys {
        info_columns.insert(key.clone(), Vec::new());
    }

    // Open the VCF for record reading
    let file = File::open(path)?;
    let mut buf_reader = BufReader::new(bgzf::Reader::new(file));
    let mut line = String::new();
    loop {
        line.clear();
        let n = buf_reader.read_line(&mut line)?;
        if n == 0 || line.starts_with("#CHROM") { break; }
    }
    let mut reader = vcf::io::Reader::new(buf_reader);

    for result in reader.records() {
        let record = result?;

        // CHROM and POS
        let c = record.reference_sequence_name().to_string();
        let p = match record.variant_start() {
            Some(Ok(p)) => usize::from(p) as i64,
            _ => 0,
        };

        // REF and ALT
        let r = record.reference_bases().to_string();
        let alt_bases = record.alternate_bases();
        let a = {
            let mut iter = alt_bases.iter();
            match iter.next() {
                Some(Ok(s)) => s.to_string(),
                _ => ".".to_string(),
            }
        };

        // FILTER
        let record_filters = record.filters();
        let f = {
            let filters_iter = record_filters.iter(&header);
            let ids: Vec<&str> = filters_iter.filter_map(|r| r.ok()).collect();
            if ids.is_empty() { "PASS".to_string() } else { ids.join(";") }
        };

        // Derived columns
        let vt = derive_variant_type(r.len(), a.len());
        let tt = is_transition(&r, &a);

        chrom.push(c);
        pos.push(p);
        ref_base.push(r);
        alt.push(a);
        filter.push(f);
        variant_type.push(vt.to_string());
        ti_tv.push(tt);

        // INFO columns: push Some(value) or None for each declared key
        for key in &info_keys {
            let value = match record.info().get(&header, key.as_str()) {
                Some(Ok(Some(val))) => Some(format_info_value(&val)),
                _ => None,
            };
            info_columns.get_mut(key).unwrap().push(value);
        }
    }

    Ok(RescueColumns {
        chrom,
        pos,
        ref_base,
        alt,
        filter,
        info: info_columns,
        info_keys,
        variant_type,
        ti_tv,
    })
}

/// Parse a rescue VCF file and return all variant records with INFO fields (legacy row-oriented).
pub fn parse_rescue_vcf(path: &Path) -> Result<Vec<RescueRecord>, Box<dyn std::error::Error>> {
    // Build a clean header by deduplicating meta lines
    let (header, info_keys) = build_header(path)?;

    // Open the VCF for record reading — skip the raw header (has duplicates)
    let file = File::open(path)?;
    let mut buf_reader = BufReader::new(bgzf::Reader::new(file));
    // Read and discard header lines until #CHROM
    let mut line = String::new();
    loop {
        line.clear();
        let n = buf_reader.read_line(&mut line)?;
        if n == 0 || line.starts_with("#CHROM") { break; }
    }
    // Now buf_reader is positioned after #CHROM line, rest is records
    let mut reader = vcf::io::Reader::new(buf_reader);

    let mut records = Vec::new();

    for result in reader.records() {
        let record = result?;

        // CHROM and POS
        let chrom = record.reference_sequence_name().to_string();
        let pos = match record.variant_start() {
            Some(Ok(p)) => usize::from(p) as i64,
            _ => 0,
        };

        // REF and ALT
        let ref_base = record.reference_bases().to_string();
        let alt_bases = record.alternate_bases();
        let alt = {
            let mut iter = alt_bases.iter();
            match iter.next() {
                Some(Ok(s)) => s.to_string(),
                _ => ".".to_string(),
            }
        };

        // FILTER
        let filters = record.filters();
        let filter = {
            let filters_iter = filters.iter(&header);
            let ids: Vec<&str> = filters_iter.filter_map(|r| r.ok()).collect();
            if ids.is_empty() { "PASS".to_string() } else { ids.join(";") }
        };

        // INFO — parse lazily, collect as strings
        let info: Vec<(String, String)> = info_keys
            .iter()
            .filter_map(|key| {
                match record.info().get(&header, key.as_str()) {
                    Some(Ok(Some(value))) => {
                        let s = format_info_value(&value);
                        Some((key.clone(), s))
                    }
                    _ => None,
                }
            })
            .collect();

        records.push(RescueRecord { chrom, pos, ref_base, alt, filter, info });
    }

    Ok(records)
}

/// Format an INFO field Value to a string.
fn format_info_value(value: &vcf::variant::record::info::field::Value<'_>) -> String {
    use vcf::variant::record::info::field::Value;
    use vcf::variant::record::info::field::value::Array;

    match value {
        Value::Integer(n) => n.to_string(),
        Value::Float(f) => f.to_string(),
        Value::Flag => "true".to_string(),
        Value::String(s) => s.to_string(),
        Value::Character(c) => c.to_string(),
        Value::Array(arr) => {
            let strings: Vec<String> = match arr {
                Array::Integer(vals) => vals.iter().filter_map(|r| r.ok()).flatten().map(|v| v.to_string()).collect(),
                Array::Float(vals) => vals.iter().filter_map(|r| r.ok()).flatten().map(|v| v.to_string()).collect(),
                Array::String(vals) => vals.iter().filter_map(|r| r.ok()).flatten().map(|v| v.to_string()).collect(),
                Array::Character(vals) => vals.iter().filter_map(|r| r.ok()).flatten().map(|v| v.to_string()).collect(),
            };
            strings.join(",")
        },
    }
}
