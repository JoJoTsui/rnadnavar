//! CxDy tiering engine — Rust implementation with parity to Python TieringEngine.
//!
//! Computes CxDy tiers from rescue VCF fields: FILTER, FILTERS_NORMALIZED,
//! GNOMAD_AF, COSMIC_CNT, REDI_EVIDENCE, N_DNA_CALLERS_SUPPORT, N_RNA_CALLERS_SUPPORT.

use std::collections::HashMap;

/// Tiering result per variant.
#[derive(Debug, Clone)]
pub struct TierResult {
    pub final_tier: String,
    pub caller_tier: String,
    pub database_tier: String,
    pub dna_caller_count: i64,
    pub rna_caller_count: i64,
    pub tier_quality: i64,
}

/// Quality scores from tier_config.py TIER_QUALITY_SCORES.
fn tier_quality(final_tier: &str) -> i64 {
    match final_tier {
        "C1D1" => 140, "C1D0" => 130,
        "C2D1" => 120, "C2D0" => 110,
        "C3D1" => 100, "C3D0" => 90,
        "C4D1" => 80,  "C4D0" => 70,
        "C5D1" => 60,  "C5D0" => 50,
        "C6D1" => 40,  "C6D0" => 30,
        "C7D1" => 20,  "C7D0" => 10,
        _ => 0,
    }
}

/// Compute caller tier C1-C7 from (dna_count, rna_count).
/// Matches CALLER_TIER_RULES from tier_config.py exactly.
fn compute_caller_tier(dna_count: i64, rna_count: i64) -> String {
    // Check in priority order C1 → C7. First match wins.
    if dna_count >= 2 && rna_count >= 2 { return "C1".into(); }
    if dna_count >= 2 && rna_count <= 1 { return "C2".into(); }
    if rna_count >= 2 && dna_count <= 1 { return "C3".into(); }
    if dna_count == 1 && rna_count == 1  { return "C4".into(); }
    if dna_count == 1 && rna_count == 0  { return "C5".into(); }
    if dna_count == 0 && rna_count == 1  { return "C6".into(); }
    "C7".into() // dna_count == 0 && rna_count == 0
}

/// Parse FILTERS_NORMALIZED string into per-caller fields.
/// Format: "DNA_strelka:Somatic|RNA_mutect2:Germline|..."
fn parse_filters_normalized(filters_str: &str) -> HashMap<String, String> {
    let mut result = HashMap::new();
    if filters_str.is_empty() { return result; }

    let delim = if filters_str.contains('|') { '|' } else { ';' };
    for entry in filters_str.split(delim) {
        let entry = entry.trim();
        if entry.is_empty() || !entry.contains(':') { continue; }

        let parts: Vec<&str> = entry.splitn(2, ':').collect();
        if parts.len() != 2 { continue; }

        let caller_spec = parts[0].trim();
        let category = parts[1].trim().to_string();

        // Parse "DNA_strelka" or "RNA_mutect2" format
        let lower = caller_spec.to_lowercase();
        let (modality, caller) = if let Some(c) = lower.strip_prefix("dna_") {
            ("DNA", c)
        } else if let Some(c) = lower.strip_prefix("rna_") {
            ("RNA", c)
        } else {
            continue;
        };

        // Normalize caller name
        let caller_normalized = match caller {
            "strelka" => "Strelka",
            "deepsomatic" => "DeepSomatic",
            "mutect2" => "Mutect2",
            _ => continue,
        };

        let field_name = format!("FILTER_NORMALIZED_{}_{}_TUMOR", caller_normalized, modality);
        result.insert(field_name, category);
    }
    result
}

/// Count category-concordant DNA/RNA callers.
/// Matches Python category_matcher.count_concordant_callers() exactly.
fn count_concordant_callers(
    final_filter: &str,
    filter_normalized_fields: &HashMap<String, String>,
) -> (i64, i64) {
    let mut dna_count: i64 = 0;
    let mut rna_count: i64 = 0;
    let mut counted_dna = std::collections::HashSet::new();
    let mut counted_rna = std::collections::HashSet::new();

    for (field_name, field_value) in filter_normalized_fields {
        if field_value.is_empty() || field_value != final_filter { continue; }

        // Parse FILTER_NORMALIZED_<Caller>_<Modality>_TUMOR
        if !field_name.starts_with("FILTER_NORMALIZED_") { continue; }
        let rest = &field_name["FILTER_NORMALIZED_".len()..];
        let parts: Vec<&str> = rest.rsplitn(3, '_').collect();
        // parts: [TUMOR, RNA|DNA, Caller] or similar
        if parts.len() < 3 { continue; }

        // Last segment should be "TUMOR"
        if parts[0] != "TUMOR" { continue; }
        let modality = parts[1];
        let caller = parts[2..].join("_"); // Handle multi-word caller names

        match modality {
            "DNA" => {
                if !counted_dna.contains(&caller) {
                    dna_count += 1;
                    counted_dna.insert(caller.clone());
                }
            }
            "RNA" => {
                if !counted_rna.contains(&caller) {
                    rna_count += 1;
                    counted_rna.insert(caller.clone());
                }
            }
            _ => {}
        }
    }
    (dna_count, rna_count)
}

/// Check database evidence support.
/// Matches Python database_checker logic.
fn has_database_support(
    gnomad_af: Option<f64>,
    cosmic_cnt: Option<i64>,
    redi_evidence: Option<&str>,
    info_dict: &HashMap<String, String>,
) -> bool {
    // gnomAD: AF > 0.001
    if let Some(af) = gnomad_af {
        if af > 0.001 { return true; }
    } else if let Some(af_str) = info_dict.get("GNOMAD_AF") {
        if let Ok(af) = af_str.parse::<f64>() {
            if af > 0.001 { return true; }
        }
    }

    // COSMIC: CNT > 0 or ID present
    if let Some(cnt) = cosmic_cnt {
        if cnt > 0 { return true; }
    } else if let Some(cs_str) = info_dict.get("COSMIC_CNT") {
        if let Ok(cnt) = cs_str.parse::<i64>() {
            if cnt > 0 { return true; }
        }
    }
    if let Some(cs_id) = info_dict.get("COSMIC_ID") {
        if !cs_id.is_empty() && cs_id != "." { return true; }
    }

    // REDIportal: check REDI_EVIDENCE
    if let Some(ev) = redi_evidence {
        match ev {
            "HIGH" | "MEDIUM" | "LOW" => return true,
            _ => {}
        }
    }

    false
}

/// Compute tier for a single variant.
pub fn compute_tier(
    final_filter: &str,
    filters_normalized: &str,
    gnomad_af: Option<f64>,
    cosmic_cnt: Option<i64>,
    redi_evidence: Option<&str>,
    dna_support: Option<i64>,
    rna_support: Option<i64>,
    info_dict: &HashMap<String, String>,
) -> TierResult {
    let final_filter = if final_filter.is_empty() || final_filter == "." { "PASS" } else { final_filter };

    // Parse FILTERS_NORMALIZED into per-caller fields
    let mut filter_normalized = parse_filters_normalized(filters_normalized);

    // Also check for individual FILTER_NORMALIZED_* columns in info_dict
    for (key, value) in info_dict {
        if key.starts_with("FILTER_NORMALIZED_") {
            filter_normalized.entry(key.clone()).or_insert_with(|| value.clone());
        }
    }

    // Count concordant callers
    let (mut dna_count, mut rna_count) = count_concordant_callers(final_filter, &filter_normalized);

    // If no concordant callers found via FILTERS_NORMALIZED, use fallback counts
    if dna_count == 0 && rna_count == 0 {
        dna_count = dna_support.unwrap_or(0);
        rna_count = rna_support.unwrap_or(0);
    }

    // Database tier
    let db_support = has_database_support(gnomad_af, cosmic_cnt, redi_evidence, info_dict);
    let database_tier = if db_support { "D1" } else { "D0" };

    // Caller tier
    let caller_tier = compute_caller_tier(dna_count, rna_count);

    // Final tier
    let final_tier = format!("{}{}", caller_tier, database_tier);
    let quality = tier_quality(&final_tier);

    TierResult {
        final_tier,
        caller_tier,
        database_tier: database_tier.to_string(),
        dna_caller_count: dna_count,
        rna_caller_count: rna_count,
        tier_quality: quality,
    }
}

/// Batch compute tiers for all variants.
pub fn compute_tiers_batch(
    filters: &[String],
    filters_normalized: &[String],
    gnomad_af: &[Option<f64>],
    cosmic_cnt: &[Option<i64>],
    redi_evidence: &[Option<String>],
    dna_support: &[Option<i64>],
    rna_support: &[Option<i64>],
    info_keys: &[Vec<String>],   // per-variant: list of additional INFO key-value pairs
    info_vals: &[Vec<String>],   // per-variant: matching values
) -> Vec<TierResult> {
    let n = filters.len();
    let mut results = Vec::with_capacity(n);

    for i in 0..n {
        // Build info_dict from keys/vals
        let mut info_dict = HashMap::new();
        if i < info_keys.len() && i < info_vals.len() {
            for (k, v) in info_keys[i].iter().zip(info_vals[i].iter()) {
                info_dict.insert(k.clone(), v.clone());
            }
        }

        let result = compute_tier(
            &filters[i],
            &filters_normalized.get(i).map(|s| s.as_str()).unwrap_or(""),
            gnomad_af.get(i).and_then(|v| *v),
            cosmic_cnt.get(i).and_then(|v| *v),
            redi_evidence.get(i).and_then(|o| o.as_deref()),
            dna_support.get(i).and_then(|v| *v),
            rna_support.get(i).and_then(|v| *v),
            &info_dict,
        );
        results.push(result);
    }

    results
}
