mod bam;
mod caller;
mod pileup;
mod tier;
mod vcf;

use std::collections::HashSet;
use std::path::PathBuf;

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

/// Parse a rescue VCF and return a list of dicts (CHROM, POS, REF, ALT, FILTER, INFO...).
///
/// Rust parsing is done with the GIL released so multiple samples' VCFs can
/// be parsed concurrently. Python object construction re-acquires the GIL.
#[pyfunction]
fn parse_rescue(py: Python<'_>, path: String) -> PyResult<Bound<'_, PyList>> {
    let path_buf = PathBuf::from(&path);
    // detach() releases the GIL during Rust VCF parsing (the CPU-intensive part).
    // Must convert errors to String because Box<dyn Error> is not Ungil.
    let records = py.detach(|| {
        vcf::parse_rescue_vcf(&path_buf)
            .map_err(|e| e.to_string())
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

    // Re-acquire GIL for Python object construction
    let list = PyList::empty(py);
    for rec in &records {
        let d = PyDict::new(py);
        d.set_item("CHROM", &rec.chrom)?;
        d.set_item("POS", rec.pos)?;
        d.set_item("REF", &rec.ref_base)?;
        d.set_item("ALT", &rec.alt)?;
        d.set_item("FILTER", &rec.filter)?;
        for (key, value) in &rec.info {
            d.set_item(key.as_str(), value.as_str())?;
        }
        list.append(d)?;
    }

    Ok(list)
}

/// Parse a rescue VCF into column-oriented data.
///
/// Returns a dict of column_name → list_of_values instead of a list of per-record
/// dicts. Eliminates the intermediate 7M-PyDict construction for large samples.
///
/// Rust parsing is done with the GIL released so multiple samples' VCFs can
/// be parsed concurrently. Python object construction re-acquires the GIL.
#[pyfunction]
fn parse_rescue_columns(py: Python<'_>, path: String) -> PyResult<Bound<'_, PyDict>> {
    let path_buf = PathBuf::from(&path);
    let columns = py.detach(|| {
        vcf::parse_rescue_columns(&path_buf)
            .map_err(|e| e.to_string())
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

    // Re-acquire GIL for Python object construction
    let d = PyDict::new(py);

    fn push_str_list(py: Python<'_>, dict: &Bound<'_, PyDict>, key: &str, vals: &[String]) -> PyResult<()> {
        let list = PyList::empty(py);
        for v in vals { list.append(v.as_str())?; }
        dict.set_item(key, list)
    }
    fn push_i64_list(py: Python<'_>, dict: &Bound<'_, PyDict>, key: &str, vals: &[i64]) -> PyResult<()> {
        let list = PyList::empty(py);
        for v in vals { list.append(*v)?; }
        dict.set_item(key, list)
    }
    fn push_opt_str_list(py: Python<'_>, dict: &Bound<'_, PyDict>, key: &str, vals: &[Option<String>]) -> PyResult<()> {
        let list = PyList::empty(py);
        for v in vals { list.append(v.as_deref())?; }
        dict.set_item(key, list)
    }
    fn push_opt_bool_list(py: Python<'_>, dict: &Bound<'_, PyDict>, key: &str, vals: &[Option<bool>]) -> PyResult<()> {
        let list = PyList::empty(py);
        for v in vals { list.append(*v)?; }
        dict.set_item(key, list)
    }

    push_str_list(py, &d, "CHROM", &columns.chrom)?;
    push_i64_list(py, &d, "POS", &columns.pos)?;
    push_str_list(py, &d, "REF", &columns.ref_base)?;
    push_str_list(py, &d, "ALT", &columns.alt)?;
    push_str_list(py, &d, "FILTER", &columns.filter)?;
    push_str_list(py, &d, "variant_type", &columns.variant_type)?;
    push_opt_bool_list(py, &d, "ti_tv", &columns.ti_tv)?;

    // Push each INFO column in header order
    for key in &columns.info_keys {
        if let Some(vals) = columns.info.get(key) {
            push_opt_str_list(py, &d, key, vals)?;
        }
    }

    Ok(d)
}

/// Compute whole-genome BAM statistics.
///
/// Releases the GIL during the BAM scan so multiple threads can process
/// different BAM files concurrently.
#[pyfunction]
fn bam_stats(py: Python<'_>, path: String, max_reads: u64) -> PyResult<Bound<'_, PyDict>> {
    let path_buf = PathBuf::from(&path);
    let stats = py.detach(|| {
        bam::whole_genome_stats(&path_buf, max_reads)
            .map_err(|e| e.to_string())
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

    let d = PyDict::new(py);
    d.set_item("total_reads", stats.total_reads)?;
    d.set_item("mapped_reads", stats.mapped_reads)?;
    d.set_item("mapping_rate", stats.mapping_rate)?;
    d.set_item("mean_coverage", stats.mean_coverage)?;
    d.set_item("mean_insert_size", stats.mean_insert_size)?;
    d.set_item("mean_mapq", stats.mean_mapq)?;
    d.set_item("duplication_rate_pct", stats.duplication_rate_pct)?;
    d.set_item("properly_paired_pct", stats.properly_paired_pct)?;
    d.set_item("insert_size_stddev", stats.insert_size_stddev)?;
    Ok(d)
}

/// Compute BAM statistics with BED-guided on-target coverage.
///
/// When BED regions are provided, only bases overlapping BED intervals
/// are counted toward mean_coverage (WES mode). Other metrics (reads,
/// mapping rate, insert size, MAPQ) remain whole-genome.
#[pyfunction]
fn bam_stats_bed(
    py: Python<'_>,
    path: String,
    max_reads: u64,
    bed_chroms: Vec<String>,
    bed_starts: Vec<i64>,
    bed_ends: Vec<i64>,
) -> PyResult<Bound<'_, PyDict>> {
    let path_buf = PathBuf::from(&path);
    let stats = py.detach(|| {
        bam::whole_genome_stats_bed(&path_buf, max_reads, &bed_chroms, &bed_starts, &bed_ends)
            .map_err(|e| e.to_string())
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

    let d = PyDict::new(py);
    d.set_item("total_reads", stats.total_reads)?;
    d.set_item("mapped_reads", stats.mapped_reads)?;
    d.set_item("mapping_rate", stats.mapping_rate)?;
    d.set_item("mean_coverage", stats.mean_coverage)?;
    d.set_item("mean_insert_size", stats.mean_insert_size)?;
    d.set_item("mean_mapq", stats.mean_mapq)?;
    d.set_item("duplication_rate_pct", stats.duplication_rate_pct)?;
    d.set_item("properly_paired_pct", stats.properly_paired_pct)?;
    d.set_item("insert_size_stddev", stats.insert_size_stddev)?;
    Ok(d)
}

/// Parse a normalized caller VCF and return FORMAT fields at target positions.
///
/// Releases the GIL during parsing. Returns a dict of column_name → list_of_values.
#[pyfunction]
fn parse_caller_vcf(
    py: Python<'_>,
    path: String,
    target_chroms: Vec<String>,
    target_positions: Vec<i64>,
    target_refs: Vec<String>,
    target_alts: Vec<String>,
    sample_suffix: String,
    caller_name: String,
) -> PyResult<Bound<'_, PyDict>> {
    let path_buf = PathBuf::from(&path);
    let targets: HashSet<(String, i64, String, String)> = target_chroms.into_iter()
        .zip(target_positions.into_iter())
        .zip(target_refs.into_iter())
        .zip(target_alts.into_iter())
        .map(|(((c, p), r), a)| (c, p, r, a))
        .collect();
    let kind = caller::CallerKind::from_name(&caller_name);

    let results = py.detach(|| {
        caller::parse_caller_vcf(&path_buf, targets, &sample_suffix, kind)
            .map_err(|e| e.to_string())
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

    let d = PyDict::new(py);

    fn push_str_list(py: Python<'_>, dict: &Bound<'_, PyDict>, key: &str, vals: &[String]) -> PyResult<()> {
        let list = PyList::empty(py);
        for v in vals { list.append(v.as_str())?; }
        dict.set_item(key, list)
    }
    fn push_i64_list(py: Python<'_>, dict: &Bound<'_, PyDict>, key: &str, vals: &[i64]) -> PyResult<()> {
        let list = PyList::empty(py);
        for v in vals { list.append(*v)?; }
        dict.set_item(key, list)
    }
    fn push_opt_i64_list(py: Python<'_>, dict: &Bound<'_, PyDict>, key: &str, vals: &[Option<i64>]) -> PyResult<()> {
        let list = PyList::empty(py);
        for v in vals { list.append(*v)?; }
        dict.set_item(key, list)
    }
    fn push_opt_f64_list(py: Python<'_>, dict: &Bound<'_, PyDict>, key: &str, vals: &[Option<f64>]) -> PyResult<()> {
        let list = PyList::empty(py);
        for v in vals { list.append(*v)?; }
        dict.set_item(key, list)
    }
    fn push_opt_str_list(py: Python<'_>, dict: &Bound<'_, PyDict>, key: &str, vals: &[Option<String>]) -> PyResult<()> {
        let list = PyList::empty(py);
        for v in vals { list.append(v.as_deref())?; }
        dict.set_item(key, list)
    }

    push_str_list(py, &d, "CHROM", &results.chroms)?;
    push_i64_list(py, &d, "POS", &results.positions)?;
    push_str_list(py, &d, "REF", &results.refs)?;
    push_str_list(py, &d, "ALT", &results.alts)?;
    push_opt_i64_list(py, &d, "DP", &results.dp)?;
    push_opt_i64_list(py, &d, "AD_REF", &results.ad_ref)?;
    push_opt_i64_list(py, &d, "AD_ALT", &results.ad_alt)?;
    push_opt_str_list(py, &d, "GT", &results.gt)?;
    push_opt_f64_list(py, &d, "VAF_CALLER", &results.vaf_caller)?;
    push_opt_i64_list(py, &d, "TAR", &results.tar)?;
    push_opt_i64_list(py, &d, "TIR", &results.tir)?;
    push_opt_i64_list(py, &d, "TOR", &results.tor)?;
    push_opt_i64_list(py, &d, "AU", &results.au)?;
    push_opt_i64_list(py, &d, "CU", &results.cu)?;
    push_opt_i64_list(py, &d, "GU", &results.gu)?;
    push_opt_i64_list(py, &d, "TU", &results.tu)?;
    push_opt_str_list(py, &d, "SB", &results.sb)?;
    push_opt_str_list(py, &d, "FAD", &results.fad)?;
    Ok(d)
}

/// Compute CxDy tiers for all variants.
///
/// Releases the GIL during computation. Returns a dict of column_name → list.
#[pyfunction]
fn compute_tiers(
    py: Python<'_>,
    filters: Vec<String>,
    filters_normalized: Vec<String>,
    gnomad_af: Vec<Option<f64>>,
    cosmic_cnt: Vec<Option<i64>>,
    redi_evidence: Vec<Option<String>>,
    dna_support: Vec<Option<i64>>,
    rna_support: Vec<Option<i64>>,
) -> PyResult<Bound<'_, PyDict>> {
    let results = py.detach(|| {
        Ok::<_, String>(tier::compute_tiers_batch(
            &filters,
            &filters_normalized,
            &gnomad_af,
            &cosmic_cnt,
            &redi_evidence,
            &dna_support,
            &rna_support,
            &vec![Vec::new(); filters.len()],  // info_keys: empty
            &vec![Vec::new(); filters.len()],  // info_vals: empty
        ))
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

    let n = results.len();
    let d = PyDict::new(py);

    let mut final_tiers = Vec::with_capacity(n);
    let mut caller_tiers = Vec::with_capacity(n);
    let mut database_tiers = Vec::with_capacity(n);
    let mut dna_counts = Vec::with_capacity(n);
    let mut rna_counts = Vec::with_capacity(n);
    let mut qualities = Vec::with_capacity(n);

    for r in &results {
        final_tiers.push(r.final_tier.clone());
        caller_tiers.push(r.caller_tier.clone());
        database_tiers.push(r.database_tier.clone());
        dna_counts.push(r.dna_caller_count);
        rna_counts.push(r.rna_caller_count);
        qualities.push(r.tier_quality);
    }

    {
        let list = PyList::empty(py);
        for v in &final_tiers { list.append(v.as_str())?; }
        d.set_item("final_tier", list)?;
    }
    {
        let list = PyList::empty(py);
        for v in &caller_tiers { list.append(v.as_str())?; }
        d.set_item("caller_tier", list)?;
    }
    {
        let list = PyList::empty(py);
        for v in &database_tiers { list.append(v.as_str())?; }
        d.set_item("database_tier", list)?;
    }
    {
        let list = PyList::empty(py);
        for v in &dna_counts { list.append(*v)?; }
        d.set_item("dna_caller_count", list)?;
    }
    {
        let list = PyList::empty(py);
        for v in &rna_counts { list.append(*v)?; }
        d.set_item("rna_caller_count", list)?;
    }
    {
        let list = PyList::empty(py);
        for v in &qualities { list.append(*v)?; }
        d.set_item("tier_quality", list)?;
    }
    Ok(d)
}

/// Perform BAM pileup at variant positions.
///
/// Releases the GIL during pileup.
#[pyfunction]
fn pileup_variants(
    py: Python<'_>,
    bam_path: String,
    chroms: Vec<String>,
    positions: Vec<i64>,
    ref_bases: Vec<String>,
    alt_bases: Vec<String>,
) -> PyResult<Bound<'_, PyDict>> {
    let path_buf = PathBuf::from(&bam_path);
    let results = py.detach(|| {
        pileup::pileup_variants(&path_buf, &chroms, &positions, &ref_bases, &alt_bases)
            .map_err(|e| e.to_string())
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

    pileup_results_to_pydict(py, &results)
}

/// Perform combined BAM pileup across multiple BAMs in a single FFI call.
///
/// Opens all BAMs once, shares position grouping and window iteration,
/// and uses binary search to match reads to positions within each region.
/// When bed_regions is provided, uses BED intervals instead of 1Mb windows.
///
/// Returns a dict mapping bam_label → dict of pileup columns.
#[pyfunction]
fn pileup_variants_multi(
    py: Python<'_>,
    bam_paths: Vec<String>,
    bam_labels: Vec<String>,
    chroms: Vec<String>,
    positions: Vec<i64>,
    ref_bases: Vec<String>,
    alt_bases: Vec<String>,
    bed_chroms: Vec<String>,
    bed_starts: Vec<i64>,
    bed_ends: Vec<i64>,
) -> PyResult<Bound<'_, PyDict>> {
    let bed_regions: Option<Vec<(String, i64, i64)>> = if bed_chroms.is_empty() {
        None
    } else {
        let n = bed_chroms.len();
        Some((0..n).map(|i| (bed_chroms[i].clone(), bed_starts[i], bed_ends[i])).collect())
    };

    let multi_results = py.detach(|| {
        pileup::pileup_variants_multi(
            &bam_paths, &bam_labels,
            &chroms, &positions, &ref_bases, &alt_bases,
            bed_regions.as_deref(),
        )
        .map_err(|e| e.to_string())
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

    let d = PyDict::new(py);
    for (label, results) in &multi_results {
        d.set_item(label.as_str(), pileup_results_to_pydict(py, results)?)?;
    }
    Ok(d)
}

/// Compute per-base coverage depth bins across BED regions.
///
/// Returns a dict mapping "cov_Nx_pct" → percentage of on-target bases with
/// depth >= N, or None when no BED regions are provided.
#[pyfunction]
fn coverage_bins(py: Python<'_>, bam_path: String, bed_regions: Vec<(String, u32, u32)>) -> PyResult<Option<Bound<'_, PyDict>>> {
    let path_buf = PathBuf::from(&bam_path);
    let result = py.detach(|| {
        bam::coverage_bins(&path_buf, &bed_regions)
            .map_err(|e| e.to_string())
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

    match result {
        Some(map) => {
            let dict = PyDict::new(py);
            for (k, v) in map { dict.set_item(k, v)?; }
            Ok(Some(dict))
        }
        None => Ok(None),
    }
}

/// Compute whole-genome coverage depth bins for all reference sequences.
///
/// Returns a dict mapping "cov_Nx_pct" → percentage of genome bases with
/// depth >= N, or None when no BAI index is available.
#[pyfunction]
fn wg_coverage_bins(py: Python<'_>, bam_path: String) -> PyResult<Option<Bound<'_, PyDict>>> {
    let path_buf = PathBuf::from(&bam_path);
    let result = py.detach(|| {
        bam::wg_coverage_bins(&path_buf)
            .map_err(|e| e.to_string())
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e))?;

    match result {
        Some(map) => {
            let dict = PyDict::new(py);
            for (k, v) in map { dict.set_item(k, v)?; }
            Ok(Some(dict))
        }
        None => Ok(None),
    }
}

/// Convert a Vec<PileupResult> to a Python dict of column lists.
fn pileup_results_to_pydict<'a>(py: Python<'a>, results: &[pileup::PileupResult]) -> PyResult<Bound<'a, PyDict>> {
    let d = PyDict::new(py);
    let dp: Vec<Option<i64>> = results.iter().map(|r| r.dp).collect();
    let ref_dp: Vec<Option<i64>> = results.iter().map(|r| r.ref_dp).collect();
    let alt_dp: Vec<Option<i64>> = results.iter().map(|r| r.alt_dp).collect();
    let f1r2_ref: Vec<Option<i64>> = results.iter().map(|r| r.f1r2_ref).collect();
    let f2r1_ref: Vec<Option<i64>> = results.iter().map(|r| r.f2r1_ref).collect();
    let f1r2_alt: Vec<Option<i64>> = results.iter().map(|r| r.f1r2_alt).collect();
    let f2r1_alt: Vec<Option<i64>> = results.iter().map(|r| r.f2r1_alt).collect();
    let mean_bq: Vec<Option<f64>> = results.iter().map(|r| r.mean_bq).collect();
    let mean_mq: Vec<Option<f64>> = results.iter().map(|r| r.mean_mq).collect();

    {
        let lst = PyList::empty(py); for v in &dp { lst.append(*v)?; } d.set_item("DP", lst)?;
    }
    {
        let lst = PyList::empty(py); for v in &ref_dp { lst.append(*v)?; } d.set_item("REF_DP", lst)?;
    }
    {
        let lst = PyList::empty(py); for v in &alt_dp { lst.append(*v)?; } d.set_item("ALT_DP", lst)?;
    }
    {
        let lst = PyList::empty(py); for v in &f1r2_ref { lst.append(*v)?; } d.set_item("F1R2_ref", lst)?;
    }
    {
        let lst = PyList::empty(py); for v in &f2r1_ref { lst.append(*v)?; } d.set_item("F2R1_ref", lst)?;
    }
    {
        let lst = PyList::empty(py); for v in &f1r2_alt { lst.append(*v)?; } d.set_item("F1R2_alt", lst)?;
    }
    {
        let lst = PyList::empty(py); for v in &f2r1_alt { lst.append(*v)?; } d.set_item("F2R1_alt", lst)?;
    }
    {
        let lst = PyList::empty(py); for v in &mean_bq { lst.append(*v)?; } d.set_item("mean_BQ", lst)?;
    }
    {
        let lst = PyList::empty(py); for v in &mean_mq { lst.append(*v)?; } d.set_item("mean_MQ", lst)?;
    }
    Ok(d)
}

/// stats_core — Rust-accelerated VCF/BAM parsing for seq2neo variant statistics.
#[pymodule]
fn stats_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_rescue, m)?)?;
    m.add_function(wrap_pyfunction!(parse_rescue_columns, m)?)?;
    m.add_function(wrap_pyfunction!(bam_stats, m)?)?;
    m.add_function(wrap_pyfunction!(bam_stats_bed, m)?)?;
    m.add_function(wrap_pyfunction!(parse_caller_vcf, m)?)?;
    m.add_function(wrap_pyfunction!(compute_tiers, m)?)?;
    m.add_function(wrap_pyfunction!(pileup_variants, m)?)?;
    m.add_function(wrap_pyfunction!(pileup_variants_multi, m)?)?;
    m.add_function(wrap_pyfunction!(coverage_bins, m)?)?;
    m.add_function(wrap_pyfunction!(wg_coverage_bins, m)?)?;
    Ok(())
}
