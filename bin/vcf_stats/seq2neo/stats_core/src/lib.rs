mod bam;
mod vcf;

use std::path::PathBuf;

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

/// Parse a rescue VCF and return a list of dicts (CHROM, POS, REF, ALT, FILTER, INFO...).
#[pyfunction]
fn parse_rescue(py: Python<'_>, path: String) -> PyResult<Bound<'_, PyList>> {
    let records = vcf::parse_rescue_vcf(&PathBuf::from(&path))
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
            format!("VCF parse error: {e}")
        ))?;

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

/// Compute whole-genome BAM statistics.
#[pyfunction]
fn bam_stats(py: Python<'_>, path: String, max_reads: u64) -> PyResult<Bound<'_, PyDict>> {
    let stats = bam::whole_genome_stats(&PathBuf::from(&path), max_reads)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
            format!("BAM stats error: {e}")
        ))?;

    let d = PyDict::new(py);
    d.set_item("total_reads", stats.total_reads)?;
    d.set_item("mapped_reads", stats.mapped_reads)?;
    d.set_item("mapping_rate", stats.mapping_rate)?;
    d.set_item("mean_insert_size", stats.mean_insert_size)?;
    d.set_item("mean_mapq", stats.mean_mapq)?;
    Ok(d)
}

/// stats_core — Rust-accelerated VCF/BAM parsing for seq2neo variant statistics.
#[pymodule]
fn stats_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_rescue, m)?)?;
    m.add_function(wrap_pyfunction!(bam_stats, m)?)?;
    Ok(())
}
