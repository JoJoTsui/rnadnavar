## Context

The current Python implementation (`rescue_parser.py`, `caller_parser.py`, `bam_stats.py`) uses cyvcf2/pysam for VCF/BAM parsing. Known issues: cyvcf2 + htslib fork-safety problems (already worked around with ThreadPoolExecutor), pysam pileup is too slow for variant-wise analysis (~9 hr for 65 samples), and the statistics are computed ad-hoc rather than systematically across 4 levels.

The pipeline needs to process 65 samples × (6 caller VCFs + 1 rescue VCF + 3 BAMs) = 650 files, with per-variant detail for 3.3M+ variants. The Rust+noodles approach gives 10-100× speedup for VCF/BAM I/O and enables single-pass statistics aggregation.

## Goals / Non-Goals

**Goals:**
- Replace cyvcf2 VCF parsing with noodles-vcf via maturin Python bindings
- Replace pysam BAM pileup with noodles-bam for variant-wise strand/quality analysis
- Implement 4-level statistics: dataset, per-tier, per-sample, per-sample×per-tier
- 3 BAM types per sample: DNA normal (DN), DNA tumor (DT), RNA tumor (RT)
- Two pileup modes: all variants (default), exclude NoConsensus FILTER
- All existing tests pass with the new implementation

**Non-Goals:**
- Replacing polars DataFrame operations (already fast)
- Replacing altair visualization (render-only, not a bottleneck)
- Replacing the tiering engine (CPU-light, already correct)
- Streaming/real-time processing (batch is sufficient)
- Cross-sample normalization or batch effects

## Decisions

### D0: Python/Rust boundary
**Choice**: Python handles orchestration (manifest loading, argument parsing, tiering integration), output generation (altair visualization, CSV/Parquet writing), and polars aggregation. Rust handles only compute-intensive VCF/BAM I/O — reading VCF INFO/FORMAT fields and BAM pileup at variant positions. The Rust crate returns Arrow record batches to Python; polars ingests them for the 4-level aggregation.
**Rationale**: Keeps visualization and tiering in Python where they're easy to maintain. Rust accelerates only the bottleneck (file I/O). Clear boundary: Rust reads files, Python does everything else.

### D1: Rust crate structure — single `stats_core` module
**Choice**: Single Rust crate `stats_core` (not `_core` or `seq2neo_core`) with submodules for VCF, BAM, and stats. All exposed via pyo3 as a single Python extension module. Import in Python as `from vcf_stats.seq2neo.stats_core import ...`.
**Dependency versions**: pyo3 0.28.3, pyo3-arrow 0.17.0, noodles 0.111, noodles-vcf 0.88, noodles-bam 0.90 — all latest stable releases, mutually compatible.
**Alternative**: Separate crates for VCF and BAM. Rejected — shared types (Position, VariantRecord) and BAM pileup needs VCF position lists.
**Rationale**: Keeps build simple (one `maturin build`), shared data structures, single import. Latest versions avoid dep conflicts that plagued earlier attempts.

### D2: VCF parsing — noodles-vcf single-pass
**Choice**: Rescue VCF INFO fields + caller FORMAT fields extracted in a single coordinated pass. Rescue VCF is scanned first for positions + INFO, then each caller VCF is scanned for FORMAT at those positions.
**Alternative**: Parse all VCFs independently then join. Rejected — caller VCFs have different variant sets; position-targeted extraction is more efficient.
**Rationale**: noodles-vcf provides streaming record iteration. We build a `HashSet<(String, u32)>` of target positions from the rescue VCF, then scan each caller VCF filtering to those positions.

### D3: BAM pileup — noodles-bam per-position fetch
**Choice**: For each (chrom, pos) from the variant list, fetch reads overlapping that position, then count allele support, strand orientation (F1R2/F2R1), base qualities, and mapping qualities.
**Alternative**: samtools mpileup subprocess. Rejected — process overhead, parsing text output, less control.
**Rationale**: noodles-bam's `QueryInterval` and read iterator give direct access to alignment data. The Rust implementation iterates reads at each position, computes per-allele strand counts and quality distributions, returns structured arrays.

### D4: Statistics engine — polars integration
**Choice**: Rust produces per-variant columns as Arrow arrays (via pyarrow interop), which polars ingests directly. Aggregation (4 levels) uses polars `group_by` + agg expressions in Python.
**Alternative**: Rust does all aggregation. Rejected — polars is already fast for group-by operations; rewriting in Rust adds complexity with little gain.
**Rationale**: The bottleneck is VCF/BAM I/O, not aggregation. polars can handle 3.3M rows × 200 columns efficiently. Rust produces clean columnar data; polars groups and aggregates.

### D5: 4-level statistics via polars (Python)
**Choice**: The Rust core produces per-variant DataFrames (via Arrow). The 4-level aggregation is done in Python using polars `group_by` + agg expressions. Level 1 (dataset) uses no grouping, Level 2 (tier) groups by final_tier, Level 3 (sample) groups by sample_id, Level 4 (sample-tier) groups by (sample_id, final_tier).
**Rationale**: polars is already fast for aggregation (C++/Rust backed). The bottleneck is VCF/BAM I/O, not group-by. Keeping aggregation in Python avoids duplicating metric definitions in Rust.

### D7: uv-managed environment
**Choice**: `uv` manages the entire Python environment — venv creation, package installation, and maturin integration. `maturin develop` runs within the uv-managed venv. The `pyproject.toml` lists maturin as a build dependency.
**Rationale**: Single tool for Python dependency management. uv is already the project standard. Maturin integrates with uv's venv for `maturin develop` builds.

### D6: Python fallback for VCF/BAM parsing
**Choice**: The Rust _core module is optional at import time. If `import _core` fails, fall back to existing cyvcf2/pysam-based parsers with a warning.
**Alternative**: Make Rust mandatory. Rejected — development convenience and CI environments may not have Rust.
**Rationale**: `try: from ._core import ... except ImportError: from ._py_fallback import ...` pattern. Production runs use Rust; development can use Python fallback.

## Risks / Trade-offs

- **maturin build complexity**: Devs need Rust+cargo installed. Mitigation: Python fallback, CI installs Rust via rustup.
- **noodles API stability**: noodles is pre-1.0, API may change. Mitigation: pin version in Cargo.toml, wrapper layer isolates noodles types.
- **Memory: 3.3M variants × ~100 columns**: ~3 GB in polars. Mitigation: process samples in batches, write parquet incrementally.
- **BAM index requirement**: noodles-bam needs .bai indices. Some BAM files may have stale indices (observed: "index older than data"). Mitigation: rebuild index if missing/stale (samtools index).

## Open Questions

- noodles-vcf 0.88 API for reading INFO and FORMAT fields — need to verify field access patterns with latest API
- noodles-bam 0.90 API for per-position read fetching — validate against the 0.78-era API used in examples
- BAM pileup strand counting algorithm — validate F1R2/F2R1 logic against samtools mpileup on test dataset
- PYO3_PYTHON env var must be set to the uv-managed venv Python path for cargo builds
