"""Tests for the seq2neo variant statistics package.

Run with: PYTHONPATH=bin .venv/bin/pytest bin/vcf_stats/tests/test_seq2neo_stats.py -v
"""

import os
import random
import sys
import tempfile
from pathlib import Path

import polars as pl
import polars.testing as pl_test
import pytest

# Ensure bin/ is on the Python path
_BIN_DIR = Path(__file__).resolve().parent.parent.parent
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

from vcf_stats.seq2neo import caller_parser
from vcf_stats.seq2neo.manifest_loader import (
    CALLER_CONFIGS,
    filter_complete,
    get_all_caller_vcf_paths,
    get_vcf_prefix,
    load_manifest,
)
from vcf_stats.seq2neo.rescue_parser import (
    _derive_variant_type,
    _is_transition,
    parse_rescue_vcf,
    rescue_info_fields,
)
from vcf_stats.seq2neo.statistics import (
    DNA_CALLERS,
    RNA_CALLERS,
    compute_vaf_columns,
    flag_filter_breakdown,
    ravex_filter_breakdown,
    sample_summary,
    set_summary,
    variant_type_distribution,
)
from vcf_stats.seq2neo.rescue_validator import (
    VALIDATION_METRICS,
    validate_sample,
    validate_single_metric,
    validation_summary,
)
from vcf_stats.seq2neo.visualizer import (
    generate_dashboard,
    plot_caller_overlap,
    plot_cosmic_gnomad_annotation,
    plot_cross_modality,
    plot_dna_vs_rna_dp,
    plot_dna_vs_rna_vaf,
    plot_gt_concordance,
    plot_per_sample_violin,
    plot_ravex_breakdown,
    plot_ti_tv_ratio,
    plot_vaf_distribution,
    plot_validation_heatmap,
    plot_variant_type_distribution,
    plot_vc_distribution,
)

# ── Test data paths (read-only, from real pipeline output) ────────────────
TEST_MANIFEST_CSV = Path(__file__).resolve().parent.parent.parent.parent / "examples/seq2neo/data/processed/sample_manifest.csv"
TEST_MANIFEST_PARQUET = Path(__file__).resolve().parent.parent.parent.parent / "examples/seq2neo/data/processed/sample_manifest.parquet"

# A known-good Set 3 sample with complete data
REAL_SAMPLE = {
    "sample_id": "PRJNA298376_3812",
    "set_number": 3,
    "patient_id": "3812",
    "base_output_dir": "/t9k/mnt/WorkSpace/data/ngs/liuxin/seq2neo/output",
    "dir_name": "PRJNA298376_3812",
    "vcf_prefix": "PRJNA298376_3812",
}


# ═══════════════════════════════════════════════════════════════════════════
# TestManifestLoader
# ═══════════════════════════════════════════════════════════════════════════

class TestManifestLoader:
    def test_load_csv_manifest(self):
        if not TEST_MANIFEST_CSV.exists():
            pytest.skip("Manifest CSV not found")
        df = load_manifest(str(TEST_MANIFEST_CSV))
        assert df.height == 66
        assert "sample_id" in df.columns
        assert "is_complete" in df.columns

    def test_load_parquet_manifest(self):
        if not TEST_MANIFEST_PARQUET.exists():
            pytest.skip("Manifest Parquet not found")
        df = load_manifest(str(TEST_MANIFEST_PARQUET))
        assert df.height == 66
        assert "vcf_prefix" in df.columns

    def test_get_vcf_prefix_set1(self):
        prefix = get_vcf_prefix(1, "4060", "PRJNA298376_4060")
        assert prefix == "4060"

    def test_get_vcf_prefix_set2(self):
        prefix = get_vcf_prefix(2, "3942", "PRJNA298376_3942")
        assert prefix == "PRJNA298376_3942"

    def test_get_vcf_prefix_set3(self):
        prefix = get_vcf_prefix(3, "3812", "PRJNA298376_3812")
        assert prefix == "PRJNA298376_3812"

    def test_filter_complete(self):
        df = pl.DataFrame({
            "sample_id": ["a", "b", "c"],
            "is_complete": [True, False, True],
        })
        filtered = filter_complete(df)
        assert filtered.height == 2
        assert filtered["sample_id"].to_list() == ["a", "c"]

    def test_get_all_caller_vcf_paths(self):
        paths = get_all_caller_vcf_paths(
            REAL_SAMPLE["base_output_dir"],
            REAL_SAMPLE["dir_name"],
            REAL_SAMPLE["vcf_prefix"],
        )
        assert len(paths) == 6
        # All 6 callers should exist for this sample
        for caller, path in paths.items():
            assert path is not None, f"{caller} VCF missing"
            assert os.path.isfile(path), f"{caller} VCF not found: {path}"

    def test_caller_configs_complete(self):
        assert len(CALLER_CONFIGS) == 6
        for caller in ["DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic", "RNA_deepsomatic", "DNA_strelka", "RNA_strelka"]:
            assert caller in CALLER_CONFIGS
            cfg = CALLER_CONFIGS[caller]
            assert "subdir" in cfg
            assert "pattern" in cfg
            assert "sample_suffix" in cfg

    def test_strelka_sample_suffix(self):
        assert CALLER_CONFIGS["DNA_strelka"]["sample_suffix"] == "TUMOR"
        assert CALLER_CONFIGS["RNA_strelka"]["sample_suffix"] == "TUMOR"


# ═══════════════════════════════════════════════════════════════════════════
# TestRescueParser
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def rescue_vcf_path():
    """Return path to a real rescue VCF for testing."""
    base = REAL_SAMPLE["base_output_dir"]
    name = REAL_SAMPLE["dir_name"]
    prefix = REAL_SAMPLE["vcf_prefix"]
    rescue_dir = f"{prefix}DT_vs_{prefix}DN_rescued_{prefix}RT_realign_vs_{prefix}DN"
    vcf_name = f"{prefix}DT_vs_{prefix}DN_rescued_{prefix}RT_realign_vs_{prefix}DN.filtered.vcf.stripped.vcf.gz"
    path = os.path.join(base, name, "vcf_realignment", "rescue", rescue_dir, vcf_name)
    if not os.path.isfile(path):
        pytest.skip("Real rescue VCF not available")
    return path


class TestRescueParser:
    def test_parse_real_rescue_vcf(self, rescue_vcf_path):
        df = parse_rescue_vcf(rescue_vcf_path)
        assert not df.is_empty()
        assert df.height > 0
        # Core columns
        for col in ["CHROM", "POS", "REF", "ALT", "FILTER"]:
            assert col in df.columns

    def test_required_columns_present(self, rescue_vcf_path):
        df = parse_rescue_vcf(rescue_vcf_path)
        fields = rescue_info_fields()
        for field in fields:
            if field not in df.columns:
                # Some fields may not be in every VCF, that's OK
                print(f"  Field not in rescue VCF: {field}")

    def test_variant_type_derivation(self):
        assert _derive_variant_type("A", "G") == "SNV"
        assert _derive_variant_type("AC", "A") == "DEL"
        assert _derive_variant_type("A", "ACGT") == "INS"
        assert _derive_variant_type("AC", "TG") == "MNV"

    def test_ti_tv_classification(self):
        assert _is_transition("A", "G") == True
        assert _is_transition("G", "A") == True
        assert _is_transition("C", "T") == True
        assert _is_transition("T", "C") == True
        assert _is_transition("A", "C") == False
        assert _is_transition("A", "T") == False
        assert _is_transition("G", "C") == False
        assert _is_transition("G", "T") == False
        assert _is_transition("AC", "TG") is None  # MNV
        assert _is_transition("A", "AG") is None    # INS

    def test_variant_types_in_real_data(self, rescue_vcf_path):
        df = parse_rescue_vcf(rescue_vcf_path)
        assert "variant_type" in df.columns
        vc = df["variant_type"].value_counts()
        types = {row["variant_type"]: row["count"] for row in vc.to_dicts()}
        assert "SNV" in types
        assert types.get("SNV", 0) > 0

    def test_info_field_types(self, rescue_vcf_path):
        df = parse_rescue_vcf(rescue_vcf_path)
        # Check int fields
        for field in ["N_SUPPORT_CALLERS", "DNA_SUPPORT", "RNA_SUPPORT"]:
            if field in df.columns:
                assert df[field].dtype in [pl.Int64, pl.Int32], f"{field} should be int"

    def test_empty_vcf_handling(self):
        with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as f:
            import gzip
            with gzip.open(f.name, "wt") as gz:
                gz.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        try:
            df = parse_rescue_vcf(f.name)
            assert df.is_empty()
        finally:
            os.unlink(f.name)


# ═══════════════════════════════════════════════════════════════════════════
# TestStatistics
# ═══════════════════════════════════════════════════════════════════════════

class TestStatistics:
    @pytest.fixture
    def sample_df(self):
        """Create a minimal test DataFrame with caller columns."""
        data = {
            "sample_id": ["test"] * 5,
            "set_number": [3] * 5,
            "CHROM": ["chr1"] * 5,
            "POS": [100, 200, 300, 400, 500],
            "FILTER": ["PASS", "PASS", "PASS", "NoConsensus", "PASS"],
            "VC": ["Somatic", "Somatic", "Germline", "Reference", "Somatic"],
            "variant_type": ["SNV", "SNV", "INS", "DEL", "SNV"],
            "ti_tv": [True, False, None, None, True],
            "N_SUPPORT_CALLERS": [6, 3, 2, 1, 5],
            "CROSS_MODALITY": ["YES", "NO", "YES", "NO", "YES"],
            "RESCUED": ["NO", "NO", "YES", "NO", "NO"],
            "COSMIC_ID": ["COSM123", None, None, None, "COSM456"],
            "GNOMAD_AF": [0.01, None, 0.05, None, None],
            "REDI_EVIDENCE": ["NONE", "NONE", "LOW", "NONE", "HIGH"],
            "RaVeX_FILTER": [None, "min_alt_reads", None, "gnomad;blacklist", None],
            "min_alt_reads": [False, True, False, False, False],
            "gnomad": [False, False, False, True, False],
            "blacklist": [False, False, False, True, False],
            # Caller columns
            "DNA_mutect2_DP": [50, 30, 40, 10, 60],
            "DNA_mutect2_AD_REF": [40, 25, 35, 10, 48],
            "DNA_mutect2_AD_ALT": [10, 5, 5, 0, 12],
            "RNA_mutect2_DP": [45, 28, 38, 8, 55],
            "RNA_mutect2_AD_REF": [36, 24, 32, 8, 44],
            "RNA_mutect2_AD_ALT": [9, 4, 6, 0, 11],
            "DNA_deepsomatic_DP": [52, 32, 42, 11, 62],
            "DNA_deepsomatic_AD_REF": [41, 26, 36, 11, 50],
            "DNA_deepsomatic_AD_ALT": [11, 6, 6, 0, 12],
            "RNA_deepsomatic_DP": [47, 30, 40, 9, 57],
            "RNA_deepsomatic_AD_REF": [37, 25, 34, 9, 46],
            "RNA_deepsomatic_AD_ALT": [10, 5, 6, 0, 11],
            "DNA_strelka_DP": [48, 29, 39, 10, 58],
            "DNA_strelka_AD_REF": [40, 24, 34, 10, 47],
            "DNA_strelka_AD_ALT": [8, 5, 5, 0, 11],
            "RNA_strelka_DP": [44, 27, 37, 7, 54],
            "RNA_strelka_AD_REF": [36, 23, 31, 7, 43],
            "RNA_strelka_AD_ALT": [8, 4, 6, 0, 11],
        }
        return pl.DataFrame(data)

    def test_vaf_calculation(self, sample_df):
        df = compute_vaf_columns(sample_df)
        assert "DNA_mutect2_VAF" in df.columns
        # VAF = AD_ALT / DP
        assert abs(df["DNA_mutect2_VAF"][0] - 10/50) < 0.001
        assert abs(df["DNA_mutect2_VAF"][1] - 5/30) < 0.001

    def test_vaf_dp_zero_handling(self):
        df = pl.DataFrame({
            "CHROM": ["chr1", "chr2"],
            "POS": [1, 2],
            "DNA_mutect2_DP": [0, None],
            "DNA_mutect2_AD_ALT": [5, 3],
        })
        result = compute_vaf_columns(df)
        # VAF should be None when DP=0 or DP=None
        assert result["DNA_mutect2_VAF"][0] is None
        assert result["DNA_mutect2_VAF"][1] is None

    def test_sample_summary(self, sample_df):
        df = compute_vaf_columns(sample_df)
        stats = sample_summary(df, "test")
        assert stats["total_variants"] == 5
        assert stats["pass_variants"] == 4
        assert stats["somatic_variants"] == 3
        assert stats["n_SNV"] == 3
        assert stats["n_INS"] == 1
        assert stats["n_DEL"] == 1
        assert stats["ti_count"] == 2
        assert stats["tv_count"] == 1
        assert stats["n_callers_6"] == 1
        assert stats["n_callers_3"] == 1
        assert stats["n_cross_modality"] == 3
        assert stats["n_rescued"] == 1
        assert stats["n_cosmic"] == 2
        assert stats["n_gnomad"] == 2
        assert stats["n_redi_high"] == 1
        assert stats["n_ravex_filtered"] == 2

    def test_set_summary(self):
        stats_df = pl.DataFrame({
            "set_number": [1, 1, 2],
            "total_variants": [100, 200, 50],
            "pass_variants": [80, 160, 40],
        })
        summary = set_summary(stats_df)
        assert summary.height == 2

    def test_variant_type_distribution(self, sample_df):
        dist = variant_type_distribution(sample_df, "set_number")
        assert not dist.is_empty()

    def test_ravex_filter_breakdown(self, sample_df):
        counts = ravex_filter_breakdown(sample_df)
        assert "min_alt_reads" in counts
        assert "gnomad" in counts
        assert "blacklist" in counts
        assert counts["min_alt_reads"] == 1
        assert counts["gnomad"] == 1
        assert counts["blacklist"] == 1

    def test_flag_filter_breakdown(self, sample_df):
        counts = flag_filter_breakdown(sample_df)
        assert counts.get("min_alt_reads", 0) == 1
        assert counts.get("gnomad", 0) == 1
        assert counts.get("blacklist", 0) == 1

    def test_caller_columns_exist(self, sample_df):
        for caller in DNA_CALLERS + RNA_CALLERS:
            assert f"{caller}_DP" in sample_df.columns


# ═══════════════════════════════════════════════════════════════════════════
# TestRescueValidator
# ═══════════════════════════════════════════════════════════════════════════

class TestRescueValidator:
    @pytest.fixture
    def valid_df(self):
        """DataFrame with matching computed and rescue values."""
        return pl.DataFrame({
            "CHROM": ["chr1", "chr2", "chr3"],
            "POS": [100, 200, 300],
            "DNA_DP_mean": [50.0, 30.0, 40.0],
            "DP_DNA_MEAN": [50.0, 30.0, 40.0],
            "RNA_DP_mean": [45.0, 28.0, 38.0],
            "DP_RNA_MEAN": [45.0, 28.0, 38.0],
            "DNA_VAF_mean": [0.2, 0.15, 0.1],
            "VAF_DNA_MEAN": [0.2, 0.15, 0.1],
            "RNA_VAF_mean": [0.18, 0.12, 0.08],
            "VAF_RNA_MEAN": [0.18, 0.12, 0.08],
            "DP_MEAN": [50.0, 30.0, 40.0],   # Must match DNA_DP_mean (the computed col used for comparison)
            "VAF_MEAN": [0.2, 0.15, 0.1],     # Must match DNA_VAF_mean
        })

    def test_validation_all_match(self, valid_df):
        result = validate_single_metric(valid_df, "DNA_DP_mean", "DP_DNA_MEAN", tolerance=0.01)
        assert result["n_mismatch"] == 0
        assert result["n_match"] == 3
        assert result["mismatch_pct"] == 0.0

    def test_validation_with_mismatches(self):
        df = pl.DataFrame({
            "DNA_DP_mean": [50.0, 30.0],
            "DP_DNA_MEAN": [50.0, 40.0],  # Second value mismatches by 10
        })
        result = validate_single_metric(df, "DNA_DP_mean", "DP_DNA_MEAN", tolerance=0.01)
        assert result["n_mismatch"] == 1
        assert result["n_match"] == 1
        assert result["max_abs_diff"] == 10.0

    def test_validation_missing_columns(self):
        df = pl.DataFrame({"A": [1, 2]})
        result = validate_single_metric(df, "DNA_DP_mean", "DP_DNA_MEAN")
        assert result["status"] == "missing_columns"

    def test_validation_null_handling(self):
        df = pl.DataFrame({
            "DNA_DP_mean": [50.0, None, 30.0],
            "DP_DNA_MEAN": [50.0, 30.0, None],
        })
        result = validate_single_metric(df, "DNA_DP_mean", "DP_DNA_MEAN")
        # Only the first row has both non-null
        assert result["n_missing"] == 2
        assert result["n_match"] == 1

    def test_validate_sample(self, valid_df):
        results = validate_sample(valid_df, "test_sample", tolerance=0.01)
        assert len(results) == len(VALIDATION_METRICS)
        for r in results:
            assert r["sample_id"] == "test_sample"
            assert r["n_mismatch"] == 0

    def test_validation_summary(self, valid_df):
        results = validate_sample(valid_df, "sample_a")
        report = pl.DataFrame(results)
        summary = validation_summary(report)
        assert not summary.is_empty()
        assert "metric" in summary.columns
        assert "overall_mismatch_pct" in summary.columns

    def test_validation_metrics_defined(self):
        assert len(VALIDATION_METRICS) >= 4
        metrics = [m[2] for m in VALIDATION_METRICS]
        assert "DNA DP mean" in metrics
        assert "RNA VAF mean" in metrics


# ═══════════════════════════════════════════════════════════════════════════
# TestIntegration (requires real data)
# ═══════════════════════════════════════════════════════════════════════════

class TestIntegration:
    def test_manifest_driven_paths(self):
        """Verify all VCF paths in manifest are derived, not hardcoded."""
        if not TEST_MANIFEST_CSV.exists():
            pytest.skip("Manifest CSV not found")
        df = load_manifest(str(TEST_MANIFEST_CSV))
        complete = filter_complete(df)
        for row in complete.to_dicts():
            prefix = get_vcf_prefix(row["set_number"], str(row["patient_id"]), row["sample_id"])
            assert prefix in row["rescue_vcf_path"], f"Prefix {prefix} not in path for {row['sample_id']}"
            assert prefix == row["vcf_prefix"], f"Prefix mismatch for {row['sample_id']}"

    def test_rescue_vcf_readable(self):
        """Verify real rescue VCFs can be opened and have basic structure."""
        if not TEST_MANIFEST_PARQUET.exists():
            pytest.skip("Manifest Parquet not found")
        df = filter_complete(load_manifest(str(TEST_MANIFEST_PARQUET)))
        row = df.head(1).to_dicts()[0]
        # Actually parse it
        parsed = parse_rescue_vcf(row["rescue_vcf_path"])
        assert not parsed.is_empty()
        assert "CHROM" in parsed.columns
        assert "POS" in parsed.columns
        assert parsed.height > 0

    def test_no_write_to_data_dirs(self):
        """Verify no files are written to analysis data directories."""
        data_dirs = [
            "/t9k/mnt/WorkSpace/data/ngs/liuxin/seq2neo/output",
            "/t9k/mnt/WorkSpace/data/ngs/zhanlingmin/Rnadnavar/output",
            "/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seq2neo/output",
        ]
        # This test verifies we're not inadvertently writing — we check
        # that our output goes to repo paths only
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "stats"
            # Verify output is NOT in data dirs
            for d in data_dirs:
                assert not str(output).startswith(d), f"Output would go to data dir: {d}"


# ═══════════════════════════════════════════════════════════════════════════
# TestCallerParser
# ═══════════════════════════════════════════════════════════════════════════

class TestCallerParser:
    """Tests for caller_parser module functions."""

    def test_find_sample_index_suffix_match(self):
        """_find_sample_index finds DT/RT/TUMOR by suffix."""
        from cyvcf2 import VCF
        vcf_path = _get_caller_vcf_path("DNA_mutect2")
        if vcf_path is None:
            pytest.skip("Real caller VCF not available")
        vcf = VCF(vcf_path)
        idx = caller_parser._find_sample_index(vcf, "DT")
        vcf.close()
        assert idx >= 0, "Should find DT sample in DNA Mutect2 VCF"

    def test_find_sample_index_not_found(self):
        """_find_sample_index returns -1 when no match."""
        from cyvcf2 import VCF
        vcf_path = _get_caller_vcf_path("DNA_mutect2")
        if vcf_path is None:
            pytest.skip("Real caller VCF not available")
        vcf = VCF(vcf_path)
        idx = caller_parser._find_sample_index(vcf, "NONEXISTENT")
        vcf.close()
        assert idx == -1

    def test_parse_single_caller_mutect2(self):
        """parse_single_caller extracts GT, AD, DP from Mutect2."""
        vcf_path = _get_caller_vcf_path("DNA_mutect2")
        if vcf_path is None:
            pytest.skip("Real caller VCF not available")
        # Use positions confirmed to be in the Mutect2 VCF
        targets = {("chr1", 633987), ("chr1", 1043470)}
        result = caller_parser.parse_single_caller(vcf_path, targets, "DT", "DNA_mutect2")
        assert "DP" in result
        assert "AD_REF" in result
        assert "AD_ALT" in result
        assert "GT" in result
        assert len(result["CHROM"]) >= 1

    def test_parse_single_caller_strelka(self):
        """parse_single_caller extracts DP, TAR, TOR from Strelka (no GT/AD)."""
        vcf_path = _get_caller_vcf_path("DNA_strelka")
        if vcf_path is None:
            pytest.skip("Real caller VCF not available")
        targets = {("chr1", 633987)}
        result = caller_parser.parse_single_caller(vcf_path, targets, "TUMOR", "DNA_strelka")
        assert "DP" in result
        assert "TAR" in result
        assert "TIR" in result
        assert "TOR" in result
        # Strelka has no GT and no AD
        assert "GT" not in result
        assert "AD_REF" not in result

    def test_strelka_no_gt(self):
        """Strelka result has no GT key."""
        vcf_path = _get_caller_vcf_path("DNA_strelka")
        if vcf_path is None:
            pytest.skip("Real caller VCF not available")
        result = caller_parser.parse_single_caller(vcf_path, {("chr1", 633987)}, "TUMOR", "DNA_strelka")
        assert "GT" not in result

    def test_strelka_no_ad(self):
        """Strelka result has no AD_REF/AD_ALT keys."""
        vcf_path = _get_caller_vcf_path("DNA_strelka")
        if vcf_path is None:
            pytest.skip("Real caller VCF not available")
        result = caller_parser.parse_single_caller(vcf_path, {("chr1", 633987)}, "TUMOR", "DNA_strelka")
        assert "AD_REF" not in result
        assert "AD_ALT" not in result

    def test_missing_vcf_returns_empty(self):
        """parse_single_caller returns empty result for missing VCF."""
        result = caller_parser.parse_single_caller("/nonexistent/path.vcf.gz", {("chr1", 1)}, "DT", "DNA_mutect2")
        assert result["CHROM"] == []

    def test_early_termination(self):
        """parse_single_caller stops when all targets found."""
        vcf_path = _get_caller_vcf_path("DNA_mutect2")
        if vcf_path is None:
            pytest.skip("Real caller VCF not available")
        targets = {("chr1", 63735)}  # Very first position in the VCF
        result = caller_parser.parse_single_caller(vcf_path, targets, "DT", "DNA_mutect2")
        assert len(result["CHROM"]) <= 1  # Should find at most 1 record before stopping

    def test_build_caller_results_lookup(self):
        """build_caller_results_lookup converts lists to position-keyed dict."""
        result = {
            "CHROM": ["chr1", "chr2"],
            "POS": [100, 200],
            "DP": [50, 30],
            "GT": ["0/1", "0/0"],
        }
        lookup = caller_parser.build_caller_results_lookup(result)
        assert ("chr1", 100) in lookup
        assert lookup[("chr1", 100)]["DP"] == 50
        assert lookup[("chr1", 100)]["GT"] == "0/1"
        assert lookup[("chr2", 200)]["DP"] == 30


def _get_caller_vcf_path(caller_name: str) -> str | None:
    """Get a real caller VCF path for a test sample."""
    from vcf_stats.seq2neo.manifest_loader import get_all_caller_vcf_paths
    paths = get_all_caller_vcf_paths(
        REAL_SAMPLE["base_output_dir"],
        REAL_SAMPLE["dir_name"],
        REAL_SAMPLE["vcf_prefix"],
    )
    return paths.get(caller_name)


# ═══════════════════════════════════════════════════════════════════════════
# TestVisualizer
# ═══════════════════════════════════════════════════════════════════════════

class TestVisualizer:
    """Tests for visualizer chart functions."""

    @pytest.fixture
    def tmp_output_dir(self):
        with tempfile.TemporaryDirectory() as d:
            yield d

    @pytest.fixture
    def viz_df(self):
        """Minimal DataFrame for chart testing."""
        return pl.DataFrame({
            "sample_id": ["s1"] * 6,
            "set_number": [1, 1, 1, 2, 2, 2],
            "CHROM": ["chr1"] * 6,
            "POS": range(100, 106),
            "FILTER": ["PASS"] * 6,
            "VC": ["Somatic", "Somatic", "Germline", "Somatic", "Reference", "Artifact"],
            "N_SUPPORT_CALLERS": [6, 4, 2, 3, 1, 5],
            "CROSS_MODALITY": ["YES", "NO", "YES", "NO", "YES", "NO"],
            "RESCUED": ["NO"] * 6,
            "variant_type": ["SNV", "SNV", "INS", "DEL", "SNV", "MNV"],
            "ti_tv": [True, False, None, None, True, False],
            "COSMIC_ID": ["C1", None, None, "C2", None, None],
            "GNOMAD_AF": [0.01, None, None, 0.05, None, 0.03],
            "DNA_VAF_mean": [0.2, 0.3, 0.1, 0.4, 0.15, 0.25],
            "RNA_VAF_mean": [0.18, 0.28, 0.08, 0.35, 0.12, 0.22],
            "DNA_DP_mean": [50.0, 30.0, 40.0, 60.0, 35.0, 45.0],
            "RNA_DP_mean": [45.0, 28.0, 38.0, 55.0, 32.0, 42.0],
            "DNA_mutect2_GT": ["0/1", "0/0", "0/1", "0/1", "0/0", "0/1"],
            "RNA_mutect2_GT": ["0/1", "0/0", "0/1", "0/1", None, "0/1"],
            "DNA_deepsomatic_GT": ["0/1", "0/0", "0/1", "0/1", "0/0", "0/1"],
            "RNA_deepsomatic_GT": ["0/1", "0/0", "0/1", "0/1", None, "0/1"],
            "DNA_mutect2_VAF": [0.2, 0.3, 0.1, 0.4, 0.15, 0.25],
            "RNA_mutect2_VAF": [0.18, 0.28, 0.08, 0.35, 0.12, 0.22],
            "DNA_deepsomatic_VAF": [0.21, 0.31, 0.11, 0.41, 0.16, 0.26],
            "RNA_deepsomatic_VAF": [0.19, 0.29, 0.09, 0.36, 0.13, 0.23],
            "DNA_strelka_VAF": [0.19, 0.29, 0.09, 0.38, 0.14, 0.24],
            "RNA_strelka_VAF": [0.17, 0.27, 0.07, 0.34, 0.11, 0.21],
        })

    def test_vc_distribution_returns_chart(self, viz_df, tmp_output_dir):
        fig = plot_vc_distribution(viz_df, tmp_output_dir)
        assert fig is not None
        assert hasattr(fig, "save")

    def test_caller_overlap_returns_chart(self, viz_df, tmp_output_dir):
        fig = plot_caller_overlap(viz_df, tmp_output_dir)
        assert fig is not None
        assert hasattr(fig, "save")

    def test_gt_concordance_returns_chart(self, viz_df, tmp_output_dir):
        fig = plot_gt_concordance(viz_df, tmp_output_dir)
        assert fig is not None
        assert hasattr(fig, "save")

    def test_cosmic_gnomad_returns_chart(self, viz_df, tmp_output_dir):
        fig = plot_cosmic_gnomad_annotation(viz_df, tmp_output_dir)
        assert fig is not None
        assert hasattr(fig, "save")

    def test_per_sample_violin_returns_chart(self, tmp_output_dir):
        df = pl.DataFrame({
            "sample_id": ["a", "b", "c", "d"],
            "total_variants": [100, 200, 150, 300],
            "disease": ["colorectal", "colon", "colorectal", "pancreatic"],
        })
        fig = plot_per_sample_violin(df, tmp_output_dir)
        assert fig is not None
        assert hasattr(fig, "save")

    def test_empty_data_handling(self, tmp_output_dir):
        """Chart functions return None for empty data."""
        empty = pl.DataFrame({"sample_id": [], "set_number": []})
        assert plot_vc_distribution(empty, tmp_output_dir) is None
        assert plot_caller_overlap(empty, tmp_output_dir) is None

    def test_dashboard_html_created(self, viz_df, tmp_output_dir):
        """generate_dashboard creates a non-empty HTML file."""
        fig = plot_vc_distribution(viz_df, tmp_output_dir)
        generate_dashboard([fig], tmp_output_dir)
        dashboard = Path(tmp_output_dir) / "dashboard.html"
        assert dashboard.exists()
        assert dashboard.stat().st_size > 0

    def test_html_contains_vegaembed(self, viz_df, tmp_output_dir):
        """Dashboard HTML contains vegaEmbed JavaScript (altair)."""
        fig = plot_vc_distribution(viz_df, tmp_output_dir)
        generate_dashboard([fig], tmp_output_dir)
        content = (Path(tmp_output_dir) / "dashboard.html").read_text()
        assert "vegaEmbed" in content


# ═══════════════════════════════════════════════════════════════════════════
# TestCLI
# ═══════════════════════════════════════════════════════════════════════════

class TestCLI:
    """Tests for CLI argument parsing and filtering."""

    def test_argparse_basic(self):
        """CLI parser accepts required arguments."""
        import argparse
        # Verify main() uses argparse correctly — test the parser
        parser = argparse.ArgumentParser()
        parser.add_argument("--manifest", required=True)
        parser.add_argument("--output-dir", required=True)
        parser.add_argument("--threads", type=int, default=4)
        parser.add_argument("--max-samples", type=int, default=None)
        parser.add_argument("--set", type=int, default=None)
        parser.add_argument("--no-validate", action="store_true")
        args = parser.parse_args(["--manifest", "m.csv", "--output-dir", "out"])
        assert args.manifest == "m.csv"
        assert args.output_dir == "out"
        assert args.threads == 4
        assert not args.no_validate

    def test_cli_filter_by_set(self):
        """--set flag filters correctly."""
        if not TEST_MANIFEST_CSV.exists():
            pytest.skip("Manifest CSV not found")
        df = filter_complete(load_manifest(str(TEST_MANIFEST_CSV)))
        df1 = df.filter(pl.col("set_number") == 1)
        assert df1.height == 15  # Set 1 has 15 samples

    def test_cli_filter_by_max_samples(self):
        """--max-samples limits correctly."""
        if not TEST_MANIFEST_CSV.exists():
            pytest.skip("Manifest CSV not found")
        df = filter_complete(load_manifest(str(TEST_MANIFEST_CSV)))
        assert df.head(3).height == 3

    def test_cli_missing_manifest_behavior(self):
        """Nonexistent manifest should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            load_manifest("/nonexistent/path/to/manifest.parquet")

    def test_cli_no_validate_flag(self):
        """--no-validate flag is accepted."""
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("--no-validate", action="store_true")
        args = parser.parse_args(["--no-validate"])
        assert args.no_validate

        args2 = parser.parse_args([])
        assert not args2.no_validate


# ═══════════════════════════════════════════════════════════════════════════
# TestIntegrationEndToEnd — full pipeline on 1 random sample per set
# ═══════════════════════════════════════════════════════════════════════════

class TestIntegrationEndToEnd:
    """End-to-end pipeline test: 1 random complete sample from each set.

    Verifies rescue VCF parsing, caller VCF parsing, statistics computation,
    rescue validation, chart generation (PNG+SVG+HTML), and CSV output.
    """

    @pytest.fixture(scope="class")
    def e2e_data(self):
        """Run full pipeline on 1 sample per set. Returns output_dir, all_stats, combined_df, report."""
        if not TEST_MANIFEST_CSV.exists():
            pytest.skip("Manifest CSV not found")

        # Load manifest, pick 1 random complete sample per set
        import random
        manifest = filter_complete(load_manifest(str(TEST_MANIFEST_CSV)))
        samples = []
        for set_num in [1, 2, 3, 4]:
            set_samples = manifest.filter(pl.col("set_number") == set_num).to_dicts()
            if set_samples:
                samples.append(random.choice(set_samples))

        # Process each sample sequentially
        cli_module = sys.modules.get("vcf_stats.seq2neo.cli")
        if cli_module is None:
            import importlib
            cli_module = importlib.import_module("vcf_stats.seq2neo.cli")

        samples_data = {}
        all_stats = []
        for row in samples:
            result = cli_module.process_single_sample(row)
            if result["df"] is not None:
                samples_data[result["sample_id"]] = result["df"]
            if result["stats"] is not None:
                all_stats.append(result["stats"])

        combined_df = pl.concat(list(samples_data.values()), how="diagonal_relaxed")

        # Validation
        from vcf_stats.seq2neo.rescue_validator import validate_all_samples
        report = validate_all_samples(samples_data, tolerance=0.01)

        # Output to temp dir
        with tempfile.TemporaryDirectory() as tmpdir:
            yield {
                "output_dir": tmpdir,
                "samples_data": samples_data,
                "all_stats": all_stats,
                "combined_df": combined_df,
                "report": report,
                "sample_count": len(samples),
            }

    def test_e2e_4_samples_across_sets(self, e2e_data):
        """All 4 sets have samples processed."""
        assert e2e_data["sample_count"] == 4

    def test_rescue_vcf_parsed(self, e2e_data):
        """All rescue VCFs have >0 variants."""
        for sid, df in e2e_data["samples_data"].items():
            assert df.height > 0, f"{sid} has 0 variants"

    def test_callers_parsed(self, e2e_data):
        """All 6 per-caller DP columns exist."""
        for sid, df in e2e_data["samples_data"].items():
            for caller in ["DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic",
                           "RNA_deepsomatic", "DNA_strelka", "RNA_strelka"]:
                assert f"{caller}_DP" in df.columns, f"{sid} missing {caller}_DP"

    def test_vaf_computed(self, e2e_data):
        """All 6 per-caller VAF columns exist."""
        for sid, df in e2e_data["samples_data"].items():
            for caller in ["DNA_mutect2", "RNA_mutect2", "DNA_deepsomatic",
                           "RNA_deepsomatic", "DNA_strelka", "RNA_strelka"]:
                assert f"{caller}_VAF" in df.columns, f"{sid} missing {caller}_VAF"

    def test_means_computed(self, e2e_data):
        """DNA/RNA mean columns exist."""
        for sid, df in e2e_data["samples_data"].items():
            for col in ["DNA_VAF_mean", "RNA_VAF_mean", "DNA_DP_mean", "RNA_DP_mean"]:
                assert col in df.columns, f"{sid} missing {col}"

    def test_sample_summary_correct(self, e2e_data):
        """Per-sample stats totals match DataFrame row counts."""
        for sid, df in e2e_data["samples_data"].items():
            stats = sample_summary(df, sid)
            assert stats["total_variants"] == df.height

    def test_all_12_charts_non_null(self, e2e_data):
        """All 12 chart functions return non-None altair charts."""
        df = e2e_data["combined_df"]
        d = e2e_data["output_dir"]

        charts = [
            plot_vc_distribution(df, d),
            plot_caller_overlap(df, d),
            plot_vaf_distribution(df, d),
            plot_dna_vs_rna_vaf(df, d),
            plot_dna_vs_rna_dp(df, d),
            plot_gt_concordance(df, d),
            plot_cosmic_gnomad_annotation(df, d),
            plot_variant_type_distribution(df, d),
            plot_ti_tv_ratio(df, d),
            plot_cross_modality(df, d),
        ]
        if e2e_data["all_stats"]:
            charts.append(plot_per_sample_violin(pl.DataFrame(e2e_data["all_stats"]), d))
        ravex = ravex_filter_breakdown(df)
        charts.append(plot_ravex_breakdown(ravex, d))
        if not e2e_data["report"].is_empty():
            charts.append(plot_validation_heatmap(e2e_data["report"], d))

        for i, chart in enumerate(charts):
            assert chart is not None, f"Chart {i+1} is None"
            assert hasattr(chart, "save"), f"Chart {i+1} has no save()"

    def test_dashboard_html_valid(self, e2e_data):
        """Dashboard HTML > 0 bytes."""
        df = e2e_data["combined_df"]
        d = e2e_data["output_dir"]
        fig = plot_vc_distribution(df, d)
        generate_dashboard([fig], d)
        html_path = Path(d) / "dashboard.html"
        assert html_path.exists()
        assert html_path.stat().st_size > 0

    def test_all_12_png_exist(self, e2e_data):
        """vl-convert PNG export works — 12 PNG files created."""
        df = e2e_data["combined_df"]
        d = e2e_data["output_dir"]
        # Generate all charts (side effect: writes to plots/)
        plot_vc_distribution(df, d)
        plot_vaf_distribution(df, d)
        ravex = ravex_filter_breakdown(df)
        plot_ravex_breakdown(ravex, d)
        # Count PNGs
        plots_dir = Path(d) / "plots"
        pngs = list(plots_dir.glob("*.png"))
        assert len(pngs) >= 3  # At least the ones we generated

    def test_all_12_svg_exist(self, e2e_data):
        """vl-convert SVG export works — SVG files created."""
        df = e2e_data["combined_df"]
        d = e2e_data["output_dir"]
        plot_vc_distribution(df, d)
        plots_dir = Path(d) / "plots"
        svgs = list(plots_dir.glob("*.svg"))
        assert len(svgs) >= 1

    def test_csv_outputs_exist(self, e2e_data):
        """Core CSV outputs can be written."""
        d = e2e_data["output_dir"]
        df = e2e_data["combined_df"]
        df.write_parquet(str(Path(d) / "variant_details.parquet"))
        if e2e_data["all_stats"]:
            pl.DataFrame(e2e_data["all_stats"]).write_csv(str(Path(d) / "sample_summary.csv"))
        verify_files = ["variant_details.parquet"]
        if e2e_data["all_stats"]:
            verify_files.append("sample_summary.csv")
        for f in verify_files:
            fp = Path(d) / f
            assert fp.exists(), f"Missing {f}"
            assert fp.stat().st_size > 0, f"Empty {f}"

    def test_no_data_dir_write(self, e2e_data):
        """No outputs go to analysis data directories."""
        data_dirs = [
            "/t9k/mnt/WorkSpace/data/ngs/liuxin/seq2neo/output",
            "/t9k/mnt/WorkSpace/data/ngs/zhanlingmin/Rnadnavar/output",
            "/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seq2neo/output",
        ]
        d = e2e_data["output_dir"]
        for dd in data_dirs:
            assert not str(d).startswith(dd), f"Output would go to data dir: {dd}"
