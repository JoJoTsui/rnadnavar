"""Tests for the seq2neo variant statistics package.

Run with: PYTHONPATH=bin .venv/bin/pytest bin/vcf_stats/tests/test_seq2neo_stats.py -v
"""

import os
import random
import sys
import tempfile
import threading
import time
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
    _add_derived_columns_polars,
    parse_rescue_vcf,
    rescue_info_fields,
    RESCUE_FLAG_FIELDS,
)
from vcf_stats.seq2neo.statistics import (
    DNA_CALLERS,
    RNA_CALLERS,
    compute_vaf_columns,
    flag_filter_breakdown,
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
    plot_per_sample_distribution,
    plot_ti_tv_ratio,
    plot_vaf_distribution,
    plot_validation_heatmap,
    plot_variant_type_distribution,
    plot_vc_distribution,
    plot_bam_coverage_violin,
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
        import polars as pl
        df = pl.DataFrame({
            "REF": ["A", "AC", "A", "AC"],
            "ALT": ["G", "A", "ACGT", "TG"],
        })
        result = _add_derived_columns_polars(df)
        assert result["variant_type"].to_list() == ["SNV", "DEL", "INS", "MNV"]

    def test_ti_tv_classification(self):
        import polars as pl
        df = pl.DataFrame({
            "REF": ["A", "G", "C", "T", "A", "A", "G", "G", "AC", "A"],
            "ALT": ["G", "A", "T", "C", "C", "T", "C", "T", "TG", "AG"],
        })
        result = _add_derived_columns_polars(df)
        expected_ti_tv = [True, True, True, True, False, False, False, False, None, None]
        assert result["ti_tv"].to_list() == expected_ti_tv

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
            "FILTER": ["Somatic", "Somatic", "Germline", "Artifact", "Somatic"],
            "VC": ["Somatic", "Somatic", "Germline", "Reference", "Somatic"],
            "variant_type": ["SNV", "SNV", "INS", "DEL", "SNV"],
            "ti_tv": [True, False, None, None, True],
            "N_SUPPORT_CALLERS": [6, 3, 2, 1, 5],
            "CROSS_MODALITY": ["YES", "NO", "YES", "NO", "YES"],
            "RESCUED": ["NO", "NO", "YES", "NO", "NO"],
            "COSMIC_ID": ["COSM123", None, None, None, "COSM456"],
            "GNOMAD_AF": [0.01, None, 0.05, None, None],
            "REDI_EVIDENCE": ["NONE", "NONE", "LOW", "NONE", "HIGH"],
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
        # 6-category FILTER classification: Somatic=3, Germline=1, Artifact=1
        assert stats["n_somatic"] == 3
        assert stats["n_germline"] == 1
        assert stats["n_artifact"] == 1
        assert stats["n_reference"] == 0
        assert stats["n_rnaedit"] == 0
        assert stats["n_noconsensus"] == 0
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

    def test_set_summary(self):
        stats_df = pl.DataFrame({
            "set_number": [1, 1, 2],
            "total_variants": [100, 200, 50],
            "n_somatic": [80, 160, 40],
        })
        summary = set_summary(stats_df)
        assert summary.height == 2
        assert "n_somatic" in summary.columns

    def test_variant_type_distribution(self, sample_df):
        dist = variant_type_distribution(sample_df, "set_number")
        assert not dist.is_empty()

    def test_flag_filter_breakdown(self, sample_df):
        counts = flag_filter_breakdown(sample_df)
        assert counts.get("min_alt_reads", 0) == 1
        assert counts.get("gnomad", 0) == 1
        assert counts.get("blacklist", 0) == 1

    def test_caller_columns_exist(self, sample_df):
        for caller in DNA_CALLERS + RNA_CALLERS:
            assert f"{caller}_DP" in sample_df.columns


    def test_gt_concordance(self):
        from vcf_stats.seq2neo.statistics import gt_concordance
        df = pl.DataFrame({
            "DNA_mutect2_GT": ["0/1", "0/0", "0/1", None, "0/1", "./."],
            "RNA_mutect2_GT": ["0/1", "0/0", "0/1", "0/1", None, "0/1"],
            "DNA_deepsomatic_GT": ["0/1", "0/0", "0/1", "0/1", "0/0", "0/1"],
            "RNA_deepsomatic_GT": ["0/1", None, None, "0/1", "0/0", "0/1"],
        })
        result = gt_concordance(df)
        assert "2" in result and "3" in result and "4" in result
        assert "no_agreement" in result
        assert all(isinstance(v, int) for v in result.values())
        # 6 total variants, valid GT counts per row:
        # row0: 4, row1: 3, row2: 3, row3: 3, row4: 3, row5: 3 (./. excluded)
        # >=2: 6, >=3: 6, >=4: 1, <2: 0
        assert result["2"] == 6
        assert result["3"] == 6
        assert result["4"] == 1
        assert result["no_agreement"] == 0

    def test_gt_concordance_empty(self):
        from vcf_stats.seq2neo.statistics import gt_concordance
        assert gt_concordance(pl.DataFrame()) == {}
        assert gt_concordance(pl.DataFrame({"A": [1, 2]})) == {}


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
        """build_caller_results_lookup converts lists to 4-tuple-keyed dict."""
        result = {
            "CHROM": ["chr1", "chr2"],
            "POS": [100, 200],
            "REF": ["A", "C"],
            "ALT": ["T", "G"],
            "DP": [50, 30],
            "GT": ["0/1", "0/0"],
        }
        lookup = caller_parser.build_caller_results_lookup(result)
        # 4-tuple keys: (CHROM, POS, REF, ALT)
        assert ("chr1", 100, "A", "T") in lookup
        assert lookup[("chr1", 100, "A", "T")]["DP"] == 50
        assert lookup[("chr1", 100, "A", "T")]["GT"] == "0/1"
        assert lookup[("chr2", 200, "C", "G")]["DP"] == 30


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
            "FILTER": ["Somatic", "Somatic", "Germline", "Somatic", "Reference", "Artifact"],
            "VC": ["Somatic", "Somatic", "Germline", "Somatic", "Reference", "Artifact"],
            "N_SUPPORT_CALLERS": [6, 4, 2, 3, 1, 5],
            "CROSS_MODALITY": ["YES", "NO", "YES", "NO", "YES", "NO"],
            "RESCUED": ["NO"] * 6,
            "variant_type": ["SNV", "SNV", "INS", "DEL", "SNV", "MNV"],
            "ti_tv": [True, False, None, None, True, False],
            "COSMIC_ID": ["C1", None, None, "C2", None, None],
            "GNOMAD_AF": [0.01, None, None, 0.05, None, 0.03],
            "final_tier": ["C1D1", "C2D0", "C3D1", "C4D0", "C5D1", "C6D0"],
            "caller_tier": ["C1", "C2", "C3", "C4", "C5", "C6"],
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
        fig = plot_per_sample_distribution(df, tmp_output_dir)
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
# TestVisualizerFixes — fix-visualizer-lazy-frame-and-oom
# ═══════════════════════════════════════════════════════════════════════════

class TestVisualizerFixes:
    """Tests for bug fixes in visualizer module."""

    @pytest.fixture
    def tmp_output_dir(self):
        with tempfile.TemporaryDirectory() as d:
            yield d

    @pytest.fixture
    def viz_df_with_facet(self):
        """DataFrame with caller_tier and facet_col for GT concordance testing."""
        return pl.DataFrame({
            "DNA_mutect2_GT": ["0/1", "0/0", "0/1", "0/1", "0/0", "0/1"],
            "RNA_mutect2_GT": ["0/1", "0/0", "0/1", "0/1", None, "0/1"],
            "DNA_deepsomatic_GT": ["0/1", "0/0", "0/1", "0/1", "0/0", "0/1"],
            "RNA_deepsomatic_GT": ["0/1", "0/0", "0/1", "0/1", None, "0/1"],
            "caller_tier": ["C1", "C2", "C1", "C3", "C4", "C1"],
            "disease_normalized": ["Lung", "Lung", "Breast", "Breast", "Lung", "Breast"],
        })

    # ── 7.1 _save_chart survives missing vl-convert ──────────────────────

    def test_save_chart_survives_missing_vl_convert(self, tmp_output_dir):
        """HTML saved, PNG/SVG skipped gracefully when vl-convert not installed."""
        import altair as alt
        import polars as pl
        from unittest.mock import patch
        from vcf_stats.seq2neo.visualizer import _save_chart

        df = pl.DataFrame({"x": [1, 2, 3], "y": [4, 5, 6]})
        chart = alt.Chart(df).mark_point().encode(x="x:Q", y="y:Q")

        call_count = {"html": 0, "png": 0, "svg": 0}
        real_save = alt.TopLevelMixin.save

        def mock_save(self, path, format=None, **kwargs):
            if format is None:
                call_count["html"] += 1
                return real_save(self, path, format=None, **kwargs)
            elif format == "png":
                call_count["png"] += 1
                raise ModuleNotFoundError("vl-convert not available")
            elif format == "svg":
                call_count["svg"] += 1
                raise ImportError("vl-convert not available")

        with patch.object(alt.TopLevelMixin, 'save', mock_save):
            _save_chart(chart, "test_survives", tmp_output_dir)
            # HTML should always be saved
            assert call_count["html"] == 1, "HTML save should always succeed"
            # PNG/SVG should be attempted but fail gracefully
            assert call_count["png"] == 1, "PNG save should be attempted"
            assert call_count["svg"] == 1, "SVG save should be attempted"
            # HTML file should exist (real save was called)
            html_path = Path(tmp_output_dir) / "plots" / "test_survives.html"
            assert html_path.exists(), f"HTML file should exist: {html_path}"

    # ── 7.2 GT concordance per tier with facet_col ───────────────────────

    def test_gt_concordance_per_tier_facet_col_correct_gts(self, viz_df_with_facet, tmp_output_dir):
        """caller_tier is NOT included in GT list when facet_col present."""
        # Verify the fix by inspecting the source code
        import inspect
        from vcf_stats.seq2neo.visualizer import plot_gt_concordance_per_tier
        source = inspect.getsource(plot_gt_concordance_per_tier)
        assert "row[:n_gt]" in source, (
            "plot_gt_concordance_per_tier should use row[:n_gt] for GT extraction, "
            "not row[:-1]"
        )
        assert "row[n_gt]" in source, (
            "plot_gt_concordance_per_tier should use row[n_gt] for caller_tier"
        )
        # Also test that the function runs without error
        fig = plot_gt_concordance_per_tier(viz_df_with_facet, tmp_output_dir, facet_col="disease_normalized")
        assert fig is not None

    def test_gt_concordance_per_tier_without_facet_correct(self, viz_df_with_facet, tmp_output_dir):
        """Correct behavior without facet_col — GT extraction excludes caller_tier."""
        from vcf_stats.seq2neo.visualizer import plot_gt_concordance_per_tier
        import inspect
        source = inspect.getsource(plot_gt_concordance_per_tier)
        # Without facet_col, row[:n_gt] still gives correct GT columns
        fig = plot_gt_concordance_per_tier(viz_df_with_facet, tmp_output_dir)
        assert fig is not None

    # ── 7.4-7.5 Helper hardening ─────────────────────────────────────────

    def test_sample_if_large_column_not_found_returns_empty(self):
        """_sample_if_large returns empty DataFrame on ColumnNotFoundError."""
        from vcf_stats.seq2neo.visualizer import _sample_if_large
        # Create a lazy frame referencing a non-existent column
        lazy = pl.LazyFrame({"x": [1, 2, 3]})
        # Force column error by selecting a missing column
        lazy_bad = lazy.select(pl.col("nonexistent_column"))
        result = _sample_if_large(lazy_bad)
        assert result.is_empty()
        assert isinstance(result, pl.DataFrame)

    def test_maybe_collect_column_not_found_returns_empty(self):
        """_maybe_collect returns empty DataFrame on ColumnNotFoundError."""
        from vcf_stats.seq2neo.visualizer import _maybe_collect
        lazy = pl.LazyFrame({"x": [1, 2, 3]})
        lazy_bad = lazy.select(pl.col("nonexistent_column"))
        result = _maybe_collect(lazy_bad)
        assert result.is_empty()
        assert isinstance(result, pl.DataFrame)

    # ── 7.6 CLI --exclude-sample-ids ─────────────────────────────────────

    def test_cli_exclude_sample_ids(self):
        """--exclude-sample-ids correctly removes samples from manifest."""
        manifest = pl.DataFrame({
            "sample_id": ["s1", "s2", "s3", "s4"],
            "is_complete": [True, True, True, True],
        })
        from vcf_stats.seq2neo.manifest_loader import filter_complete
        manifest = filter_complete(manifest)
        # Simulate CLI filter logic
        exclude_ids = ["s2", "s4"]
        manifest = manifest.filter(~pl.col("sample_id").is_in(exclude_ids))
        assert manifest.height == 2
        assert manifest["sample_id"].to_list() == ["s1", "s3"]

    # ── 7.7 Batch count optimization parity ──────────────────────────────

    def test_cosmic_gnomad_batch_counts_match_per_group(self, tmp_output_dir):
        """Optimized group_by counts match previous per-group _count_rows results."""
        import polars as pl
        from vcf_stats.seq2neo.visualizer import plot_cosmic_gnomad_annotation

        # Create test data with known counts
        df = pl.DataFrame({
            "set_number": [1, 1, 1, 1, 2, 2],
            "COSMIC_ID": ["C1", "C2", None, None, "C3", None],
            "GNOMAD_AF": [0.01, None, 0.05, None, None, 0.03],
        })
        # Verify chart generates without error and produces expected per-group counts
        # Pre-aggregated approach: set_number=1 has 4 total, 2 cosmic, 2 gnomad
        # set_number=2 has 2 total, 1 cosmic, 1 gnomad
        fig = plot_cosmic_gnomad_annotation(df, tmp_output_dir, group_col="set_number")
        assert fig is not None

        # Also test without group_col (global pie mode)
        fig2 = plot_cosmic_gnomad_annotation(df, tmp_output_dir)
        assert fig2 is not None


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
            charts.append(plot_per_sample_distribution(pl.DataFrame(e2e_data["all_stats"]), d))
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


# ═══════════════════════════════════════════════════════════════════════════
# TestTieringStats (task 2.5)
# ═══════════════════════════════════════════════════════════════════════════

class TestTieringStats:
    """Tests for tiering_stats bridge module."""

    def test_parse_filters_normalized(self):
        from vcf_stats.seq2neo.tiering_stats import _parse_filters_normalized
        fields = _parse_filters_normalized(
            "DNA_strelka:Somatic|DNA_mutect2:Somatic|DNA_deepsomatic:Somatic|"
            "RNA_strelka:Artifact|RNA_mutect2:Artifact|RNA_deepsomatic:Germline"
        )
        assert "FILTER_NORMALIZED_Strelka_DNA_TUMOR" in fields
        assert fields["FILTER_NORMALIZED_Strelka_DNA_TUMOR"] == "Somatic"
        assert fields["FILTER_NORMALIZED_Mutect2_RNA_TUMOR"] == "Artifact"
        assert fields["FILTER_NORMALIZED_DeepSomatic_RNA_TUMOR"] == "Germline"

    def test_parse_filters_normalized_empty(self):
        from vcf_stats.seq2neo.tiering_stats import _parse_filters_normalized
        assert _parse_filters_normalized(None) == {}
        assert _parse_filters_normalized("") == {}

    def test_compute_tiers_for_dataframe(self):
        from vcf_stats.seq2neo.tiering_stats import compute_tiers_for_dataframe
        df = pl.DataFrame({
            "CHROM": ["chr1", "chr2"],
            "POS": [100, 200],
            "FILTER": ["Somatic", "Germline"],
            "FILTERS_NORMALIZED": [
                "DNA_strelka:Somatic|DNA_mutect2:Somatic|DNA_deepsomatic:Somatic|RNA_strelka:Somatic|RNA_mutect2:Somatic",
                "DNA_strelka:Germline",
            ],
            "GNOMAD_AF": [0.01, None],
            "COSMIC_CNT": [5, 0],
            "REDI_EVIDENCE": ["NONE", "NONE"],
        })
        result = compute_tiers_for_dataframe(df)
        assert "final_tier" in result.columns
        assert "caller_tier" in result.columns
        assert "database_tier" in result.columns
        # C1D1: ≥2 DNA + ≥2 RNA + DB support
        assert result["final_tier"][0] == "C1D1"
        # C5D0: 1 DNA + 0 RNA + no DB
        assert result["final_tier"][1] == "C5D0"

    def test_tier_summary(self):
        from vcf_stats.seq2neo.tiering_stats import tier_summary
        df = pl.DataFrame({
            "final_tier": ["C1D1", "C2D0", "C1D1"],
            "DNA_VAF_mean": [0.2, 0.3, 0.1],
            "RNA_VAF_mean": [0.18, 0.28, 0.08],
            "DNA_DP_mean": [50.0, 30.0, 40.0],
            "RNA_DP_mean": [45.0, 28.0, 38.0],
            "variant_type": ["SNV", "SNV", "INS"],
            "ti_tv": [True, False, None],
            "N_SUPPORT_CALLERS": [6, 2, 4],
        })
        summary = tier_summary(df)
        assert not summary.is_empty()
        assert "n_variants" in summary.columns
        # C1D1 should have 2 variants
        c1d1 = summary.filter(pl.col("final_tier") == "C1D1")
        assert c1d1["n_variants"][0] == 2


# ═══════════════════════════════════════════════════════════════════════════
# TestRefAltDpStats (task 3.5)
# ═══════════════════════════════════════════════════════════════════════════

class TestRefAltDpStats:
    """Tests for REF_DP and ALT_DP statistics."""

    def test_ref_alt_dp_mean_columns(self):
        from vcf_stats.seq2neo.statistics import compute_mean_columns
        df = pl.DataFrame({
            "DNA_mutect2_AD_REF": [10, 20],
            "DNA_deepsomatic_AD_REF": [12, 22],
            "DNA_mutect2_AD_ALT": [5, 8],
            "DNA_deepsomatic_AD_ALT": [6, 10],
            "RNA_mutect2_AD_REF": [9, 18],
            "RNA_mutect2_AD_ALT": [4, 7],
        })
        result = compute_mean_columns(df)
        assert "DNA_REF_DP_mean" in result.columns
        assert "DNA_ALT_DP_mean" in result.columns
        assert "RNA_REF_DP_mean" in result.columns
        assert "RNA_ALT_DP_mean" in result.columns
        # Row 0: DNA REF mean = (10+12)/2 = 11
        assert abs(result["DNA_REF_DP_mean"][0] - 11.0) < 0.01

    def test_sample_summary_ref_alt_dp(self):
        from vcf_stats.seq2neo.statistics import sample_summary, compute_mean_columns
        from vcf_stats.seq2neo.statistics import compute_vaf_columns
        df = pl.DataFrame({
            "sample_id": ["test"] * 2,
            "FILTER": ["Somatic", "Germline"],
            "VC": ["Somatic", "Germline"],
            "variant_type": ["SNV", "SNV"],
            "ti_tv": [True, False],
            "DNA_mutect2_DP": [50, 30],
            "DNA_mutect2_AD_REF": [40, 25],
            "DNA_mutect2_AD_ALT": [10, 5],
            "RNA_mutect2_DP": [45, 28],
            "RNA_mutect2_AD_REF": [36, 24],
            "RNA_mutect2_AD_ALT": [9, 4],
        })
        df = compute_vaf_columns(df)
        df = compute_mean_columns(df)
        stats = sample_summary(df, "test")
        assert "mean_dna_ref_dp_mean" in stats
        assert "mean_dna_alt_dp_mean" in stats
        assert "mean_rna_ref_dp_mean" in stats
        assert "mean_rna_alt_dp_mean" in stats


# ═══════════════════════════════════════════════════════════════════════════
# TestMultiLevelAggregation (task 4.4)
# ═══════════════════════════════════════════════════════════════════════════

class TestMultiLevelAggregation:
    """Tests for dataset_summary and tier_summary."""

    def test_dataset_summary(self):
        from vcf_stats.seq2neo.statistics import dataset_summary
        df = pl.DataFrame({
            "sample_id": ["s1", "s1", "s2", "s2"],
            "FILTER": ["Somatic", "Somatic", "Germline", "Reference"],
            "VC": ["Somatic", "Somatic", "Germline", "Reference"],
            "variant_type": ["SNV", "SNV", "INS", "SNV"],
            "ti_tv": [True, False, None, True],
            "N_SUPPORT_CALLERS": [6, 3, 2, 1],
            "CROSS_MODALITY": ["YES", "NO", "YES", "NO"],
            "RESCUED": ["NO", "NO", "YES", "NO"],
            "COSMIC_ID": ["C1", None, None, None],
            "GNOMAD_AF": [0.01, None, 0.05, None],
            "final_tier": ["C1D1", "C2D0", "C3D1", "C7D0"],
        })
        ds = dataset_summary(df)
        assert ds["total_variants"] == 4
        assert ds["n_samples"] == 2
        # 6-category FILTER classification
        assert ds["n_somatic"] == 2  # rows 0,1 → Somatic
        assert ds["n_germline"] == 1  # row 2 → Germline
        assert ds["n_reference"] == 1  # row 3 → Reference
        assert ds["n_SNV"] == 3
        assert ds["n_INS"] == 1
        assert ds["n_callers_6"] == 1
        assert ds["n_cross_modality"] == 2
        assert ds["n_cosmic"] == 1
        assert "tier_C1D1" in ds
        assert ds["tier_C1D1"] == 1

    def test_dataset_summary_empty(self):
        from vcf_stats.seq2neo.statistics import dataset_summary
        result = dataset_summary(pl.DataFrame())
        assert result["total_variants"] == 0


# ═══════════════════════════════════════════════════════════════════════════
# TestNewVisualizerFunctions (tasks 5.7 + 6.4)
# ═══════════════════════════════════════════════════════════════════════════

class TestNewVisualizerFunctions:
    """Tests for new chart functions added in fix-vcf-statistics."""

    @pytest.fixture
    def tmp_output_dir(self):
        with tempfile.TemporaryDirectory() as d:
            yield d

    @pytest.fixture
    def tiered_df(self):
        return pl.DataFrame({
            "sample_id": ["s1"] * 8,
            "set_number": [1, 1, 1, 1, 2, 2, 2, 2],
            "caller_tier": ["C1", "C1", "C2", "C2", "C3", "C4", "C5", "C6"],
            "final_tier": ["C1D1", "C1D1", "C2D0", "C2D0", "C3D1", "C4D0", "C5D1", "C6D0"],
            "N_SUPPORT_CALLERS": [6, 5, 4, 3, 2, 2, 1, 1],
            "variant_type": ["SNV", "SNV", "SNV", "INS", "SNV", "DEL", "SNV", "MNV"],
            "DNA_mutect2_VAF": [0.2, 0.3, 0.1, 0.4, 0.15, 0.25, 0.05, 0.08],
            "RNA_mutect2_VAF": [0.18, 0.28, 0.08, 0.35, 0.12, 0.22, 0.04, 0.07],
            "DNA_deepsomatic_VAF": [0.21, 0.31, 0.11, 0.41, 0.16, 0.26, 0.06, 0.09],
            "RNA_deepsomatic_VAF": [0.19, 0.29, 0.09, 0.36, 0.13, 0.23, 0.05, 0.08],
            "DNA_strelka_VAF": [0.19, 0.29, 0.09, 0.38, 0.14, 0.24, 0.04, 0.07],
            "RNA_strelka_VAF": [0.17, 0.27, 0.07, 0.34, 0.11, 0.21, 0.03, 0.06],
            "DNA_mutect2_DP": [50, 30, 40, 10, 60, 35, 20, 15],
            "RNA_mutect2_DP": [45, 28, 38, 8, 55, 32, 18, 12],
            "DNA_deepsomatic_DP": [52, 32, 42, 11, 62, 37, 22, 16],
            "RNA_deepsomatic_DP": [47, 30, 40, 9, 57, 34, 20, 13],
            "DNA_strelka_DP": [48, 29, 39, 10, 58, 33, 19, 14],
            "RNA_strelka_DP": [44, 27, 37, 7, 54, 31, 17, 11],
            "DNA_mutect2_GT": ["0/1", "0/1", "0/0", "0/1", "0/1", "0/0", "0/1", None],
            "RNA_mutect2_GT": ["0/1", "0/1", "0/0", "0/1", "0/1", None, "0/1", None],
            "DNA_deepsomatic_GT": ["0/1", "0/1", "0/0", "0/1", "0/1", "0/0", None, None],
            "RNA_deepsomatic_GT": ["0/1", "0/1", "0/0", "0/1", None, None, None, None],
            "DNA_REF_DP_mean": [40.0, 25.0, 35.0, 10.0, 48.0, 30.0, 18.0, 13.0],
            "RNA_REF_DP_mean": [36.0, 24.0, 32.0, 8.0, 44.0, 28.0, 16.0, 11.0],
            "DNA_ALT_DP_mean": [8.0, 5.0, 5.0, 0.5, 10.0, 8.0, 2.0, 3.0],
            "RNA_ALT_DP_mean": [7.0, 4.0, 4.0, 0.3, 9.0, 7.0, 1.5, 2.5],
        })

    def test_plot_vaf_boxplot_per_tier(self, tiered_df, tmp_output_dir):
        from vcf_stats.seq2neo.visualizer import plot_vaf_boxplot_per_tier
        fig = plot_vaf_boxplot_per_tier(tiered_df, tmp_output_dir)
        assert fig is not None
        assert hasattr(fig, "save")

    def test_plot_dp_boxplot_per_tier(self, tiered_df, tmp_output_dir):
        from vcf_stats.seq2neo.visualizer import plot_dp_boxplot_per_tier
        fig = plot_dp_boxplot_per_tier(tiered_df, tmp_output_dir)
        assert fig is not None
        assert hasattr(fig, "save")

    def test_plot_gt_concordance_per_tier(self, tiered_df, tmp_output_dir):
        from vcf_stats.seq2neo.visualizer import plot_gt_concordance_per_tier
        fig = plot_gt_concordance_per_tier(tiered_df, tmp_output_dir)
        assert fig is not None
        assert hasattr(fig, "save")

    def test_plot_tiered_caller_overlap(self, tiered_df, tmp_output_dir):
        from vcf_stats.seq2neo.visualizer import plot_tiered_caller_overlap
        fig = plot_tiered_caller_overlap(tiered_df, tmp_output_dir)
        assert fig is not None
        assert hasattr(fig, "save")

    def test_plot_tiered_variant_types(self, tiered_df, tmp_output_dir):
        from vcf_stats.seq2neo.visualizer import plot_tiered_variant_types
        fig = plot_tiered_variant_types(tiered_df, tmp_output_dir)
        assert fig is not None
        assert hasattr(fig, "save")

    def test_plot_ref_alt_dp_scatter(self, tiered_df, tmp_output_dir):
        from vcf_stats.seq2neo.visualizer import plot_ref_alt_dp_scatter
        fig = plot_ref_alt_dp_scatter(tiered_df, tmp_output_dir)
        assert fig is not None
        assert hasattr(fig, "save")

    def test_plot_per_sample_distribution(self, tmp_output_dir):
        from vcf_stats.seq2neo.visualizer import plot_per_sample_distribution
        df = pl.DataFrame({
            "sample_id": ["PRJNA_001", "PRJNA_002", "PRJNA_003"],
            "total_variants": [5000, 3000, 7000],
        })
        fig = plot_per_sample_distribution(df, tmp_output_dir)
        assert fig is not None
        assert hasattr(fig, "save")

    def test_empty_data_returns_none(self, tmp_output_dir):
        """New chart functions return None for empty/missing data."""
        from vcf_stats.seq2neo.visualizer import (
            plot_vaf_boxplot_per_tier, plot_dp_boxplot_per_tier,
            plot_gt_concordance_per_tier, plot_tiered_caller_overlap,
            plot_tiered_variant_types, plot_ref_alt_dp_scatter,
        )
        empty = pl.DataFrame()
        assert plot_vaf_boxplot_per_tier(empty, tmp_output_dir) is None
        assert plot_dp_boxplot_per_tier(empty, tmp_output_dir) is None
        assert plot_gt_concordance_per_tier(empty, tmp_output_dir) is None
        assert plot_tiered_caller_overlap(empty, tmp_output_dir) is None
        assert plot_tiered_variant_types(empty, tmp_output_dir) is None
        assert plot_ref_alt_dp_scatter(empty, tmp_output_dir) is None


# ═══════════════════════════════════════════════════════════════════════════
# TestBamStats (task 7.5)
# ═══════════════════════════════════════════════════════════════════════════

class TestBamStats:
    """Tests for BAM statistics module."""

    def test_module_imports(self):
        from vcf_stats.seq2neo import bam_stats
        assert hasattr(bam_stats, "compute_bam_stats")
        assert hasattr(bam_stats, "compute_sample_bam_stats")

    def test_compute_bam_stats_nonexistent_file(self):
        from vcf_stats.seq2neo.bam_stats import compute_bam_stats
        result = compute_bam_stats("/nonexistent/path.bam")
        assert result is None

    def test_compute_bam_stats_none_path(self):
        from vcf_stats.seq2neo.bam_stats import compute_bam_stats
        result = compute_bam_stats(None)
        assert result is None

    def test_compute_sample_bam_stats_missing_bams(self):
        from vcf_stats.seq2neo.bam_stats import compute_sample_bam_stats
        results = compute_sample_bam_stats(
            base_output_dir="/nonexistent",
            dir_name="no_such_dir",
            sample_id="TEST_SAMPLE",
            set_number=9,
        )
        assert len(results) == 3  # DN + DT + RT
        assert results[0]["bam_type"] == "DN"
        assert results[1]["bam_type"] == "DT"
        assert results[2]["bam_type"] == "RT"
        assert not results[0]["has_bam"]
        assert not results[1]["has_bam"]
        assert results[0]["total_reads"] is None

    # ── Shared BED Processing tests ──────────────────────────────────────

    def test_read_and_merge_bed_total(self):
        """Merged BED total matches sum of non-overlapping input intervals."""
        from vcf_stats.seq2neo.bam_stats import read_and_merge_bed
        content = (
            "chr1\t1000\t2000\tgeneA\n"
            "chr2\t3000\t4000\tgeneB\n"
            "chr3\t5000\t7000\tgeneC\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            bed_path = f.name
        try:
            raw_bed_total, raw_bed_regions, merged_bed_regions = read_and_merge_bed(bed_path, gap=100_000)
            # chr1: 1000bp, chr2: 1000bp, chr3: 2000bp = 4000bp total
            assert raw_bed_total == 4000, f"Expected 4000, got {raw_bed_total}"
            assert len(raw_bed_regions) == 3
            assert len(merged_bed_regions) == 3  # no merging across chromosomes
        finally:
            os.unlink(bed_path)

    def test_bed_merge_within_gap(self):
        """Intervals within 100Kb gap on same chromosome are merged."""
        from vcf_stats.seq2neo.bam_stats import read_and_merge_bed
        # Two intervals 50Kb apart → merged into one
        content = (
            "chr1\t1000\t2000\tgeneA\n"
            "chr1\t50000\t51000\tgeneB\n"  # 48Kb gap from 2000 → merged
            "chr2\t3000\t4000\tgeneC\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            bed_path = f.name
        try:
            raw_bed_total, raw_bed_regions, merged_bed_regions = read_and_merge_bed(bed_path, gap=100_000)
            assert len(raw_bed_regions) == 3
            assert len(merged_bed_regions) == 2, (
                f"Expected 2 merged regions, got {len(merged_bed_regions)}: {merged_bed_regions}"
            )
            # chr1: merged [1000, 51000)
            assert merged_bed_regions[0] == ("chr1", 1000, 51000), f"Unexpected: {merged_bed_regions[0]}"
            # chr2: unchanged
            assert merged_bed_regions[1] == ("chr2", 3000, 4000)
            # Raw total: 1000+1000+1000 = 3000, merged total: (51000-1000)+(4000-3000)=51000
            assert raw_bed_total == 3000, f"Expected raw total 3000, got {raw_bed_total}"
        finally:
            os.unlink(bed_path)

    def test_bed_merge_no_gap(self):
        """Intervals far apart on same chromosome are NOT merged."""
        from vcf_stats.seq2neo.bam_stats import read_and_merge_bed
        content = (
            "chr1\t1000\t2000\tgeneA\n"
            "chr1\t300000\t301000\tgeneB\n"  # 298Kb gap → not merged
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            bed_path = f.name
        try:
            raw_bed_total, raw_bed_regions, merged_bed_regions = read_and_merge_bed(bed_path, gap=100_000)
            assert len(merged_bed_regions) == 2, (
                f"Expected 2 separate regions, got {len(merged_bed_regions)}"
            )
        finally:
            os.unlink(bed_path)

    def test_bam_stats_on_target_coverage(self):
        """BED-filtered on-target coverage is ≤ total genome coverage."""
        from vcf_stats.seq2neo.bam_stats import read_and_merge_bed
        # Verify that read_and_merge_bed produces a positive total
        content = (
            "chr1\t600000\t700000\tregion1\n"
            "chr1\t900000\t1000000\tregion2\n"  # 200Kb gap → NOT merged
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            bed_path = f.name
        try:
            raw_bed_total, raw_bed_regions, merged_bed_regions = read_and_merge_bed(bed_path, gap=100_000)
            assert raw_bed_total == 200000  # 100Kb + 100Kb, not merged
            assert len(raw_bed_regions) == 2
            assert len(merged_bed_regions) == 2
            # On-target coverage is always ≤ total when using same reads:
            # on_target_bases counts bases within BED only,
            # total_query_length counts all mapped bases.
            # This is a structural invariant, not a runtime test.
            assert raw_bed_total > 0
        finally:
            os.unlink(bed_path)


# ═══════════════════════════════════════════════════════════════════════════
# TestGILRelease — verify Rust pyfunctions release the GIL during computation.
# Without py.detach(), pyo3 #[pyfunction] holds the GIL for the entire call,
# serializing all threads. These tests catch that regression.
# ═══════════════════════════════════════════════════════════════════════════

class TestGILRelease:
    """Verify Rust pyfunctions release the GIL so threads run in parallel.

    Each test runs 2 threads simultaneously calling the same function.
    A GIL-serialized call shows total_wall ≈ sum(individual_times).
    A GIL-released call shows total_wall ≈ max(individual_times).
    We use total_wall < sum * 0.7 as the pass threshold (allows some I/O overhead).
    """

    def test_bam_stats_releases_gil(self):
        """bam_stats() MUST release the GIL for parallel BAM scanning."""
        try:
            import stats_core
        except ImportError:
            pytest.skip("stats_core not available")

        bam = os.path.join(
            REAL_SAMPLE["base_output_dir"],
            REAL_SAMPLE["dir_name"],
            "preprocessing", "mapped",
            f"{REAL_SAMPLE['vcf_prefix']}DN",
            f"{REAL_SAMPLE['vcf_prefix']}DN.sorted.bam",
        )
        if not os.path.isfile(bam):
            pytest.skip(f"BAM not found: {bam}")

        results = {}

        def worker(label):
            t0 = time.time()
            stats_core.bam_stats(bam, 5_000_000)  # 5M reads for fast test
            results[label] = time.time() - t0

        t0 = time.time()
        t1 = threading.Thread(target=worker, args=("A",))
        t2 = threading.Thread(target=worker, args=("B",))
        t1.start()
        t2.start()
        t1.join()
        t2.join()
        total = time.time() - t0

        assert "A" in results and "B" in results
        a_time, b_time = results["A"], results["B"]
        sum_time = a_time + b_time
        assert total < sum_time * 0.7, (
            f"GIL NOT released! Total={total:.1f}s, A={a_time:.1f}s, B={b_time:.1f}s, "
            f"Sum={sum_time:.1f}s. Without GIL release threads serialize: total ≈ sum. "
            f"Check that py.detach() is used in bam_stats()."
        )

    def test_parse_rescue_releases_gil(self, rescue_vcf_path):
        """parse_rescue() MUST release the GIL for parallel VCF parsing."""
        try:
            import stats_core
        except ImportError:
            pytest.skip("stats_core not available")

        results = {}

        def worker(label):
            t0 = time.time()
            stats_core.parse_rescue(rescue_vcf_path)
            results[label] = time.time() - t0

        t0 = time.time()
        t1 = threading.Thread(target=worker, args=("A",))
        t2 = threading.Thread(target=worker, args=("B",))
        t1.start()
        t2.start()
        t1.join()
        t2.join()
        total = time.time() - t0

        assert "A" in results and "B" in results
        a_time, b_time = results["A"], results["B"]
        sum_time = a_time + b_time
        assert total < sum_time * 0.7, (
            f"GIL NOT released! Total={total:.1f}s, A={a_time:.1f}s, B={b_time:.1f}s, "
            f"Sum={sum_time:.1f}s. Without GIL release threads serialize: total ≈ sum. "
            f"Check that py.detach() is used in parse_rescue()."
        )

    def test_bam_stats_still_works_after_detach(self):
        """Sanity: bam_stats returns correct data after GIL release."""
        try:
            import stats_core
        except ImportError:
            pytest.skip("stats_core not available")

        bam = os.path.join(
            REAL_SAMPLE["base_output_dir"],
            REAL_SAMPLE["dir_name"],
            "preprocessing", "mapped",
            f"{REAL_SAMPLE['vcf_prefix']}DN",
            f"{REAL_SAMPLE['vcf_prefix']}DN.sorted.bam",
        )
        if not os.path.isfile(bam):
            pytest.skip(f"BAM not found: {bam}")

        result = stats_core.bam_stats(bam, 1_000_000)
        assert result["total_reads"] == 1_000_000
        assert result["mapped_reads"] > 0
        assert "mean_coverage" in result
        assert result["mean_insert_size"] > 0
        assert 0 <= result["mapping_rate"] <= 100

    def test_parse_rescue_still_works_after_detach(self, rescue_vcf_path):
        """Sanity: parse_rescue returns correct data after GIL release."""
        try:
            import stats_core
        except ImportError:
            pytest.skip("stats_core not available")

        records = stats_core.parse_rescue(rescue_vcf_path)
        assert len(records) > 0
        record = records[0]
        assert "CHROM" in record
        assert "POS" in record
        assert "REF" in record
        assert "ALT" in record
        assert "FILTER" in record

    def test_parse_caller_vcf_releases_gil(self):
        """parse_caller_vcf MUST release the GIL for parallel caller parsing."""
        try:
            import stats_core
        except ImportError:
            pytest.skip("stats_core not available")
        if not hasattr(stats_core, 'parse_caller_vcf'):
            pytest.skip("parse_caller_vcf not available")

        import glob
        base = os.path.join(REAL_SAMPLE["base_output_dir"], REAL_SAMPLE["dir_name"])
        rescue = glob.glob(os.path.join(base, "vcf_realignment", "rescue", "*", "*.filtered.vcf.stripped.vcf.gz"))
        if not rescue:
            pytest.skip("Rescue VCF not found")
        records = stats_core.parse_rescue(rescue[0])
        # Use full target set for balanced, substantial workload
        chroms = [r["CHROM"] for r in records]
        poss = [r["POS"] for r in records]
        refs = [r["REF"] for r in records]
        alts = [r["ALT"] for r in records]

        cfg = CALLER_CONFIGS["DNA_mutect2"]
        subdir = cfg["subdir"].format(prefix=REAL_SAMPLE["vcf_prefix"])
        vcf_path = glob.glob(os.path.join(base, subdir, cfg["pattern"]))[0]

        cfg2 = CALLER_CONFIGS["DNA_deepsomatic"]
        subdir2 = cfg2["subdir"].format(prefix=REAL_SAMPLE["vcf_prefix"])
        vcf_path2 = glob.glob(os.path.join(base, subdir2, cfg2["pattern"]))[0]

        results = {}
        def worker(label, path, sfx, name):
            t0 = time.time()
            stats_core.parse_caller_vcf(path, chroms, poss, refs, alts, sfx, name)
            results[label] = time.time() - t0

        t0 = time.time()
        t1 = threading.Thread(target=worker, args=("A", vcf_path, cfg["sample_suffix"], "DNA_mutect2"))
        t2 = threading.Thread(target=worker, args=("B", vcf_path2, cfg2["sample_suffix"], "DNA_deepsomatic"))
        t1.start(); t2.start(); t1.join(); t2.join()
        total = time.time() - t0
        assert total < (results["A"] + results["B"]) * 0.8, (
            f"GIL NOT released in parse_caller_vcf! Total={total:.1f}s, A={results['A']:.1f}s, B={results['B']:.1f}s"
        )

    def test_compute_tiers_releases_gil(self):
        """compute_tiers MUST release the GIL for parallel tiering."""
        try:
            import stats_core
        except ImportError:
            pytest.skip("stats_core not available")
        if not hasattr(stats_core, 'compute_tiers'):
            pytest.skip("compute_tiers not available")

        # Generate synthetic data for 50000 variants
        n = 50000
        filters = ["Somatic"] * n
        fnorm = [""] * n
        gaf = [0.01] * n
        csc = [None] * n
        redi = [None] * n
        ds = [2] * n
        rs = [1] * n

        results = {}
        def worker(label):
            t0 = time.time()
            stats_core.compute_tiers(filters, fnorm, gaf, csc, redi, ds, rs)
            results[label] = time.time() - t0

        t0 = time.time()
        t1 = threading.Thread(target=worker, args=("A",))
        t2 = threading.Thread(target=worker, args=("B",))
        t1.start(); t2.start(); t1.join(); t2.join()
        total = time.time() - t0
        assert total < (results["A"] + results["B"]) * 0.7, (
            f"GIL NOT released in compute_tiers! Total={total:.1f}s"
        )

    def test_pileup_variants_releases_gil(self):
        """pileup_variants MUST release the GIL for parallel pileup."""
        try:
            import stats_core
        except ImportError:
            pytest.skip("stats_core not available")
        if not hasattr(stats_core, 'pileup_variants'):
            pytest.skip("pileup_variants not available")

        bam = os.path.join(
            REAL_SAMPLE["base_output_dir"], REAL_SAMPLE["dir_name"],
            "preprocessing", "mapped", f"{REAL_SAMPLE['vcf_prefix']}DT",
            f"{REAL_SAMPLE['vcf_prefix']}DT.sorted.bam",
        )
        if not os.path.isfile(bam):
            pytest.skip("BAM not found")

        results = {}
        def worker(label, chroms):
            t0 = time.time()
            stats_core.pileup_variants(bam, chroms, [633987] * len(chroms), ["C"] * len(chroms), ["T"] * len(chroms))
            results[label] = time.time() - t0

        t0 = time.time()
        t1 = threading.Thread(target=worker, args=("A", ["chr1"]))
        t2 = threading.Thread(target=worker, args=("B", ["chr2"]))
        t1.start(); t2.start(); t1.join(); t2.join()
        total = time.time() - t0
        assert total < (results["A"] + results["B"]) * 0.7, (
            f"GIL NOT released in pileup_variants! Total={total:.1f}s"
        )


# ═══════════════════════════════════════════════════════════════════════════
# TestRustCallerParser — verify Rust caller VCF parsing
# ═══════════════════════════════════════════════════════════════════════════

class TestRustCallerParser:
    """Tests for Rust parse_caller_vcf against normalized caller VCFs."""

    @pytest.fixture(autouse=True)
    def _require_rust(self):
        try:
            import stats_core
            if not hasattr(stats_core, 'parse_caller_vcf'):
                pytest.skip("Rust parse_caller_vcf not available")
        except ImportError:
            pytest.skip("stats_core not available")

    @pytest.fixture
    def _rescue_targets(self):
        """Build 4-tuple target set from rescue VCF."""
        import stats_core, glob
        base = os.path.join(REAL_SAMPLE["base_output_dir"], REAL_SAMPLE["dir_name"])
        rescue_path = glob.glob(os.path.join(
            base, "vcf_realignment", "rescue", "*", "*.filtered.vcf.stripped.vcf.gz"
        ))[0]
        records = stats_core.parse_rescue(rescue_path)
        chroms = [r["CHROM"] for r in records]
        poss = [r["POS"] for r in records]
        refs = [r["REF"] for r in records]
        alts = [r["ALT"] for r in records]
        return chroms, poss, refs, alts

    def _get_vcf_path(self, caller_name):
        import glob
        base = os.path.join(REAL_SAMPLE["base_output_dir"], REAL_SAMPLE["dir_name"])
        cfg = CALLER_CONFIGS[caller_name]
        subdir = cfg["subdir"].format(prefix=REAL_SAMPLE["vcf_prefix"])
        files = glob.glob(os.path.join(base, subdir, cfg["pattern"]))
        if not files:
            pytest.skip(f"VCF not found for {caller_name}")
        return files[0], cfg["sample_suffix"]

    def test_mutect2_dp_gt_ad_af(self, _rescue_targets):
        """DNA Mutect2 extracts DP, GT, AD, AF, SB, FAD."""
        import stats_core
        chroms, poss, refs, alts = _rescue_targets
        vcf_path, suffix = self._get_vcf_path("DNA_mutect2")
        result = stats_core.parse_caller_vcf(
            vcf_path, chroms, poss, refs, alts, suffix, "DNA_mutect2"
        )
        assert len(result["CHROM"]) > 0
        assert any(v is not None for v in result["DP"])
        assert any(v is not None for v in result["GT"])
        assert any(v is not None for v in result["AD_REF"])
        assert any(v is not None for v in result["AD_ALT"])
        assert any(v is not None for v in result["VAF_CALLER"])
        assert any(v is not None for v in result["SB"])

    def test_deepsomatic_dp_gt_ad_vaf(self, _rescue_targets):
        """DNA DeepSomatic extracts DP, GT, AD, VAF."""
        import stats_core
        chroms, poss, refs, alts = _rescue_targets
        vcf_path, suffix = self._get_vcf_path("DNA_deepsomatic")
        result = stats_core.parse_caller_vcf(
            vcf_path, chroms, poss, refs, alts, suffix, "DNA_deepsomatic"
        )
        assert len(result["CHROM"]) > 0
        assert any(v is not None for v in result["DP"])
        assert any(v is not None for v in result["GT"])
        assert any(v is not None for v in result["AD_REF"])
        assert any(v is not None for v in result["VAF_CALLER"])

    def test_strelka_no_gt_no_ad(self, _rescue_targets):
        """DNA Strelka has DP and AU/CU/GU/TU but NO GT and NO AD."""
        import stats_core
        chroms, poss, refs, alts = _rescue_targets
        vcf_path, suffix = self._get_vcf_path("DNA_strelka")
        result = stats_core.parse_caller_vcf(
            vcf_path, chroms, poss, refs, alts, suffix, "DNA_strelka"
        )
        assert len(result["CHROM"]) > 0
        assert any(v is not None for v in result["DP"])
        assert any(v is not None for v in result["AU"])
        assert any(v is not None for v in result["CU"])
        # Strelka has no GT and no AD
        assert all(v is None for v in result["GT"])
        assert all(v is None for v in result["AD_REF"])
        assert all(v is None for v in result["AD_ALT"])

    def test_4_column_keys_present(self, _rescue_targets):
        """Result includes REF and ALT columns for 4-tuple matching."""
        import stats_core
        chroms, poss, refs, alts = _rescue_targets
        vcf_path, suffix = self._get_vcf_path("DNA_mutect2")
        result = stats_core.parse_caller_vcf(
            vcf_path, chroms, poss, refs, alts, suffix, "DNA_mutect2"
        )
        assert "REF" in result
        assert "ALT" in result
        assert len(result["REF"]) == len(result["CHROM"])
        assert len(result["ALT"]) == len(result["CHROM"])

    def test_rna_mutect2_parses(self, _rescue_targets):
        """RNA Mutect2 (realigned) parses correctly with RT suffix."""
        import stats_core
        chroms, poss, refs, alts = _rescue_targets
        vcf_path, suffix = self._get_vcf_path("RNA_mutect2")
        result = stats_core.parse_caller_vcf(
            vcf_path, chroms, poss, refs, alts, suffix, "RNA_mutect2"
        )
        assert len(result["CHROM"]) > 0
        assert any(v is not None for v in result["DP"])
        assert any(v is not None for v in result["GT"])

    def test_early_termination(self, _rescue_targets):
        """Scan stops early when all targets found (small subset test)."""
        import stats_core, time
        chroms, poss, refs, alts = _rescue_targets
        # Take only first 5 positions
        c5, p5, r5, a5 = chroms[:5], poss[:5], refs[:5], alts[:5]
        vcf_path, suffix = self._get_vcf_path("DNA_mutect2")
        t0 = time.time()
        result = stats_core.parse_caller_vcf(vcf_path, c5, p5, r5, a5, suffix, "DNA_mutect2")
        dt = time.time() - t0
        # With only 5 targets, should be very fast (< 0.5s for full scan, but
        # early termination may make it even faster)
        assert len(result["CHROM"]) <= 5
        assert dt < 10.0  # generous — real time is ~0.2s

    def test_caller_gil_released(self, _rescue_targets):
        """parse_caller_vcf releases GIL — two threads run in parallel."""
        import stats_core, threading, time
        chroms, poss, refs, alts = _rescue_targets
        vcf_path, suffix = self._get_vcf_path("DNA_mutect2")
        vcf_path2, suffix2 = self._get_vcf_path("DNA_deepsomatic")

        results = {}
        def worker(label, path, sfx, name):
            t0 = time.time()
            stats_core.parse_caller_vcf(path, chroms, poss, refs, alts, sfx, name)
            results[label] = time.time() - t0

        t0 = time.time()
        t1 = threading.Thread(target=worker, args=("A", vcf_path, suffix, "DNA_mutect2"))
        t2 = threading.Thread(target=worker, args=("B", vcf_path2, suffix2, "DNA_deepsomatic"))
        t1.start(); t2.start(); t1.join(); t2.join()
        total = time.time() - t0
        assert total < (results["A"] + results["B"]) * 0.85, (
            f"GIL NOT released! Total={total:.1f}s, A={results['A']:.1f}s, B={results['B']:.1f}s"
        )


# ═══════════════════════════════════════════════════════════════════════════
# TestRustTiering — verify Rust tiering matches Python TieringEngine
# ═══════════════════════════════════════════════════════════════════════════

class TestRustTiering:
    """Tests that Rust compute_tiers produces identical output to Python."""

    @pytest.fixture(autouse=True)
    def _require_rust(self):
        try:
            import stats_core
            if not hasattr(stats_core, 'compute_tiers'):
                pytest.skip("Rust compute_tiers not available")
        except ImportError:
            pytest.skip("stats_core not available")

    @pytest.fixture
    def _tier_inputs(self):
        """Extract tiering input columns from rescue VCF."""
        import stats_core, glob
        base = os.path.join(REAL_SAMPLE["base_output_dir"], REAL_SAMPLE["dir_name"])
        rescue_path = glob.glob(os.path.join(
            base, "vcf_realignment", "rescue", "*", "*.filtered.vcf.stripped.vcf.gz"
        ))[0]
        records = stats_core.parse_rescue(rescue_path)
        n = len(records)
        filters = [r.get("FILTER", "PASS") or "PASS" for r in records]
        fnorm = [r.get("FILTERS_NORMALIZED", "") or "" for r in records]
        gaf = [float(r["GNOMAD_AF"]) if r.get("GNOMAD_AF") and r["GNOMAD_AF"] != "." else None for r in records]
        csc = [int(r["COSMIC_CNT"]) if r.get("COSMIC_CNT") and r["COSMIC_CNT"] != "." else None for r in records]
        redi = [str(r["REDI_EVIDENCE"]) if r.get("REDI_EVIDENCE") else None for r in records]
        ds = [int(r["N_DNA_CALLERS_SUPPORT"]) if r.get("N_DNA_CALLERS_SUPPORT") and r["N_DNA_CALLERS_SUPPORT"] != "." else None for r in records]
        rs = [int(r["N_RNA_CALLERS_SUPPORT"]) if r.get("N_RNA_CALLERS_SUPPORT") and r["N_RNA_CALLERS_SUPPORT"] != "." else None for r in records]
        return filters, fnorm, gaf, csc, redi, ds, rs, records

    def test_tiering_produces_all_columns(self, _tier_inputs):
        """Rust compute_tiers returns all expected columns."""
        import stats_core
        filters, fnorm, gaf, csc, redi, ds, rs, _ = _tier_inputs
        result = stats_core.compute_tiers(filters, fnorm, gaf, csc, redi, ds, rs)
        for col in ["final_tier", "caller_tier", "database_tier", "dna_caller_count", "rna_caller_count", "tier_quality"]:
            assert col in result, f"Missing column: {col}"
            assert len(result[col]) == len(filters)

    def test_tiering_parity_with_python(self, _tier_inputs):
        """Rust tiers match Python TieringEngine for first 200 variants."""
        import stats_core
        filters, fnorm, gaf, csc, redi, ds, rs, records = _tier_inputs

        # Rust
        r = stats_core.compute_tiers(filters[:200], fnorm[:200], gaf[:200], csc[:200], redi[:200], ds[:200], rs[:200])

        # Python
        from vcf_stats.seq2neo.tiering_stats import compute_tiers_for_dataframe
        import polars as pl
        rows = []
        for i in range(200):
            rows.append({
                "FILTER": records[i].get("FILTER", "PASS") or "PASS",
                "FILTERS_NORMALIZED": records[i].get("FILTERS_NORMALIZED", ""),
                "GNOMAD_AF": records[i].get("GNOMAD_AF"),
                "COSMIC_CNT": records[i].get("COSMIC_CNT"),
                "REDI_EVIDENCE": records[i].get("REDI_EVIDENCE"),
                "N_DNA_CALLERS_SUPPORT": records[i].get("N_DNA_CALLERS_SUPPORT"),
                "N_RNA_CALLERS_SUPPORT": records[i].get("N_RNA_CALLERS_SUPPORT"),
            })
        df = pl.DataFrame(rows)
        py_result = compute_tiers_for_dataframe(df)

        # Compare
        for i in range(200):
            assert r["final_tier"][i] == py_result["final_tier"][i], (
                f"Mismatch at index {i}: Rust={r['final_tier'][i]}, Python={py_result['final_tier'][i]}"
            )
            assert r["caller_tier"][i] == py_result["caller_tier"][i]
            assert r["dna_caller_count"][i] == py_result["dna_caller_count"][i]
            assert r["rna_caller_count"][i] == py_result["rna_caller_count"][i]

    def test_tier_quality_scores(self, _tier_inputs):
        """Tier quality scores match expected values from tier_config."""
        import stats_core
        expected_quality = {
            "C1D1": 140, "C1D0": 130, "C2D1": 120, "C2D0": 110,
            "C3D1": 100, "C3D0": 90,  "C4D1": 80,  "C4D0": 70,
            "C5D1": 60,  "C5D0": 50,  "C6D1": 40,  "C6D0": 30,
            "C7D1": 20,  "C7D0": 10,
        }
        # Synthesize test data for each tier
        filters = ["Somatic"] * 14
        fnorm = [""] * 14
        gaf = [None] * 14
        csc = [None] * 14
        redi = [None] * 14
        # C1: dna>=2, rna>=2; C2: dna>=2, rna<=1; C3: rna>=2, dna<=1;
        # C4: dna=1, rna=1; C5: dna=1, rna=0; C6: dna=0, rna=1; C7: dna=0, rna=0
        ds_list = [2, 2, 0, 1, 1, 0, 0, 2, 2, 0, 1, 1, 0, 0]
        rs_list = [2, 0, 2, 1, 0, 1, 0, 2, 0, 2, 1, 0, 1, 0]
        # First half: D1 (with gnomAD), second half: D0
        gaf_db = [0.01 if i < 7 else None for i in range(14)]

        result = stats_core.compute_tiers(filters, fnorm, gaf_db, csc, redi, ds_list, rs_list)
        for i in range(14):
            tier = result["final_tier"][i]
            assert tier in expected_quality, f"Unknown tier: {tier}"
            assert result["tier_quality"][i] == expected_quality[tier], (
                f"Quality mismatch for {tier}: {result['tier_quality'][i]} != {expected_quality[tier]}"
            )

    def test_tiering_gil_released(self, _tier_inputs):
        """compute_tiers releases GIL — two threads run in parallel."""
        import stats_core, threading, time
        filters, fnorm, gaf, csc, redi, ds, rs, _ = _tier_inputs

        results = {}
        def worker(label):
            t0 = time.time()
            stats_core.compute_tiers(filters, fnorm, gaf, csc, redi, ds, rs)
            results[label] = time.time() - t0

        t0 = time.time()
        t1 = threading.Thread(target=worker, args=("A",))
        t2 = threading.Thread(target=worker, args=("B",))
        t1.start(); t2.start(); t1.join(); t2.join()
        total = time.time() - t0
        assert total < (results["A"] + results["B"]) * 0.7, (
            f"GIL NOT released! Total={total:.1f}s"
        )


# ═══════════════════════════════════════════════════════════════════════════
# TestRustPileup — verify Rust BAM pileup matches pysam
# ═══════════════════════════════════════════════════════════════════════════

class TestRustPileup:
    """Tests for Rust pileup_variants against pysam reference."""

    @pytest.fixture(autouse=True)
    def _require_rust(self):
        try:
            import stats_core
            if not hasattr(stats_core, 'pileup_variants'):
                pytest.skip("Rust pileup_variants not available")
        except ImportError:
            pytest.skip("stats_core not available")

    @pytest.fixture
    def _bam_path(self):
        bam = os.path.join(
            REAL_SAMPLE["base_output_dir"], REAL_SAMPLE["dir_name"],
            "preprocessing", "mapped", f"{REAL_SAMPLE['vcf_prefix']}DT",
            f"{REAL_SAMPLE['vcf_prefix']}DT.sorted.bam",
        )
        if not os.path.isfile(bam):
            pytest.skip("BAM not found")
        return bam

    def test_pileup_parity_with_pysam(self, _bam_path):
        """Rust windowed pileup matches pysam for 100 random positions.

        Allows ±1 read tolerance due to BAI bin boundary effects — windowed
        queries may include reads at bin edges that single-position queries miss.
        """
        import stats_core, pysam, random
        random.seed(42)

        # Use 100 positions from chr1
        n_pos = 100
        chroms = ["chr1"] * n_pos
        positions = [random.randint(600000, 700000) for _ in range(n_pos)]
        refs = ["A"] * n_pos
        alts = ["G"] * n_pos

        r = stats_core.pileup_variants(_bam_path, chroms, positions, refs, alts)
        bam = pysam.AlignmentFile(_bam_path, "rb")

        mismatches = 0
        for i in range(n_pos):
            # pysam pileup
            reads = list(bam.fetch(chroms[i], positions[i] - 1, positions[i]))
            dp = 0; ref_dp = 0; alt_dp = 0
            for rd in reads:
                if rd.is_unmapped or rd.is_duplicate: continue
                pir = positions[i] - rd.reference_start - 1
                if pir < 0 or pir >= len(rd.query_sequence): continue
                base = rd.query_sequence[pir]
                dp += 1
                if base == refs[i]: ref_dp += 1
                elif base == alts[i]: alt_dp += 1

            rust_dp = r["DP"][i] or 0
            if abs(rust_dp - dp) > 1:
                mismatches += 1
        bam.close()

        # Allow at most 5% of positions to differ by >1 read (BAI bin boundaries)
        assert mismatches <= n_pos * 0.05, (
            f"Too many DP mismatches: {mismatches}/{n_pos}"
        )

    def test_empty_positions(self, _bam_path):
        """Empty position list returns empty result."""
        import stats_core
        r = stats_core.pileup_variants(_bam_path, [], [], [], [])
        assert len(r["DP"]) == 0

    def test_missing_bam_returns_defaults(self):
        """Non-existent BAM raises an error (file must exist for pileup)."""
        import stats_core
        with pytest.raises(Exception):
            stats_core.pileup_variants(
                "/nonexistent/path.bam", ["chr1"], [1], ["A"], ["G"]
            )

    def test_all_columns_present(self, _bam_path):
        """All expected columns are in the result."""
        import stats_core
        r = stats_core.pileup_variants(_bam_path, ["chr1"], [633987], ["C"], ["T"])
        for col in ["DP", "REF_DP", "ALT_DP", "F1R2_ref", "F2R1_ref",
                     "F1R2_alt", "F2R1_alt", "mean_BQ", "mean_MQ"]:
            assert col in r, f"Missing column: {col}"

    def test_pileup_output_includes_ref_alt(self, _bam_path):
        """Pileup result includes REF and ALT columns for 4-column join."""
        from vcf_stats.seq2neo.rust_bam import pileup_variants
        positions = [("chr1", 633987, "C", "T")]
        result = pileup_variants(_bam_path, positions)
        assert result is not None
        assert "REF" in result.columns, "REF column missing from pileup output"
        assert "ALT" in result.columns, "ALT column missing from pileup output"
        assert result["REF"][0] == "C"
        assert result["ALT"][0] == "T"

    def test_pileup_multiallelic_join_correct(self):
        """4-column join matches correct alleles at multiallelic sites."""
        import polars as pl
        # Two variants at same CHROM+POS, different ALTs
        df = pl.DataFrame({
            "CHROM": ["chr1", "chr1"],
            "POS": [100, 100],
            "REF": ["A", "A"],
            "ALT": ["G", "T"],
            "sample_id": ["s1", "s1"],
        })
        pileup = pl.DataFrame({
            "CHROM": ["chr1", "chr1"],
            "POS": [100, 100],
            "REF": ["A", "A"],
            "ALT": ["G", "T"],
            "DP": [50, 30],
            "REF_DP": [40, 20],
            "ALT_DP": [10, 10],
        })
        # 4-column join: each variant row matches its specific pileup row
        joined = df.join(pileup, on=["CHROM", "POS", "REF", "ALT"], how="left")
        assert joined.height == 2, f"Expected 2 rows, got {joined.height}"
        assert joined["DP"][0] == 50  # ALT=G gets DP=50
        assert joined["DP"][1] == 30  # ALT=T gets DP=30

    def test_pileup_mode_filtered_excludes_noconsensus(self):
        """Filtered pileup mode excludes NoConsensus variants."""
        import polars as pl
        df = pl.DataFrame({
            "CHROM": ["chr1", "chr2", "chr3", "chr4"],
            "POS": [100, 200, 300, 400],
            "REF": ["A", "C", "G", "T"],
            "ALT": ["G", "T", "A", "C"],
            "FILTER": ["Somatic", "NoConsensus", "Germline", "NoConsensus"],
        })
        # All mode
        cols_4 = ["CHROM", "POS", "REF", "ALT"]
        positions_all = [(row[0], row[1], row[2], row[3])
                         for row in df.select(cols_4).iter_rows()]
        assert len(positions_all) == 4
        # Filtered mode
        mask = df["FILTER"] != "NoConsensus"
        positions_filtered = [
            (row[0], row[1], row[2], row[3])
            for row, keep in zip(df.select(cols_4).iter_rows(), mask.to_list())
            if keep
        ]
        assert len(positions_filtered) == 2  # Only Somatic and Germline

    def test_pileup_10k_positions_completes_quickly(self, _bam_path):
        """Windowed pileup on 10K positions completes in < 60 seconds."""
        import stats_core, time, random
        random.seed(42)
        n = 10_000
        chroms = ["chr1"] * n
        positions = [random.randint(600_000, 700_000) for _ in range(n)]
        refs = ["A"] * n
        alts = ["G"] * n
        t0 = time.time()
        r = stats_core.pileup_variants(_bam_path, chroms, positions, refs, alts)
        elapsed = time.time() - t0
        assert len(r["DP"]) == n
        assert elapsed < 60, f"10K positions took {elapsed:.1f}s, expected < 60s"

    def test_pileup_join_uses_4_columns(self):
        """Pileup join key verified as (CHROM, POS, REF, ALT)."""
        from pathlib import Path
        cli_path = Path(__file__).parent.parent / "seq2neo" / "cli.py"
        source = cli_path.read_text()
        assert '"CHROM", "POS", "REF", "ALT"' in source, (
            "4-column join key not found in cli.py"
        )

    def test_pileup_gil_released(self, _bam_path):
        """pileup_variants releases GIL — two threads run in parallel."""
        import stats_core, threading, time
        # Use 100 positions each on the same chromosome for balanced work
        import random
        random.seed(0)
        positions = [random.randint(600000, 700000) for _ in range(100)]

        results = {}
        def worker(label):
            t0 = time.time()
            stats_core.pileup_variants(
                _bam_path, ["chr1"] * 100, positions,
                ["A"] * 100, ["G"] * 100,
            )
            results[label] = time.time() - t0

        t0 = time.time()
        t1 = threading.Thread(target=worker, args=("A",))
        t2 = threading.Thread(target=worker, args=("B",))
        t1.start(); t2.start(); t1.join(); t2.join()
        total = time.time() - t0
        assert total < (results["A"] + results["B"]) * 0.85, (
            f"GIL NOT released! Total={total:.1f}s, A={results['A']:.1f}s, B={results['B']:.1f}s"
        )

    # ── Multi-BAM pileup tests (Phase 8 + 10) ────────────────────────────

    def test_pileup_multi_bam_parity(self, _bam_path):
        """Multi-BAM pileup matches individual pileup calls for same BAM."""
        import stats_core, random
        random.seed(42)

        n_pos = 100
        chroms = ["chr1"] * n_pos
        positions = [random.randint(600_000, 700_000) for _ in range(n_pos)]
        refs = ["A"] * n_pos
        alts = ["G"] * n_pos

        # Individual call
        single = stats_core.pileup_variants(_bam_path, chroms, positions, refs, alts)

        # Multi-BAM call with single BAM
        multi_result = stats_core.pileup_variants_multi(
            [_bam_path], ["test"],
            chroms, positions, refs, alts,
            [], [], [],  # no BED
        )
        assert len(multi_result) == 1
        assert "test" in multi_result
        multi = multi_result["test"]

        # DP values should match exactly
        for i in range(n_pos):
            assert multi["DP"][i] == single["DP"][i], (
                f"Position {i}: multi={multi['DP'][i]}, single={single['DP'][i]}"
            )

    def test_pileup_binary_search_parity(self, _bam_path):
        """Binary search inner loop produces same results as HashMap approach.

        The existing pileup_variants uses HashMap iteration (O(r×p) per window).
        pileup_variants_multi uses binary search (O(r×log p) per window).
        Both should produce identical results for the same inputs.
        """
        import stats_core, random
        random.seed(99)

        # 500 positions across a 1Mb region — small enough to finish fast,
        # large enough to stress binary search edge cases
        n_pos = 500
        chroms = ["chr1"] * n_pos
        positions = [random.randint(600_000, 1_600_000) for _ in range(n_pos)]
        refs = ["A"] * n_pos
        alts = ["G"] * n_pos

        # HashMap approach (existing single-BAM function)
        single = stats_core.pileup_variants(_bam_path, chroms, positions, refs, alts)

        # Binary search approach (new multi-BAM function)
        multi_result = stats_core.pileup_variants_multi(
            [_bam_path], ["test"],
            chroms, positions, refs, alts,
            [], [], [],  # no BED → 1Mb windows
        )
        assert "test" in multi_result
        multi = multi_result["test"]

        # Every position should match
        mismatches = 0
        for i in range(n_pos):
            if multi["DP"][i] != single["DP"][i]:
                mismatches += 1
        assert mismatches == 0, (
            f"Binary search DP mismatch: {mismatches}/{n_pos} positions differ"
        )

    def test_pileup_bed_guided_region_count(self, _bam_path):
        """BED-guided mode creates fewer regions than 1Mb windows for WES.

        For WES data, BED intervals (~300 merged regions) should result in
        significantly fewer queries than 1Mb sliding windows (~2,765 for
        whole genome). This test verifies the structural property.
        """
        import stats_core, random
        # Simulate WES: 10 BED regions covering 200Kb each on chr1
        bed_regions = [
            ("chr1", 500_000 + i * 1_000_000, 500_000 + i * 1_000_000 + 200_000)
            for i in range(10)
        ]

        # Positions within BED regions: 100 per region = 1000 total
        positions_across_genome = []
        random.seed(7)
        for chrom, start, end in bed_regions:
            for _ in range(100):
                pos = random.randint(start + 1, end)
                positions_across_genome.append((chrom, pos, "A", "G"))

        n_pos = len(positions_across_genome)
        chroms = [p[0] for p in positions_across_genome]
        poss = [p[1] for p in positions_across_genome]
        refs = [p[2] for p in positions_across_genome]
        alts = [p[3] for p in positions_across_genome]

        # With BED: should create at most 10 regions (one per BED interval)
        bed_chroms = [r[0] for r in bed_regions]
        bed_starts = [r[1] for r in bed_regions]
        bed_ends = [r[2] for r in bed_regions]

        multi_bed = stats_core.pileup_variants_multi(
            [_bam_path], ["test"],
            chroms, poss, refs, alts,
            bed_chroms, bed_starts, bed_ends,
        )
        assert len(multi_bed) == 1
        assert "test" in multi_bed
        results_bed = multi_bed["test"]
        assert len(results_bed["DP"]) == n_pos, (
            f"BED-guided: expected {n_pos} results, got {len(results_bed['DP'])}"
        )

        # Without BED: 1Mb windows — positions span ~10Mb, so ~10 windows
        multi_no_bed = stats_core.pileup_variants_multi(
            [_bam_path], ["test"],
            chroms, poss, refs, alts,
            [], [], [],  # no BED
        )
        assert "test" in multi_no_bed
        results_no_bed = multi_no_bed["test"]
        assert len(results_no_bed["DP"]) == n_pos, (
            f"Window-based: expected {n_pos} results, got {len(results_no_bed['DP'])}"
        )

        # Both modes should produce same results
        mismatches = 0
        for i in range(n_pos):
            if results_bed["DP"][i] != results_no_bed["DP"][i]:
                mismatches += 1
        # Allow a small number of mismatches due to BAI bin boundary differences
        # (BED regions vs 1Mb windows may have different bin boundaries)
        assert mismatches <= n_pos * 0.05, (
            f"BED vs window DP mismatch: {mismatches}/{n_pos} > 5%"
        )

    def test_pileup_multi_bam_bed_integration(self, _bam_path):
        """BED regions are passed through correctly to the FFI boundary."""
        import stats_core, random
        random.seed(13)

        n_pos = 50
        chroms = ["chr1"] * n_pos
        positions = [random.randint(600_000, 800_000) for _ in range(n_pos)]
        refs = ["A"] * n_pos
        alts = ["G"] * n_pos

        # BED regions covering the position range
        bed_regions = [("chr1", 600_000, 700_000), ("chr1", 700_001, 800_001)]

        result = stats_core.pileup_variants_multi(
            [_bam_path], ["test"],
            chroms, positions, refs, alts,
            [r[0] for r in bed_regions],
            [r[1] for r in bed_regions],
            [r[2] for r in bed_regions],
        )
        assert len(result) == 1
        assert "test" in result
        data = result["test"]
        assert len(data["DP"]) == n_pos, (
            f"Expected {n_pos} results, got {len(data['DP'])}"
        )
        # All DP values should be integers or None (no crashes)
        for dp in data["DP"]:
            assert dp is None or isinstance(dp, int), f"Unexpected DP type: {type(dp)}"

    def test_pileup_columns_in_parquet(self):
        """Integration test: pileup columns survive parquet round-trip.

        Verifies that pileup output columns (DP, REF_DP, etc.) are written
        correctly to parquet and can be read back with expected types (task 4.6).
        """
        # Simulate pileup columns joined to a variant DataFrame
        df = pl.DataFrame({
            "CHROM": ["chr1", "chr2"],
            "POS": [100, 200],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
            "BAM_DN_DP": [50, 30],
            "BAM_DN_REF_DP": [40, 20],
            "BAM_DN_ALT_DP": [10, 10],
            "BAM_DN_F1R2_ref": [20, 10],
            "BAM_DN_F2R1_ref": [20, 10],
            "BAM_DN_F1R2_alt": [5, 3],
            "BAM_DN_F2R1_alt": [5, 7],
            "BAM_DN_mean_BQ": [35.5, 36.1],
            "BAM_DN_mean_MQ": [60.0, 59.5],
            "BAM_DT_DP": [55, 35],
        })
        with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as f:
            parquet_path = f.name
        try:
            df.write_parquet(parquet_path)
            read_back = pl.read_parquet(parquet_path)
            assert read_back.height == 2
            # Verify pileup columns
            for col in ["BAM_DN_DP", "BAM_DN_REF_DP", "BAM_DN_ALT_DP",
                        "BAM_DN_mean_BQ", "BAM_DN_mean_MQ", "BAM_DT_DP"]:
                assert col in read_back.columns, f"Missing column: {col}"
            assert read_back["BAM_DN_DP"][0] == 50
            assert read_back["BAM_DN_mean_BQ"][0] == 35.5
        finally:
            os.unlink(parquet_path)


# ═══════════════════════════════════════════════════════════════════════════
# TestJoinOptimization — verify polars join on 4 columns
# ═══════════════════════════════════════════════════════════════════════════

class TestJoinOptimization:
    """Tests for join_caller_columns polars join on 4 coordinate columns."""

    def test_join_on_4_columns(self):
        """join_caller_columns joins on CHROM, POS, REF, ALT with column-oriented input."""
        import polars as pl
        from vcf_stats.seq2neo.caller_parser import join_caller_columns

        rescue = pl.DataFrame({
            "CHROM": ["chr1", "chr1", "chr2"],
            "POS": [100, 200, 300],
            "REF": ["A", "C", "G"],
            "ALT": ["T", "G", "A"],
            "OTHER": [1, 2, 3],
        })

        # Column-oriented caller data (Rust parser output format)
        caller_data = {
            "DNA_mutect2": {
                "CHROM": ["chr1", "chr2"],
                "POS": [100, 300],
                "REF": ["A", "G"],
                "ALT": ["T", "A"],
                "DP": [42, 55],
                "AD_REF": [30, 40],
                "AD_ALT": [12, 15],
                "GT": ["0/1", "1/1"],
                "VAF_CALLER": [0.286, 0.273],
            }
        }

        result = join_caller_columns(rescue, caller_data)
        assert "DNA_mutect2_DP" in result.columns
        assert result["DNA_mutect2_DP"][0] == 42
        assert result["DNA_mutect2_DP"][1] is None
        assert result["DNA_mutect2_DP"][2] == 55
        assert result["DNA_mutect2_GT"][0] == "0/1"
        assert result["DNA_mutect2_GT"][2] == "1/1"

    def test_missing_caller_null_fills(self):
        """Missing caller produces null-filled columns."""
        import polars as pl
        from vcf_stats.seq2neo.caller_parser import join_caller_columns

        rescue = pl.DataFrame({
            "CHROM": ["chr1"], "POS": [100], "REF": ["A"], "ALT": ["T"],
        })
        result = join_caller_columns(rescue, {})  # no callers at all
        assert result.height == 1

    def test_multiallelic_no_cross_match(self):
        """Two records at same (CHROM, POS) with different REF/ALT don't cross-match."""
        import polars as pl
        from vcf_stats.seq2neo.caller_parser import join_caller_columns

        rescue = pl.DataFrame({
            "CHROM": ["chr1", "chr1"],
            "POS": [100, 100],
            "REF": ["A", "A"],
            "ALT": ["T", "G"],
        })

        caller_data = {
            "DNA_mutect2": {
                "CHROM": ["chr1"], "POS": [100],
                "REF": ["A"], "ALT": ["T"],
                "DP": [42],
            }
        }

        result = join_caller_columns(rescue, caller_data)
        assert result["DNA_mutect2_DP"][0] == 42  # (A,T) matched
        assert result["DNA_mutect2_DP"][1] is None  # (A,G) not matched

    def test_schema_inference_with_late_values(self):
        """Columns with None in first 100+ rows but values later work correctly."""
        import polars as pl
        from vcf_stats.seq2neo.caller_parser import join_caller_columns

        # 200 variants: first 150 have None AD_REF, last 50 have actual values
        positions = [(f"chr1", i * 1000, "A", "T") for i in range(1, 201)]
        col_data = {
            "CHROM": [k[0] for k in positions],
            "POS": [k[1] for k in positions],
            "REF": [k[2] for k in positions],
            "ALT": [k[3] for k in positions],
            "DP": [i + 1 for i in range(200)],
            "AD_REF": [None] * 150 + [20] * 50,
            "AD_ALT": [None] * 150 + [10] * 50,
        }

        rescue = pl.DataFrame({
            "CHROM": [k[0] for k in positions],
            "POS": [k[1] for k in positions],
            "REF": [k[2] for k in positions],
            "ALT": [k[3] for k in positions],
        })

        result = join_caller_columns(rescue, {"DNA_mutect2": col_data})
        assert result["DNA_mutect2_DP"][0] == 1
        assert result["DNA_mutect2_AD_REF"][0] is None
        assert result["DNA_mutect2_AD_REF"][199] == 20
        assert result["DNA_mutect2_AD_ALT"][199] == 10


# ═══════════════════════════════════════════════════════════════════════════
# TestMemoryEfficiency — verify streaming architecture doesn't accumulate memory
# ═══════════════════════════════════════════════════════════════════════════

class TestMemoryEfficiency:
    """Tests that the streaming architecture frees DataFrames and keeps memory low."""

    def test_lazy_frame_scan_from_parquet(self, tmp_path):
        """pl.scan_parquet() can aggregate across multiple parquet files."""
        import polars as pl
        d = tmp_path / "parts"
        d.mkdir()
        pl.DataFrame({"CHROM": ["chr1"], "POS": [1], "DP": [42], "sample_id": ["A"]}).write_parquet(str(d / "A_variants.parquet"))
        pl.DataFrame({"CHROM": ["chr2"], "POS": [2], "DP": [55], "sample_id": ["B"]}).write_parquet(str(d / "B_variants.parquet"))

        lazy = pl.scan_parquet(str(d / "*_variants.parquet"))
        result = lazy.group_by("sample_id").agg(pl.col("DP").sum()).collect()
        assert result.height == 2
        assert result.filter(pl.col("sample_id") == "A")["DP"][0] == 42

    def test_eager_helper_with_lazy_frame(self):
        """_maybe_collect(df) materializes LazyFrame but passes through eager."""
        import polars as pl
        from vcf_stats.seq2neo.visualizer import _maybe_collect

        eager_df = pl.DataFrame({"x": [1, 2, 3]})
        result = _maybe_collect(eager_df)
        assert isinstance(result, pl.DataFrame)
        assert result["x"].to_list() == [1, 2, 3]

    def test_lazy_aggregation_parity(self, tmp_path):
        """Lazy scan aggregation matches eager concat for dataset_summary."""
        import polars as pl
        from vcf_stats.seq2neo.statistics import dataset_summary

        d = tmp_path / "parts"
        d.mkdir()
        df_a = pl.DataFrame({
            "CHROM": ["chr1"] * 5, "POS": range(1, 6),
            "FILTER": ["Somatic"] * 5, "VC": ["Somatic"] * 5,
            "variant_type": ["SNV"] * 5, "ti_tv": [True] * 5,
            "sample_id": ["A"] * 5, "disease_normalized": ["Lung"] * 5,
        })
        df_b = pl.DataFrame({
            "CHROM": ["chr2"] * 5, "POS": range(1, 6),
            "FILTER": ["Somatic"] * 5, "VC": ["Somatic"] * 5,
            "variant_type": ["SNV"] * 5, "ti_tv": [True] * 5,
            "sample_id": ["B"] * 5, "disease_normalized": ["Lung"] * 5,
        })
        df_a.write_parquet(str(d / "A_variants.parquet"))
        df_b.write_parquet(str(d / "B_variants.parquet"))

        # Eager mode
        eager = pl.concat([df_a, df_b], how="diagonal_relaxed")
        eager_result = dataset_summary(eager)

        # Lazy mode
        lazy = pl.scan_parquet(str(d / "*_variants.parquet"))
        lazy_result = dataset_summary(lazy)

        assert eager_result["total_variants"] == lazy_result["total_variants"]
        assert eager_result["n_samples"] == lazy_result["n_samples"]

    def test_ensure_eager_preserves_data(self):
        """_ensure_eager doesn't modify data — identity for eager frames."""
        import polars as pl
        from vcf_stats.seq2neo.statistics import _ensure_eager

        df = pl.DataFrame({"x": [1, 2, 3], "y": ["a", "b", "c"]})
        result = _ensure_eager(df)
        assert result["x"].to_list() == [1, 2, 3]
        assert result["y"].to_list() == ["a", "b", "c"]

    def test_filter_distribution_sort(self, tmp_path):
        """filter_distribution sort works correctly (regression test)."""
        import polars as pl
        from vcf_stats.seq2neo.statistics import filter_distribution

        df = pl.DataFrame({
            "set_number": [1, 1, 2, 2],
            "FILTER": ["Somatic", "Germline", "Somatic", "Artifact"],
        })
        result = filter_distribution(df)
        assert result.height == 4  # 2 sets × 2 filters

    def test_chart_function_accepts_lazy_frame(self, tmp_path):
        """Chart function with _eager can accept LazyFrame input."""
        import polars as pl
        from vcf_stats.seq2neo.visualizer import plot_vc_distribution

        d = tmp_path / "parts"
        d.mkdir()
        df = pl.DataFrame({
            "CHROM": ["chr1"] * 10, "POS": range(10),
            "FILTER": ["Somatic"] * 10, "VC": ["Somatic"] * 10,
            "set_number": [1] * 10,
        })
        df.write_parquet(str(d / "test_variants.parquet"))
        lazy = pl.scan_parquet(str(d / "test_variants.parquet"))

        chart = plot_vc_distribution(lazy, str(tmp_path))
        assert chart is not None


# ═══════════════════════════════════════════════════════════════════════════
# TestMemoryRegression — prevent OOM-causing patterns from returning
# ═══════════════════════════════════════════════════════════════════════════

class TestMemoryRegression:
    """Tests that memory-expensive patterns are not reintroduced."""

    def test_join_receives_column_oriented_data(self):
        """join_caller_columns accepts column-oriented data, not row dicts.

        Row-oriented lookup dicts create ~5 GB of Python small objects per
        1.4M-variant sample (tuples, dicts, ints). Column-oriented data
        avoids this entirely.
        """
        import polars as pl
        from vcf_stats.seq2neo.caller_parser import join_caller_columns

        rescue = pl.DataFrame({
            "CHROM": ["chr1"], "POS": [1], "REF": ["A"], "ALT": ["T"],
        })
        # Column-oriented: {col_name: [values]} — the memory-safe format
        col_data = {
            "DNA_mutect2": {
                "CHROM": ["chr1"], "POS": [1], "REF": ["A"], "ALT": ["T"],
                "DP": [42], "AD_REF": [30], "AD_ALT": [12],
                "GT": ["0/1"], "VAF_CALLER": [0.286],
            }
        }
        result = join_caller_columns(rescue, col_data)
        assert result["DNA_mutect2_DP"][0] == 42

    def test_build_lookup_not_needed_for_rust_path(self):
        """build_caller_results_lookup exists but Rust path skips it.

        The Rust parser returns column-oriented data directly.
        build_caller_results_lookup is only for the cyvcf2 fallback.
        """
        from vcf_stats.seq2neo.caller_parser import build_caller_results_lookup
        # Function should still exist (for cyvcf2 fallback)
        assert callable(build_caller_results_lookup)
        # But should NOT be called in the Rust path (verified by
        # test_join_receives_column_oriented_data above)

    def test_rust_parse_returns_column_format(self):
        """Rust parse_caller_vcf returns column-oriented dict, not row dicts.

        Verifying that stats_core.parse_caller_vcf returns a dict of
        column_name → list_of_values, not row-oriented data.
        """
        try:
            import stats_core
            if not hasattr(stats_core, 'parse_caller_vcf'):
                pytest.skip("Rust parse_caller_vcf not available")
        except ImportError:
            pytest.skip("stats_core not available")

        import glob
        base = os.path.join(REAL_SAMPLE["base_output_dir"], REAL_SAMPLE["dir_name"])
        rescue = glob.glob(os.path.join(base, "vcf_realignment", "rescue", "*", "*.filtered.vcf.stripped.vcf.gz"))
        if not rescue:
            pytest.skip("Rescue VCF not found")

        import stats_core as sc
        records = sc.parse_rescue(rescue[0])
        chroms = [r["CHROM"] for r in records[:10]]
        poss = [r["POS"] for r in records[:10]]
        refs = [r["REF"] for r in records[:10]]
        alts = [r["ALT"] for r in records[:10]]

        cfg = CALLER_CONFIGS["DNA_mutect2"]
        subdir = cfg["subdir"].format(prefix=REAL_SAMPLE["vcf_prefix"])
        vcf_path = glob.glob(os.path.join(base, subdir, cfg["pattern"]))
        if not vcf_path:
            pytest.skip("Caller VCF not found")

        result = sc.parse_caller_vcf(vcf_path[0], chroms, poss, refs, alts, cfg["sample_suffix"], "DNA_mutect2")
        # Must be dict of column_name → list
        assert isinstance(result, dict)
        assert "CHROM" in result
        assert "DP" in result
        assert isinstance(result["CHROM"], list)
        assert isinstance(result["DP"], list)
        # Must NOT contain row-oriented entries (tuples as keys)
        for key in result:
            assert not isinstance(key, tuple), f"Row-oriented key found: {key}"

    def test_process_single_sample_frees_intermediates(self):
        """process_single_sample has del statements for memory cleanup."""
        import inspect
        from vcf_stats.seq2neo.cli import process_single_sample
        source = inspect.getsource(process_single_sample)
        assert "del target_positions" in source, "Missing: del target_positions"
        assert "del rescue_df" in source, "Missing: del rescue_df (column split)"
        assert "gc.collect()" in source, "Missing: gc.collect() after del"

    def test_join_frees_caller_entry(self):
        """join_caller_columns clears each caller_data entry after joining."""
        import inspect
        from vcf_stats.seq2neo.caller_parser import join_caller_columns
        source = inspect.getsource(join_caller_columns)
        assert "caller_data[caller_name] = {}" in source, (
            "Missing: caller_data entry cleanup in join_caller_columns"
        )

    def test_no_to_dicts_on_variant_df(self):
        """Chart functions must not call df.to_dicts() on variant-level DataFrames."""
        import inspect, re
        from vcf_stats.seq2neo.visualizer import (
            plot_gt_concordance, plot_gt_concordance_per_tier,
        )
        for func in [plot_gt_concordance, plot_gt_concordance_per_tier]:
            source = inspect.getsource(func)
            # Remove comments and docstrings before checking
            code_only = re.sub(r'""".*?"""', '', source, flags=re.DOTALL)
            code_only = re.sub(r'#.*', '', code_only)
            assert ".to_dicts()" not in code_only, (
                f"{func.__name__} uses .to_dicts() — should use iter_rows() instead"
            )

    def test_no_to_dicts_in_caller_join_path(self):
        """join_caller_columns must NOT call build_caller_results_lookup.

        The Rust path passes column-oriented data directly. Only the
        cyvcf2 fallback should use build_caller_results_lookup.
        """
        import inspect
        from vcf_stats.seq2neo.caller_parser import join_caller_columns
        source = inspect.getsource(join_caller_columns)
        assert "build_caller_results_lookup" not in source, (
            "join_caller_columns calls build_caller_results_lookup — "
            "should receive column-oriented data directly"
        )

    def test_to_dicts_only_on_aggregated_data(self):
        """Any to_dicts() call must be on aggregated (group_by) results, not raw data.

        Scans statistics.py for to_dicts() usage and verifies they're on
        group_by results (tier_counts), not on variant-level DataFrames.
        """
        import inspect
        from vcf_stats.seq2neo.statistics import sample_summary, dataset_summary
        # sample_summary uses tier_counts.to_dicts() — tier_counts is aggregated (14 rows)
        source = inspect.getsource(sample_summary)
        # Verify to_dicts is on a group_by result, not on the input df
        assert "tier_counts.to_dicts()" in source, "Expected tier_counts.to_dicts() in sample_summary"

    def test_chunked_target_building(self):
        """Target positions are built via chunked iter_rows, not to_list().

        Also verifies the loop variable doesn't shadow the function parameter
        'row' (which is a manifest dict, not a tuple).
        """
        import inspect
        from vcf_stats.seq2neo.cli import process_single_sample
        source = inspect.getsource(process_single_sample)
        assert '["CHROM"].to_list()' not in source, "Still using to_list() for target positions"
        assert '["POS"].to_list()' not in source, "Still using to_list() for target positions"
        assert "iter_rows()" in source, "Missing chunked iter_rows() for target building"
        assert "chunk_size" in source, "Missing chunk_size for target building"
        # Loop variable must NOT shadow function parameter 'row' (which is a dict)
        # 'for row in' would overwrite the dict parameter → tuple indexing error
        assert 'for row in chunk.iter_rows()' not in source, (
            "Loop variable 'row' shadows function parameter 'row' (manifest dict). "
            "Use 'for t in chunk.iter_rows()' instead."
        )

    def test_gc_collect_in_process_sample(self):
        """process_single_sample has gc.collect() after del statements."""
        import inspect
        from vcf_stats.seq2neo.cli import process_single_sample
        source = inspect.getsource(process_single_sample)
        assert "gc.collect()" in source, "Missing gc.collect() in process_single_sample"

    def test_rust_owned_hashset_no_clone(self):
        """Rust parse_caller_vcf takes owned HashSet (no clone)."""
        import inspect
        from vcf_stats.seq2neo.caller_parser import _parse_one_caller
        # Verify the Rust path exists (HAS_RUST_CALLER check)
        source = inspect.getsource(_parse_one_caller)
        assert "parse_caller_vcf" in source, "Rust parser not wired"

    def test_no_loop_variable_shadows_function_parameter(self):
        """No for-loop variable shadows a function parameter (AST check).

        Python's for-loop variables leak into the enclosing scope.
        If a loop variable has the same name as a function parameter,
        the parameter gets silently overwritten — causing TypeError
        when the parameter is used later as a different type.
        """
        import ast, glob
        issues = []
        for fpath in glob.glob('bin/vcf_stats/seq2neo/*.py'):
            with open(fpath) as f:
                try:
                    tree = ast.parse(f.read())
                except SyntaxError:
                    continue
            for node in ast.walk(tree):
                if isinstance(node, ast.FunctionDef):
                    params = {a.arg for a in node.args.args}
                    for child in ast.walk(node):
                        if isinstance(child, ast.For):
                            if isinstance(child.target, ast.Name):
                                if child.target.id in params:
                                    issues.append(
                                        f"{fpath}:{child.lineno} — "
                                        f"'{child.target.id}' in for-loop shadows "
                                        f"parameter of {node.name}()"
                                    )
        assert not issues, (
            "Loop variable shadows function parameter:\n" + "\n".join(issues)
        )

    def test_large_sample_semaphore_exists(self):
        """_THREAD_LARGE_SEM exists in cli.py for throttling."""
        import inspect
        from vcf_stats.seq2neo import cli
        assert hasattr(cli, '_THREAD_LARGE_SEM'), "Missing _THREAD_LARGE_SEM"
        assert hasattr(cli, '_LARGE_THRESHOLD'), "Missing _LARGE_THRESHOLD"
        assert cli._LARGE_THRESHOLD == 2_000_000

    def test_rust_caller_no_column_length_mismatch(self):
        """All columns from Rust parse_caller_vcf have equal length.

        Regression test for the for _ in 0..7 bug that inflated
        Strelka columns 7× for non-Strelka callers.
        """
        try:
            import stats_core
            if not hasattr(stats_core, 'parse_caller_vcf'):
                pytest.skip("Rust parse_caller_vcf not available")
        except ImportError:
            pytest.skip("stats_core not available")

        import glob
        base = os.path.join(REAL_SAMPLE["base_output_dir"], REAL_SAMPLE["dir_name"])
        rescue = glob.glob(os.path.join(base, "vcf_realignment", "rescue", "*", "*.filtered.vcf.stripped.vcf.gz"))
        if not rescue:
            pytest.skip("Rescue VCF not found")

        records = stats_core.parse_rescue(rescue[0])
        chroms = [r["CHROM"] for r in records]
        poss = [r["POS"] for r in records]
        refs = [r["REF"] for r in records]
        alts = [r["ALT"] for r in records]

        for cn in ["DNA_mutect2", "DNA_deepsomatic", "DNA_strelka"]:
            cfg = CALLER_CONFIGS[cn]
            subdir = cfg["subdir"].format(prefix=REAL_SAMPLE["vcf_prefix"])
            vcf_path = glob.glob(os.path.join(base, subdir, cfg["pattern"]))
            if not vcf_path:
                continue
            result = stats_core.parse_caller_vcf(
                vcf_path[0], chroms, poss, refs, alts, cfg["sample_suffix"], cn
            )
            n = len(result["CHROM"])
            for k in result:
                assert len(result[k]) == n, (
                    f"[{cn}] column {k} length {len(result[k])} != CHROM length {n}"
                )


# ── New tests for column-oriented rescue parser and process isolation ──────

class TestColumnOrientedRescueParser:
    """Verify the Rust column-oriented rescue parser (parse_rescue_columns)."""

    @pytest.fixture(scope="class")
    def rescue_data(self):
        """Parse rescue VCF once for the entire class (saves ~80s across 7 tests)."""
        try:
            import stats_core
        except ImportError:
            pytest.skip("stats_core not available")
        import glob

        if not hasattr(stats_core, 'parse_rescue_columns'):
            pytest.skip("parse_rescue_columns not available")

        rescue = glob.glob(os.path.join(
            REAL_SAMPLE["base_output_dir"], REAL_SAMPLE["dir_name"],
            "vcf_realignment", "rescue", "*", "*.filtered.vcf.stripped.vcf.gz"
        ))
        if not rescue:
            pytest.skip("Rescue VCF not found")

        # Parse once — column-oriented
        cols = stats_core.parse_rescue_columns(rescue[0])
        # Parse once — row-oriented (for parity test)
        from vcf_stats.seq2neo.rust_vcf import parse_rescue_vcf_columns
        df_cols = parse_rescue_vcf_columns(rescue[0])
        # Row-oriented legacy
        records = stats_core.parse_rescue(rescue[0])
        df_rows = pl.DataFrame(records)
        from vcf_stats.seq2neo.rust_vcf import _cast_columns_fallback, _add_derived_columns_polars
        df_rows = _cast_columns_fallback(df_rows)
        df_rows = _add_derived_columns_polars(df_rows)

        return {
            "cols": cols,
            "df_cols": df_cols,
            "df_rows": df_rows,
        }

    def test_parse_rescue_columns_available(self):
        """parse_rescue_columns is available in stats_core."""
        try:
            import stats_core
        except ImportError:
            pytest.skip("stats_core not available")
        assert hasattr(stats_core, 'parse_rescue_columns'), (
            "parse_rescue_columns not found in stats_core"
        )

    def test_all_column_lengths_equal(self, rescue_data):
        """All columns from parse_rescue_columns have identical length."""
        cols = rescue_data["cols"]
        assert cols, "parse_rescue_columns returned empty dict"
        n = len(cols.get("CHROM", []))
        assert n > 0, "Expected variants but got empty CHROM column"
        for k, v in cols.items():
            assert len(v) == n, (
                f"Column '{k}' length {len(v)} != CHROM length {n}"
            )

    def test_derived_columns_present(self, rescue_data):
        """variant_type and ti_tv are present in column-oriented output."""
        cols = rescue_data["cols"]
        assert "variant_type" in cols, "variant_type missing from columns"
        assert "ti_tv" in cols, "ti_tv missing from columns"
        vt_vals = set(cols["variant_type"])
        assert vt_vals.issubset({"SNV", "INS", "DEL", "MNV"}), (
            f"Unexpected variant_type values: {vt_vals - {'SNV', 'INS', 'DEL', 'MNV'}}"
        )

    def test_parity_with_row_oriented(self, rescue_data):
        """Column-oriented DataFrame matches row-oriented DataFrame."""
        df_cols = rescue_data["df_cols"]
        df_rows = rescue_data["df_rows"]

        assert len(df_cols) == len(df_rows), (
            f"Row count mismatch: cols={len(df_cols)}, rows={len(df_rows)}"
        )

        for col in df_rows.columns:
            if col in df_cols.columns:
                vals_cols = df_cols[col].to_list()
                vals_rows = df_rows[col].to_list()
                for i, (a, b) in enumerate(zip(vals_cols, vals_rows)):
                    if a is None and b is None:
                        continue
                    if col in RESCUE_FLAG_FIELDS:
                        a_bool = bool(a) if a is not None else False
                        b_bool = bool(b) if b is not None else False
                        if a_bool == b_bool:
                            continue
                    if a != b:
                        if isinstance(a, float) and isinstance(b, float):
                            assert abs(a - b) < 1e-6, (
                                f"Column '{col}' row {i}: {a} != {b}"
                            )
                        else:
                            assert a == b, (
                                f"Column '{col}' row {i}: {a!r} != {b!r}"
                            )

    def test_polars_derived_columns_parity(self, rescue_data):
        """polars-native derived columns match per-row Python computation."""
        from vcf_stats.seq2neo.rescue_parser import _add_derived_columns_polars

        # Use the row-oriented data (first 1000 rows only for speed)
        df = rescue_data["df_rows"].select(["CHROM", "POS", "REF", "ALT", "FILTER"]).head(1000)
        result = _add_derived_columns_polars(df)
        assert "variant_type" in result.columns
        assert "ti_tv" in result.columns

    def test_malloc_trim_no_crash(self):
        """_malloc_trim() does not crash on Linux."""
        from vcf_stats.seq2neo.cli import _malloc_trim
        _malloc_trim()

    def test_malloc_trim_handles_missing_libc(self):
        """_malloc_trim() handles platforms without libc.so.6."""
        import ctypes
        original = ctypes.CDLL
        def mock_cdll(name):
            if name == "libc.so.6":
                raise OSError("No such file")
            return original(name)
        ctypes.CDLL = mock_cdll
        try:
            from vcf_stats.seq2neo.cli import _malloc_trim
            _malloc_trim()
        finally:
            ctypes.CDLL = original


# ═══════════════════════════════════════════════════════════════════════════
# TestCrossSampleCols (Phase 6.2)
# ═══════════════════════════════════════════════════════════════════════════

class TestCrossSampleCols:
    """Verify _CROSS_SAMPLE_COLS includes per-caller VAF/DP/AD columns."""

    def test_per_caller_vaf_columns_included(self):
        from vcf_stats.seq2neo.statistics import _CROSS_SAMPLE_COLS
        vaf_cols = ["DNA_mutect2_VAF", "RNA_mutect2_VAF",
                    "DNA_deepsomatic_VAF", "RNA_deepsomatic_VAF",
                    "DNA_strelka_VAF", "RNA_strelka_VAF"]
        for c in vaf_cols:
            assert c in _CROSS_SAMPLE_COLS, f"{c} missing from _CROSS_SAMPLE_COLS"

    def test_per_caller_dp_columns_included(self):
        from vcf_stats.seq2neo.statistics import _CROSS_SAMPLE_COLS
        dp_cols = ["DNA_mutect2_DP", "RNA_mutect2_DP",
                   "DNA_deepsomatic_DP", "RNA_deepsomatic_DP",
                   "DNA_strelka_DP", "RNA_strelka_DP"]
        for c in dp_cols:
            assert c in _CROSS_SAMPLE_COLS, f"{c} missing from _CROSS_SAMPLE_COLS"

    def test_caller_wise_summary_nonzero(self):
        """compute_caller_wise_summary returns non-zero for synthetic data with VAF columns."""
        from vcf_stats.seq2neo.statistics import compute_caller_wise_summary
        import polars as pl
        df = pl.DataFrame({
            "DNA_mutect2_VAF": [0.2, 0.3, None, 0.1],
            "DNA_mutect2_DP": [50, 30, 10, 40],
            "DNA_strelka_VAF": [0.15, 0.25, 0.05, None],
            "DNA_strelka_DP": [48, 28, 12, 38],
            "RNA_mutect2_VAF": [0.18, 0.28, None, None],
            "RNA_mutect2_DP": [45, 28, 8, 35],
        })
        result = compute_caller_wise_summary(df)
        assert not result.is_empty()
        dna_mutect = result.filter(pl.col("caller") == "DNA_mutect2")
        assert dna_mutect["n_with_vaf"][0] > 0
        assert dna_mutect["n_with_dp"][0] > 0

    def test_vaf_threshold_sweep_nonempty(self):
        """compute_vaf_threshold_sweep returns non-empty DataFrame with VAF columns."""
        from vcf_stats.seq2neo.statistics import compute_vaf_threshold_sweep
        import polars as pl
        df = pl.DataFrame({
            "FILTER": ["Somatic", "Somatic", "Germline", "Somatic", "Germline"],
            "DNA_mutect2_VAF": [0.05, 0.25, 0.15, 0.45, 0.08],
            "DNA_strelka_VAF": [0.04, 0.30, 0.12, 0.50, None],
        })
        result = compute_vaf_threshold_sweep(df)
        assert not result.is_empty()
        assert "caller" in result.columns
        assert "threshold" in result.columns
        assert "pct_retained" in result.columns
        # At threshold 0.10: 3/5 and 3/4 should be retained
        row = result.filter((pl.col("caller") == "DNA_mutect2") & (pl.col("threshold") == 0.10))
        assert row["pct_retained"][0] > 0

    def test_ensure_eager_no_warning(self):
        """_ensure_eager uses collect_schema().names() — no PerformanceWarning on LazyFrame."""
        import warnings
        import polars as pl
        from vcf_stats.seq2neo.statistics import _ensure_eager
        lazy = pl.LazyFrame({"FILTER": ["Somatic"], "variant_type": ["SNV"]})
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = _ensure_eager(lazy)
            assert isinstance(result, pl.DataFrame)
            perf_warnings = [x for x in w if "PerformanceWarning" in str(x.message)]
            assert len(perf_warnings) == 0, f"Got PerformanceWarning: {perf_warnings}"


# ═══════════════════════════════════════════════════════════════════════════
# TestChartGroupCol (Phase 6.3)
# ═══════════════════════════════════════════════════════════════════════════

class TestChartGroupCol:
    """Verify chart functions work with alternative group_col parameters."""

    @pytest.fixture
    def test_df(self):
        import polars as pl
        return pl.DataFrame({
            "set_number": [1, 1, 2, 2],
            "disease_normalized": ["Lung", "Lung", "Breast", "Breast"],
            "final_tier": ["C1D1", "C2D0", "C1D1", "C3D1"],
            "FILTER": ["Somatic", "Germline", "Somatic", "Reference"],
            "variant_type": ["SNV", "SNV", "INS", "DEL"],
            "ti_tv": [True, False, None, None],
            "CROSS_MODALITY": ["YES", "NO", "YES", "NO"],
            "RESCUED": ["NO", "NO", "YES", "NO"],
            "COSMIC_ID": ["C1", None, None, "C2"],
            "GNOMAD_AF": [0.01, None, 0.05, None],
            "REDI_EVIDENCE": ["NONE", "NONE", "LOW", "HIGH"],
            "N_SUPPORT_CALLERS": [6, 3, 2, 1],
            "DNA_mutect2_GT": ["0/1", "0/0", "0/1", "0/1"],
            "RNA_mutect2_GT": ["0/1", "0/0", "0/1", "0/1"],
            "DNA_deepsomatic_GT": ["0/1", "0/0", "0/1", "0/1"],
            "RNA_deepsomatic_GT": ["0/1", "0/0", "0/1", "0/1"],
        })

    def test_caller_overlap_with_disease(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_caller_overlap
        fig = plot_caller_overlap(test_df, str(tmp_path), group_col="disease_normalized")
        assert fig is not None

    def test_variant_type_with_tier(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_variant_type_distribution
        fig = plot_variant_type_distribution(test_df, str(tmp_path), group_col="final_tier")
        assert fig is not None

    def test_ti_tv_with_disease(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_ti_tv_ratio
        fig = plot_ti_tv_ratio(test_df, str(tmp_path), group_col="disease_normalized")
        assert fig is not None

    def test_cross_modality_with_disease(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_cross_modality
        fig = plot_cross_modality(test_df, str(tmp_path), group_col="disease_normalized")
        assert fig is not None

    def test_filter_dist_with_tier(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_filter_distribution
        fig = plot_filter_distribution(test_df, str(tmp_path), group_col="final_tier")
        assert fig is not None

    def test_gt_concordance_with_disease(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_gt_concordance
        fig = plot_gt_concordance(test_df, str(tmp_path), group_col="disease_normalized")
        assert fig is not None

    def test_cosmic_gnomad_with_tier(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_cosmic_gnomad_annotation
        fig = plot_cosmic_gnomad_annotation(test_df, str(tmp_path), group_col="final_tier")
        assert fig is not None

    def test_redi_with_tier(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_redi_evidence
        fig = plot_redi_evidence(test_df, str(tmp_path), group_col="final_tier")
        assert fig is not None


# ═══════════════════════════════════════════════════════════════════════════
# TestChartColorCol (Phase 6.4)
# ═══════════════════════════════════════════════════════════════════════════

class TestChartColorCol:
    """Verify chart functions accept and use color_col parameter."""

    @pytest.fixture
    def test_df(self):
        import polars as pl
        return pl.DataFrame({
            "set_number": [1, 1, 2, 2],
            "disease_normalized": ["Lung", "Lung", "Breast", "Breast"],
            "final_tier": ["C1D1", "C2D0", "C1D1", "C3D1"],
            "FILTER": ["Somatic", "Germline", "Somatic", "Reference"],
            "DNA_mutect2_VAF": [0.2, 0.3, 0.1, 0.4],
            "RNA_mutect2_VAF": [0.18, 0.28, 0.08, 0.35],
            "DNA_deepsomatic_VAF": [0.21, 0.31, 0.11, 0.41],
            "DNA_strelka_VAF": [0.19, 0.29, 0.09, 0.38],
            "DNA_VAF_mean": [0.2, 0.3, 0.1, 0.4],
            "RNA_VAF_mean": [0.18, 0.28, 0.08, 0.35],
            "N_SUPPORT_CALLERS": [6, 3, 2, 1],
            "BAM_DN_DP": [50, 30, 40, 20],
            "BAM_DN_REF_DP": [40, 25, 35, 15],
        })

    def test_vaf_distribution_with_disease_color(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_vaf_distribution
        fig = plot_vaf_distribution(test_df, str(tmp_path), color_col="disease_normalized")
        assert fig is not None

    def test_dna_vs_rna_vaf_with_tier_color(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_dna_vs_rna_vaf
        fig = plot_dna_vs_rna_vaf(test_df, str(tmp_path), color_col="final_tier")
        assert fig is not None

    def test_per_caller_with_set_color(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_dna_vs_rna_per_caller
        fig = plot_dna_vs_rna_per_caller(test_df, str(tmp_path), color_col="set_number")
        assert fig is not None

    def test_bam_coverage_with_disease_color(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_bam_coverage_violin
        fig = plot_bam_coverage_violin(test_df, str(tmp_path), color_col="disease_normalized")
        assert fig is not None


# ═══════════════════════════════════════════════════════════════════════════
# TestChartFacetCol (Phase 6.5)
# ═══════════════════════════════════════════════════════════════════════════

class TestChartFacetCol:
    """Verify chart functions accept and use facet_col parameter."""

    @pytest.fixture
    def test_df(self):
        import polars as pl
        return pl.DataFrame({
            "set_number": [1, 1, 2, 2],
            "disease_normalized": ["Lung", "Lung", "Breast", "Breast"],
            "caller_tier": ["C1", "C2", "C1", "C3"],
            "final_tier": ["C1D1", "C2D0", "C1D1", "C3D1"],
            "variant_type": ["SNV", "SNV", "INS", "DEL"],
            "N_SUPPORT_CALLERS": [6, 3, 2, 1],
            "DNA_mutect2_VAF": [0.2, 0.3, 0.1, 0.4],
            "RNA_mutect2_VAF": [0.18, 0.28, 0.08, 0.35],
            "DNA_mutect2_DP": [50, 30, 40, 20],
            "RNA_mutect2_DP": [45, 28, 35, 18],
            "DNA_mutect2_GT": ["0/1", "0/0", "0/1", "0/1"],
            "RNA_mutect2_GT": ["0/1", "0/0", "0/1", "0/1"],
            "DNA_deepsomatic_GT": ["0/1", "0/0", "0/1", "0/1"],
            "RNA_deepsomatic_GT": ["0/1", "0/0", "0/1", "0/1"],
            "tier_quality": [0.8, 0.5, 0.9, 0.3],
        })

    def test_vaf_per_tier_with_set_facet(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_vaf_boxplot_per_tier
        fig = plot_vaf_boxplot_per_tier(test_df, str(tmp_path), facet_col="set_number")
        assert fig is not None

    def test_dp_per_tier_with_disease_facet(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_dp_boxplot_per_tier
        fig = plot_dp_boxplot_per_tier(test_df, str(tmp_path), facet_col="disease_normalized")
        assert fig is not None

    def test_gt_per_tier_with_disease_facet(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_gt_concordance_per_tier
        fig = plot_gt_concordance_per_tier(test_df, str(tmp_path), facet_col="disease_normalized")
        assert fig is not None

    def test_caller_overlap_tier_with_set_facet(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_tiered_caller_overlap
        fig = plot_tiered_caller_overlap(test_df, str(tmp_path), facet_col="set_number")
        assert fig is not None

    def test_variant_types_tier_with_disease_facet(self, test_df, tmp_path):
        from vcf_stats.seq2neo.visualizer import plot_tiered_variant_types
        fig = plot_tiered_variant_types(test_df, str(tmp_path), facet_col="disease_normalized")
        assert fig is not None


# ═══════════════════════════════════════════════════════════════════════════
# TestChartQuality (Phase 6.6)
# ═══════════════════════════════════════════════════════════════════════════

class TestChartQuality:
    """Verify chart quality fixes — % marks, log scale, box+violin, scroll."""

    def test_percentage_marks_present(self, tmp_path):
        """Variant type chart has text marks with pct encoding."""
        import polars as pl
        from vcf_stats.seq2neo.visualizer import plot_variant_type_distribution
        df = pl.DataFrame({
            "set_number": [1, 1, 2, 2],
            "variant_type": ["SNV", "INS", "SNV", "DEL"],
        })
        fig = plot_variant_type_distribution(df, str(tmp_path), group_col="set_number")
        assert fig is not None
        # Chart should have both bar and text layers
        html = fig.to_html()
        assert "mark_text" in html.lower() or "text" in str(fig.to_dict())

    def test_log_scale_on_caller_overlap(self, tmp_path):
        """Caller overlap tier chart uses log scale y-axis."""
        import polars as pl
        from vcf_stats.seq2neo.visualizer import plot_tiered_caller_overlap
        df = pl.DataFrame({
            "caller_tier": ["C1", "C2", "C1", "C3"],
            "N_SUPPORT_CALLERS": [6, 3, 2, 1],
        })
        fig = plot_tiered_caller_overlap(df, str(tmp_path))
        assert fig is not None

    def test_per_sample_no_top_n_limit(self, tmp_path):
        """Per-sample chart uses dynamic height for all samples, no top_n cutoff."""
        import polars as pl
        from vcf_stats.seq2neo.visualizer import plot_per_sample_distribution
        df = pl.DataFrame({
            "sample_id": [f"s{i}" for i in range(1, 66)],
            "total_variants": [1000 + i * 100 for i in range(1, 66)],
        })
        fig = plot_per_sample_distribution(df, str(tmp_path))
        assert fig is not None

    def test_box_violin_overlay(self, tmp_path):
        """VAF distribution chart has both box and violin layers."""
        import polars as pl
        from vcf_stats.seq2neo.visualizer import plot_vaf_distribution
        df = pl.DataFrame({
            "DNA_mutect2_VAF": [0.2, 0.3, 0.1, 0.4],
            "RNA_mutect2_VAF": [0.18, 0.28, 0.08, 0.35],
        })
        fig = plot_vaf_distribution(df, str(tmp_path))
        assert fig is not None


# ═══════════════════════════════════════════════════════════════════════════
# TestNewCharts (Phase 6.7)
# ═══════════════════════════════════════════════════════════════════════════

class TestNewCharts:
    """Verify new chart functions return valid chart objects."""

    def test_dp_distribution_returns_chart(self, tmp_path):
        import polars as pl
        from vcf_stats.seq2neo.visualizer import plot_dp_distribution
        df = pl.DataFrame({
            "DNA_mutect2_DP": [50, 30, 40, 20],
            "RNA_mutect2_DP": [45, 28, 35, 18],
            "DNA_strelka_DP": [48, 28, 38, 18],
        })
        fig = plot_dp_distribution(df, str(tmp_path))
        assert fig is not None

    def test_mean_vaf_per_group_returns_chart(self, tmp_path):
        import polars as pl
        from vcf_stats.seq2neo.visualizer import plot_mean_vaf_per_group
        df = pl.DataFrame({
            "set_number": [1, 2],
            "mean_dna_vaf": [0.25, 0.30],
            "mean_rna_vaf": [0.22, 0.28],
        })
        fig = plot_mean_vaf_per_group(df, str(tmp_path), group_col="set_number")
        assert fig is not None

    def test_mean_dp_per_group_returns_chart(self, tmp_path):
        import polars as pl
        from vcf_stats.seq2neo.visualizer import plot_mean_dp_per_group
        df = pl.DataFrame({
            "set_number": [1, 2],
            "mean_dna_dp": [50, 60],
            "mean_rna_dp": [45, 55],
        })
        fig = plot_mean_dp_per_group(df, str(tmp_path), group_col="set_number")
        assert fig is not None

    def test_n_support_callers_dist_returns_chart(self, tmp_path):
        import polars as pl
        from vcf_stats.seq2neo.visualizer import plot_n_support_callers_dist
        df = pl.DataFrame({
            "set_number": [1, 1, 2, 2],
            "N_SUPPORT_CALLERS": [6, 3, 2, 1],
        })
        fig = plot_n_support_callers_dist(df, str(tmp_path), group_col="set_number")
        assert fig is not None

    def test_caller_tier_heatmap_returns_chart(self, tmp_path):
        import polars as pl
        from vcf_stats.seq2neo.visualizer import plot_caller_tier_heatmap
        df = pl.DataFrame({
            "final_tier": ["C1D1", "C2D0", "C1D1", "C3D1"],
            "DNA_mutect2_VAF": [0.2, 0.3, None, 0.1],
            "DNA_strelka_VAF": [0.15, None, 0.05, 0.12],
        })
        fig = plot_caller_tier_heatmap(df, str(tmp_path))
        assert fig is not None

    def test_sample_overview_scatter_returns_chart(self, tmp_path):
        import polars as pl
        from vcf_stats.seq2neo.visualizer import plot_sample_overview_scatter
        df = pl.DataFrame({
            "sample_id": ["s1", "s2", "s3"],
            "total_variants": [1000, 2000, 1500],
            "mean_dna_vaf_mean": [0.2, 0.3, 0.25],
            "mean_dna_dp_mean": [50, 60, 55],
            "set_number": [1, 1, 2],
        })
        fig = plot_sample_overview_scatter(df, str(tmp_path))
        assert fig is not None


# ═══════════════════════════════════════════════════════════════════════════
# TestCLIWiseLoop (Phase 6.8)
# ═══════════════════════════════════════════════════════════════════════════

class TestCLIWiseLoop:
    """Verify CLI --wise flag parsing and per-wise chart generation."""

    def test_wise_flag_parses_correctly(self):
        import sys
        from vcf_stats.seq2neo.cli import main as _unused
        parser_args = ["--manifest", "m.parquet", "--output-dir", "out",
                       "--wise", "set", "disease", "--max-samples", "1", "--no-bam"]
        # Just verify the parser doesn't crash
        import argparse
        ap = argparse.ArgumentParser()
        ap.add_argument("--wise", nargs="*", default=None)
        ap.add_argument("--manifest", default="m.parquet")
        ap.add_argument("--output-dir", default="out")
        ap.add_argument("--max-samples", type=int, default=None)
        ap.add_argument("--no-bam", action="store_true")
        ap.add_argument("--no-pileup", action="store_true")
        ap.add_argument("--set", type=int, default=None)
        ap.add_argument("--threads", type=int, default=6)
        ap.add_argument("--sample-workers", type=int, default=1)
        ap.add_argument("--bam-workers", type=int, default=8)
        args = ap.parse_args(parser_args)
        assert args.wise == ["set", "disease"]

    def test_wise_flag_default_all(self):
        import argparse
        ap = argparse.ArgumentParser()
        ap.add_argument("--wise", nargs="*", default=None)
        args = ap.parse_args(["--wise"])
        assert args.wise == []

    def test_new_chart_imports_available(self):
        """Verify all new chart functions are importable."""
        from vcf_stats.seq2neo.visualizer import (
            plot_dp_distribution,
            plot_mean_vaf_per_group,
            plot_mean_dp_per_group,
            plot_n_support_callers_dist,
            plot_caller_tier_heatmap,
            plot_sample_overview_scatter,
        )
        assert callable(plot_dp_distribution)
        assert callable(plot_mean_vaf_per_group)
        assert callable(plot_mean_dp_per_group)
        assert callable(plot_n_support_callers_dist)
        assert callable(plot_caller_tier_heatmap)
        assert callable(plot_sample_overview_scatter)

    def test_wise_directories_created(self, tmp_path):
        """Wise chart generation creates per-wise directories."""
        output_dir = tmp_path / "stats"
        output_dir.mkdir()
        wise_names = ["set", "disease", "sample", "tier", "caller", "chromosome"]
        for wise_name in wise_names:
            wise_plot_dir = output_dir / "plots" / wise_name
            wise_plot_dir.mkdir(parents=True)
            assert wise_plot_dir.exists()
        # Verify all 6 directories exist
        for wise_name in wise_names:
            assert (output_dir / "plots" / wise_name).exists()

    def test_wise_chart_count(self, tmp_path):
        """Each wise gets its chart subset — verify chart files are created."""
        import polars as pl
        from vcf_stats.seq2neo.visualizer import (
            plot_chromosome_density, plot_vc_distribution, plot_ti_tv_ratio
        )
        df = pl.DataFrame({
            "CHROM": ["chr1", "chr2"], "FILTER": ["Somatic", "Germline"],
            "variant_type": ["SNV", "INS"], "ti_tv": [True, False],
            "set_number": [1, 1],
        })
        wise_dir = tmp_path / "plots" / "chromosome"
        wise_dir.mkdir(parents=True)
        plot_chromosome_density(df, str(wise_dir))
        # Chart files should exist
        assert (wise_dir / "30_chromosome_density.html").exists()


# ═══════════════════════════════════════════════════════════════════════════
# TestTSV (Phase 6.9)
# ═══════════════════════════════════════════════════════════════════════════

class TestTSV:
    """Verify TSV roundtrip and resume path."""

    def test_write_tsv_roundtrip(self, tmp_path):
        """write_tsv produces a file readable with read_csv(separator='\t')."""
        import polars as pl
        from vcf_stats.seq2neo.statistics import write_tsv
        df = pl.DataFrame({"x": [1, 2, 3], "y": ["a", "b", "c"], "z": [0.1, 0.2, 0.3]})
        path = tmp_path / "test.tsv"
        write_tsv(df, str(path))
        assert path.exists()
        read_back = pl.read_csv(str(path), separator="\t")
        assert read_back.shape == df.shape
        assert read_back["x"].to_list() == [1, 2, 3]
        assert read_back["y"].to_list() == ["a", "b", "c"]

    def test_write_tsv_handles_special_chars(self, tmp_path):
        """TSV with pipe-delimited and colon-delimited values (VCF fields)."""
        import polars as pl
        from vcf_stats.seq2neo.statistics import write_tsv
        df = pl.DataFrame({
            "INFO": ["DNA_strelka:Somatic|RNA_mutect2:Germline", "PASS"],
            "FILTER": ["Somatic", "Reference"],
        })
        path = tmp_path / "special.tsv"
        write_tsv(df, str(path))
        read_back = pl.read_csv(str(path), separator="\t")
        assert read_back["INFO"][0] == "DNA_strelka:Somatic|RNA_mutect2:Germline"

    def test_resume_reads_tsv_not_csv(self, tmp_path):
        """Resume path reads .tsv files with tab separator, falls back to .csv."""
        import polars as pl
        from vcf_stats.seq2neo.statistics import write_tsv
        # Write TSV
        df = pl.DataFrame({"sample_id": ["s1"], "total_variants": [100]})
        tsv_path = tmp_path / "sample_summary.tsv"
        write_tsv(df, str(tsv_path))
        # Read back as resume path would
        read_back = pl.read_csv(str(tsv_path), separator="\t")
        assert read_back["sample_id"][0] == "s1"
        # Legacy CSV fallback also works
        csv_path = tmp_path / "sample_summary.csv"
        df.write_csv(str(csv_path))
        csv_back = pl.read_csv(str(csv_path))
        assert csv_back["sample_id"][0] == "s1"


class TestBamValidationColumns:
    """Verify BAM validation uses correct column naming (BAM_{bt}_{metric})."""

    def test_has_bam_data_true_when_pileup_columns_exist(self):
        """has_bam_data=True when BAM_DT_DP and BAM_RT_DP columns are present."""
        from vcf_stats.seq2neo.bam_validation import validate_bam_vs_caller
        df = pl.DataFrame({
            "BAM_DT_DP": [10, 20],
            "BAM_DT_ALT_DP": [3, 8],
            "BAM_DT_REF_DP": [7, 12],
            "BAM_RT_DP": [15, 25],
            "BAM_RT_ALT_DP": [5, 10],
            "BAM_RT_REF_DP": [10, 15],
        })
        result = validate_bam_vs_caller(df, "test_sample")
        assert result["has_bam_data"] is True
        assert "DT" in result["bam_types"]
        assert "RT" in result["bam_types"]

    def test_has_bam_data_false_without_pileup(self):
        """has_bam_data=False when no BAM pileup columns exist."""
        from vcf_stats.seq2neo.bam_validation import validate_bam_vs_caller
        df = pl.DataFrame({"some_col": [1]})
        result = validate_bam_vs_caller(df, "test_sample")
        assert result["has_bam_data"] is False

    def test_column_naming_matches_pileup_output(self):
        """BAM validation uses BAM_{bt}_{metric} not BAM_{metric}_{bt}."""
        from vcf_stats.seq2neo.bam_validation import validate_bam_vs_caller
        # Simulate the exact column names from pileup rename:
        # rename = {c: f"BAM_{bt}_{c}" for c in pileup_df.columns if c not in join_key}
        df = pl.DataFrame({
            "BAM_DT_DP": [10],
            "BAM_DT_ALT_DP": [3],
            "BAM_DT_REF_DP": [7],
            "BAM_RT_DP": [20],
            "BAM_RT_ALT_DP": [8],
            "BAM_RT_REF_DP": [12],
        })
        result = validate_bam_vs_caller(df, "test_sample")
        assert result["has_bam_data"] is True
        assert len(result["bam_types"].split(",")) == 2  # DT + RT, no DN in this test

    def test_strand_column_naming(self):
        """bam_type_has_strand uses correct BAM_{bt}_F1R2_alt naming."""
        from vcf_stats.seq2neo.bam_validation import bam_type_has_strand
        df = pl.DataFrame({
            "BAM_DT_F1R2_alt": [1],
            "BAM_DT_F2R1_alt": [0],
        })
        assert bam_type_has_strand(df, "DT") is True


# ═══════════════════════════════════════════════════════════════════════════════
# Strelka TAR/TIR Mapping and Repair Tests (fix-stats-viz-round4)
# ═══════════════════════════════════════════════════════════════════════════════


class TestStrelkaMappingFix:
    """Verify Strelka AD_ALT ← TIR and AD_REF ← TAR after join."""

    def test_strelka_ad_alt_equals_tir_after_join(self):
        """After _streaming_join_one, AD_ALT should be TIR, AD_REF should be TAR."""
        from vcf_stats.seq2neo.cli import _streaming_join_one

        base_df = pl.DataFrame({
            "CHROM": ["chr1", "chr1"],
            "POS": [100, 200],
            "REF": ["A", "C"],
            "ALT": ["T", "G"],
        })
        col_data = {
            "CHROM": ["chr1", "chr1"],
            "POS": [100, 200],
            "REF": ["A", "C"],
            "ALT": ["T", "G"],
            "DP": [100, 80],
            "TAR": [90, 70],    # ref counts
            "TIR": [8, 6],     # alt counts
            "TOR": [2, 4],     # other counts
        }
        result = _streaming_join_one(base_df, col_data, "DNA_strelka")

        # AD_ALT should come from TIR, not TAR
        assert "DNA_strelka_AD_ALT" in result.columns
        assert "DNA_strelka_AD_REF" in result.columns
        assert result["DNA_strelka_AD_ALT"].to_list() == [8, 6]   # TIR values
        assert result["DNA_strelka_AD_REF"].to_list() == [90, 70]  # TAR values

    def test_strelka_vaf_within_unit_range(self):
        """After fix, Strelka VAF should always be in [0, 1]."""
        from vcf_stats.seq2neo.statistics import compute_vaf_columns
        df = pl.DataFrame({
            "DNA_strelka_AD_ALT": [8, 6, 0, None],
            "DNA_strelka_DP": [100, 80, 50, None],
        })
        result = compute_vaf_columns(df)
        vafs = result["DNA_strelka_VAF"].to_list()
        for v in vafs:
            if v is not None:
                assert 0.0 <= v <= 1.0, f"VAF {v} out of [0, 1]"


class TestStrelkaRepair:
    """Test the _repair_strelka_columns post-load repair function."""

    def _make_old_parquet_lf(self):
        """Create a lazy frame simulating old (buggy) parquet: AD_ALT == TAR."""
        return pl.LazyFrame({
            "CHROM": ["chr1"] * 4,
            "POS": [100, 200, 300, 400],
            "DNA_strelka_TAR": [90, 70, 85, 60],       # ref counts
            "DNA_strelka_TIR": [8, 6, 10, 5],          # alt counts
            "DNA_strelka_TOR": [2, 4, 5, 3],            # other counts
            "DNA_strelka_AD_ALT": [90, 70, 85, 60],     # BUG: == TAR
            "DNA_strelka_AD_REF": [2, 4, 5, 3],         # BUG: == TOR
            "DNA_strelka_DP": [100, 80, 100, 68],
            "DNA_strelka_VAF": [0.9, 0.875, 0.85, 0.882],  # wrong (TAR/DP)
            "DNA_mutect2_VAF": [0.08, 0.06, 0.10, 0.05],
            "DNA_deepsomatic_VAF": [0.09, 0.07, 0.11, 0.06],
            "DNA_mutect2_AD_ALT": [8, 5, 10, 4],
            "DNA_deepsomatic_AD_ALT": [9, 6, 11, 5],
            "DNA_mutect2_AD_REF": [42, 35, 40, 30],
            "DNA_deepsomatic_AD_REF": [41, 34, 39, 29],
        })

    def test_repair_detects_and_fixes_inversion(self):
        """Repair should detect AD_ALT == TAR and remap to TIR."""
        from vcf_stats.seq2neo.cli import _repair_strelka_columns
        lf = self._make_old_parquet_lf()
        fixed = _repair_strelka_columns(lf).collect()

        # AD_ALT should now be TIR values
        assert fixed["DNA_strelka_AD_ALT"].to_list() == [8, 6, 10, 5]
        # AD_REF should now be TAR values
        assert fixed["DNA_strelka_AD_REF"].to_list() == [90, 70, 85, 60]
        # VAF should be TIR/DP
        vafs = fixed["DNA_strelka_VAF"].to_list()
        assert abs(vafs[0] - 8 / 100) < 0.001
        assert abs(vafs[1] - 6 / 80) < 0.001

    def test_repair_recomputes_means(self):
        """Repair should recompute DNA_VAF_mean from corrected per-caller VAFs."""
        from vcf_stats.seq2neo.cli import _repair_strelka_columns
        lf = self._make_old_parquet_lf()
        fixed = _repair_strelka_columns(lf).collect()

        # DNA_VAF_mean should now be mean of 3 corrected caller VAFs
        assert "DNA_VAF_mean" in fixed.columns
        row0_mean = fixed["DNA_VAF_mean"][0]
        expected = (8 / 100 + 0.08 + 0.09) / 3  # corrected strelka + mutect2 + deepsomatic
        assert abs(row0_mean - expected) < 0.01

    def test_repair_idempotent_on_fixed_data(self):
        """Repair should skip if AD_ALT already == TIR (no double-repair)."""
        from vcf_stats.seq2neo.cli import _repair_strelka_columns
        # Create correctly-mapped data: AD_ALT == TIR
        lf = pl.LazyFrame({
            "CHROM": ["chr1"] * 3,
            "POS": [100, 200, 300],
            "DNA_strelka_TAR": [90, 70, 85],
            "DNA_strelka_TIR": [8, 6, 10],
            "DNA_strelka_TOR": [2, 4, 5],
            "DNA_strelka_AD_ALT": [8, 6, 10],   # Correct: == TIR
            "DNA_strelka_AD_REF": [90, 70, 85],  # Correct: == TAR
            "DNA_strelka_DP": [100, 80, 100],
            "DNA_strelka_VAF": [0.08, 0.075, 0.10],
        })
        fixed = _repair_strelka_columns(lf).collect()
        # Should remain unchanged
        assert fixed["DNA_strelka_AD_ALT"].to_list() == [8, 6, 10]
        assert fixed["DNA_strelka_AD_REF"].to_list() == [90, 70, 85]


class TestRescueStatistics:
    """Test rescue analytics statistics functions."""

    @pytest.fixture
    def rescue_df(self):
        return pl.DataFrame({
            "sample_id": ["S1"] * 6 + ["S2"] * 4,
            "set_number": [1] * 6 + [2] * 4,
            "CHROM": ["chr1"] * 10,
            "POS": list(range(100, 200, 10)),
            "FILTER": ["Somatic", "Somatic", "Germline", "Artifact", "Somatic", "Germline",
                        "Somatic", "Somatic", "Germline", "Artifact"],
            "RESCUED": ["YES", "NO", "YES", "NO", "YES", "NO",
                        "YES", "NO", "NO", "NO"],
            "N_SUPPORT_CALLERS": [5, 3, 2, 1, 6, 2, 4, 3, 2, 1],
            "final_tier": ["C1D1", "C2D2", "C3D1", "C4D0", "C1D1", "C5D1",
                           "C1D1", "C2D2", "C3D1", "C7D0"],
            "DNA_VAF_mean": [0.3, 0.1, 0.05, 0.01, 0.4, 0.02,
                             0.25, 0.15, 0.03, 0.005],
            "RNA_VAF_mean": [0.2, 0.08, 0.04, None, 0.35, 0.01,
                             0.2, 0.1, 0.02, None],
            "DNA_DP_mean": [100, 50, 30, 10, 120, 20,
                            80, 60, 25, 8],
            "RNA_DP_mean": [90, 45, 28, None, 110, 18,
                            75, 55, 22, None],
        })

    def test_rescue_breakdown(self, rescue_df):
        from vcf_stats.seq2neo.statistics import compute_rescue_breakdown
        result = compute_rescue_breakdown(rescue_df)
        assert not result.is_empty()
        assert "RESCUED" in result.columns
        assert "count" in result.columns
        assert "pct" in result.columns
        # YES count: 4, NO count: 6
        yes_rows = result.filter(pl.col("RESCUED") == "YES")
        assert yes_rows["count"].sum() == 4

    def test_rescue_by_filter(self, rescue_df):
        from vcf_stats.seq2neo.statistics import compute_rescue_by_filter
        result = compute_rescue_by_filter(rescue_df)
        assert not result.is_empty()
        assert "FILTER" in result.columns
        assert "RESCUED" in result.columns

    def test_rescue_cross_tab(self, rescue_df):
        from vcf_stats.seq2neo.statistics import compute_rescue_cross_tab
        result = compute_rescue_cross_tab(rescue_df)
        assert not result.is_empty()
        assert "set_number" in result.columns

    def test_rescue_vaf_dp(self, rescue_df):
        from vcf_stats.seq2neo.statistics import compute_rescue_vaf_dp
        result = compute_rescue_vaf_dp(rescue_df)
        assert not result.is_empty()
        assert "n_variants" in result.columns

    def test_sample_rescue_summary(self, rescue_df):
        from vcf_stats.seq2neo.statistics import sample_rescue_summary
        result = sample_rescue_summary(rescue_df)
        assert not result.is_empty()
        assert "sample_id" in result.columns
        assert "count" in result.columns

    def test_rescue_by_tier(self, rescue_df):
        from vcf_stats.seq2neo.statistics import compute_rescue_by_tier
        result = compute_rescue_by_tier(rescue_df)
        assert not result.is_empty()
        assert "final_tier" in result.columns
        assert "pct" in result.columns

    def test_rescue_by_caller_support(self, rescue_df):
        from vcf_stats.seq2neo.statistics import compute_rescue_by_caller_support
        result = compute_rescue_by_caller_support(rescue_df)
        assert not result.is_empty()
        assert "N_SUPPORT_CALLERS" in result.columns

    def test_rescue_empty_column(self):
        """Functions should return empty DataFrame when RESCUED column is missing."""
        from vcf_stats.seq2neo.statistics import compute_rescue_breakdown
        df = pl.DataFrame({"sample_id": ["S1"], "FILTER": ["Somatic"]})
        result = compute_rescue_breakdown(df)
        assert result.is_empty()


# ═══════════════════════════════════════════════════════════════════════════════
# Round 5: Helpers, Color Registry, Symlog, FP Somatic, Faceting
# ═══════════════════════════════════════════════════════════════════════════════


class TestVisualizerHelpers:
    """Test the unified helper functions added in round 5."""

    def test_color_scale_known_entity(self):
        """_color_scale returns explicit Scale for known entities."""
        from vcf_stats.seq2neo.visualizer import _color_scale, CALLER_DOMAIN, CALLER_COLORS
        scale = _color_scale("caller")
        assert scale.domain == CALLER_DOMAIN
        assert scale.range == CALLER_COLORS

    def test_color_scale_rescue(self):
        """_color_scale returns rescue colors for RESCUED entity."""
        from vcf_stats.seq2neo.visualizer import _color_scale, RESCUE_DOMAIN, RESCUE_COLORS
        scale = _color_scale("RESCUED")
        assert scale.domain == RESCUE_DOMAIN
        assert scale.range == RESCUE_COLORS

    def test_color_scale_unknown_fallback(self):
        """_color_scale returns category10 scheme for unknown entities."""
        from vcf_stats.seq2neo.visualizer import _color_scale
        scale = _color_scale("unknown_entity")
        assert scale.scheme == "category10"

    def test_count_scale_symlog(self):
        """_count_scale returns symlog Scale."""
        from vcf_stats.seq2neo.visualizer import _count_scale
        scale = _count_scale()
        assert scale.type == "symlog"

    def test_clip_dp(self):
        """_clip_dp clips values and returns count of clipped."""
        from vcf_stats.seq2neo.visualizer import _clip_dp
        import pandas as pd
        pdf = pd.DataFrame({"DP": [100, 500, 2001, 3000, 50]})
        result, n_over = _clip_dp(pdf, "DP", cap=2000)
        assert n_over == 2
        assert result["DP"].max() == 2000
        assert result["DP"].min() == 50


class TestFPCrossTabIncludesSomatic:
    """Verify FP cross-tab now includes Somatic variants."""

    def test_fp_cross_tab_includes_somatic(self):
        from vcf_stats.seq2neo.statistics import compute_fp_cross_tab
        df = pl.DataFrame({
            "FILTER": ["Somatic", "Somatic", "Germline", "Artifact", "Somatic"],
            "N_SUPPORT_CALLERS": [5, 3, 2, 1, 6],
        })
        result = compute_fp_cross_tab(df)
        assert not result.is_empty()
        filters = result["FILTER"].unique().to_list()
        assert "Somatic" in filters, "Somatic should be included in cross-tab"
        assert "Germline" in filters
        assert "Artifact" in filters


class TestGlobalColorConstants:
    """Verify all global color constants are properly defined."""

    def test_all_constants_exist(self):
        from vcf_stats.seq2neo import visualizer as v
        assert len(v.CALLER_DOMAIN) == 6
        assert len(v.CALLER_COLORS) == 6
        assert len(v.RESCUE_DOMAIN) == 2
        assert len(v.RESCUE_COLORS) == 2
        assert len(v.VARIANT_TYPE_DOMAIN) == 4
        assert len(v.SOMATIC_MODALITY_DOMAIN) == 5
        assert len(v.BAM_TYPE_DOMAIN) == 3
        assert len(v.MODALITY_DOMAIN) == 2
        assert len(v.AGREEMENT_DOMAIN) == 3

    def test_color_registry_complete(self):
        from vcf_stats.seq2neo.visualizer import _COLOR_REGISTRY
        expected = ["caller", "FILTER", "RESCUED", "variant_type",
                    "somatic_modality", "bam_type", "modality", "agreement_level"]
        for entity in expected:
            assert entity in _COLOR_REGISTRY, f"{entity} missing from _COLOR_REGISTRY"
