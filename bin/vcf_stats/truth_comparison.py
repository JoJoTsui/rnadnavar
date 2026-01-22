#!/usr/bin/env python3
"""
Truth Set Comparison Module

Provides functionality for comparing variant caller VCFs against truth sets
to calculate precision, recall, and F1 metrics. Supports both direct bcftools
isec comparison and parsing of som.py output metrics.

Key Components:
    - TruthComparisonResult: Dataclass containing comparison metrics
    - TruthSetComparator: Main class for truth set comparison operations

Usage Example:
    >>> from vcf_stats.truth_comparison import TruthSetComparator, TruthComparisonResult
    >>>
    >>> # Initialize comparator with truth VCF and high-confidence regions
    >>> comparator = TruthSetComparator(
    ...     truth_vcf=Path("/path/to/truth.vcf.gz"),
    ...     high_confidence_bed=Path("/path/to/hc_regions.bed")
    ... )
    >>>
    >>> # Compare query VCF against truth set
    >>> result = comparator.compare(Path("/path/to/query.vcf.gz"))
    >>> print(f"Precision: {result.precision:.3f}")
    >>> print(f"Recall: {result.recall:.3f}")
    >>> print(f"F1 Score: {result.f1_score:.3f}")
    >>>
    >>> # Or parse existing som.py metrics
    >>> result = comparator.compare_from_som_py_metrics(
    ...     Path("/path/to/DS.metrics.json")
    ... )

Comparison Methodology:
    This module follows the same comparison methodology as Illumina hap.py som.py
    tool for somatic variant benchmarking:
    
    1. Filter query VCF to include only PASS (Somatic) variants
    2. Restrict both VCFs to high-confidence regions using BED file
    3. Use bcftools isec to compute intersection statistics
    4. Calculate TP, FP, FN based on variant overlap
    5. Compute Precision, Recall, and F1 Score

Metrics Definitions:
    - True Positives (TP): Variants present in both truth and query
    - False Positives (FP): Variants present in query but not in truth
    - False Negatives (FN): Variants present in truth but not in query
    - Precision: TP / (TP + FP) - proportion of called variants that are true
    - Recall: TP / (TP + FN) - proportion of true variants that were called
    - F1 Score: 2 * (Precision * Recall) / (Precision + Recall) - harmonic mean

Design Principles:
    - Flexible: Support both direct comparison and som.py metrics parsing
    - Robust: Handle edge cases (empty VCFs, zero denominators)
    - Consistent: Follow established benchmarking methodology
    - Informative: Provide detailed comparison results
"""

import json
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd


@dataclass
class TruthComparisonResult:
    """
    Result of truth set comparison.
    
    This dataclass encapsulates all metrics from comparing a query VCF
    against a truth set, including raw counts and calculated metrics.
    
    Attributes:
        tp: True Positives - variants present in both truth and query
        fp: False Positives - variants present in query but not in truth
        fn: False Negatives - variants present in truth but not in query
        precision: TP / (TP + FP) - proportion of called variants that are true
        recall: TP / (TP + FN) - proportion of true variants that were called
        f1_score: Harmonic mean of precision and recall
        query_total: Total number of variants in the query VCF (after filtering)
        truth_total: Total number of variants in the truth VCF (after filtering)
        source: Source of the comparison results ("bcftools_isec" or "som_py_metrics")
    
    Example:
        >>> result = TruthComparisonResult(
        ...     tp=1234,
        ...     fp=56,
        ...     fn=78,
        ...     precision=0.956,
        ...     recall=0.940,
        ...     f1_score=0.948,
        ...     query_total=1290,
        ...     truth_total=1312,
        ...     source="bcftools_isec"
        ... )
        >>> print(f"F1: {result.f1_score:.3f}")
        F1: 0.948
    """
    
    tp: int  # True Positives
    fp: int  # False Positives
    fn: int  # False Negatives
    precision: float
    recall: float
    f1_score: float
    query_total: int
    truth_total: int
    source: str  # "bcftools_isec" or "som_py_metrics"
    
    def __repr__(self) -> str:
        """Return a detailed string representation."""
        return (
            f"TruthComparisonResult("
            f"tp={self.tp}, fp={self.fp}, fn={self.fn}, "
            f"precision={self.precision:.4f}, recall={self.recall:.4f}, "
            f"f1_score={self.f1_score:.4f}, "
            f"query_total={self.query_total}, truth_total={self.truth_total}, "
            f"source='{self.source}')"
        )
    
    def to_dict(self) -> dict:
        """
        Convert result to dictionary format.
        
        Returns:
            Dictionary containing all comparison metrics
        """
        return {
            "tp": self.tp,
            "fp": self.fp,
            "fn": self.fn,
            "precision": self.precision,
            "recall": self.recall,
            "f1_score": self.f1_score,
            "query_total": self.query_total,
            "truth_total": self.truth_total,
            "source": self.source,
        }
    
    def to_series(self) -> pd.Series:
        """
        Convert result to pandas Series.
        
        Returns:
            pandas Series containing all comparison metrics
        """
        return pd.Series(self.to_dict())


class TruthSetComparator:
    """
    Compare variant caller VCFs against truth sets.
    
    This class provides functionality for comparing query VCFs against truth sets
    to calculate precision, recall, and F1 metrics. It supports both direct
    bcftools isec comparison and parsing of som.py output metrics.
    
    The comparison methodology follows Illumina hap.py som.py tool:
    1. Filter query VCF to PASS variants only
    2. Restrict to high-confidence regions
    3. Compute intersection using bcftools isec
    4. Calculate TP, FP, FN and derived metrics
    
    Attributes:
        truth_vcf: Path to the truth VCF file (optional)
        hc_bed: Path to high-confidence regions BED file (optional)
    
    Usage Example:
        >>> # Initialize with truth VCF and high-confidence regions
        >>> comparator = TruthSetComparator(
        ...     truth_vcf=Path("/path/to/truth.vcf.gz"),
        ...     high_confidence_bed=Path("/path/to/hc_regions.bed")
        ... )
        >>>
        >>> # Compare query VCF
        >>> result = comparator.compare(Path("/path/to/query.vcf.gz"))
        >>> print(f"Precision: {result.precision:.3f}")
        >>> print(f"Recall: {result.recall:.3f}")
        >>>
        >>> # Or parse som.py metrics directly
        >>> result = comparator.compare_from_som_py_metrics(
        ...     Path("/path/to/DS.metrics.json")
        ... )
    """
    
    def __init__(
        self,
        truth_vcf: Optional[Path] = None,
        high_confidence_bed: Optional[Path] = None
    ):
        """
        Initialize with optional truth VCF and high-confidence regions BED.
        
        Args:
            truth_vcf: Path to the truth VCF file containing validated variants.
                      Required for direct bcftools isec comparison.
            high_confidence_bed: Path to BED file defining high-confidence regions.
                                Used to restrict comparison to reliable regions.
        
        Raises:
            FileNotFoundError: If truth_vcf is provided but does not exist
            FileNotFoundError: If high_confidence_bed is provided but does not exist
        """
        self.truth_vcf = Path(truth_vcf) if truth_vcf else None
        self.hc_bed = Path(high_confidence_bed) if high_confidence_bed else None
        
        # Validate paths if provided
        if self.truth_vcf and not self.truth_vcf.exists():
            raise FileNotFoundError(f"Truth VCF not found: {self.truth_vcf}")
        if self.hc_bed and not self.hc_bed.exists():
            raise FileNotFoundError(f"High-confidence BED not found: {self.hc_bed}")
    
    def compare(
        self,
        query_vcf: Path,
        filter_pass_only: bool = True,
        som_py_metrics_path: Optional[Path] = None
    ) -> TruthComparisonResult:
        """
        Compare query VCF against truth set with fallback logic.
        
        This method implements a fallback strategy for comparison:
        
        1. First check for som.py metrics files (if som_py_metrics_path provided)
        2. If not found and bcftools tools available, use bcftools isec
        3. Handle edge case: F1 = 0 when precision + recall = 0
        
        The bcftools isec comparison follows the som.py methodology:
        1. Filter query VCF to PASS variants only (if filter_pass_only=True)
        2. Restrict both VCFs to high-confidence regions (if BED provided)
        3. Run bcftools isec to compute intersection
        4. Calculate TP, FP, FN and derived metrics
        
        Args:
            query_vcf: Path to the query VCF file to compare
            filter_pass_only: If True, filter query VCF to include only PASS
                            variants before comparison (default: True)
            som_py_metrics_path: Optional path to som.py metrics.json file.
                                If provided and exists, will use this instead
                                of running bcftools isec.
        
        Returns:
            TruthComparisonResult containing all comparison metrics
        
        Raises:
            FileNotFoundError: If query_vcf does not exist
            RuntimeError: If bcftools is not installed and no som.py metrics available
            ValueError: If truth_vcf was not provided and no som.py metrics available
            subprocess.CalledProcessError: If bcftools isec fails
        
        Example:
            >>> comparator = TruthSetComparator(
            ...     truth_vcf=Path("/path/to/truth.vcf.gz"),
            ...     high_confidence_bed=Path("/path/to/hc.bed")
            ... )
            >>> # Try som.py metrics first, fallback to bcftools
            >>> result = comparator.compare(
            ...     Path("/path/to/query.vcf.gz"),
            ...     som_py_metrics_path=Path("/path/to/DS.metrics.json")
            ... )
            >>> print(f"TP: {result.tp}, FP: {result.fp}, FN: {result.fn}")
        """
        # Validate query VCF exists
        query_vcf = Path(query_vcf)
        if not query_vcf.exists():
            raise FileNotFoundError(f"Query VCF not found: {query_vcf}")
        
        # FALLBACK STRATEGY:
        # Step 1: First check for som.py metrics files
        if som_py_metrics_path is not None:
            som_py_metrics_path = Path(som_py_metrics_path)
            if som_py_metrics_path.exists():
                return self.compare_from_som_py_metrics(som_py_metrics_path)
        
        # Step 2: If no som.py metrics, check if we can use bcftools isec
        if self.truth_vcf is None:
            raise ValueError(
                "Truth VCF must be provided during initialization for direct comparison. "
                "Alternatively, provide a som.py metrics file via som_py_metrics_path."
            )
        
        if not self._check_tools_available():
            raise RuntimeError(
                "bcftools is not installed or not accessible. "
                "Please install bcftools: conda install -c bioconda bcftools "
                "Or provide pre-computed som.py metrics via som_py_metrics_path."
            )
        
        # Step 3: Use bcftools isec for comparison
        return self._compare_with_bcftools(query_vcf, filter_pass_only)
    
    def _compare_with_bcftools(
        self,
        query_vcf: Path,
        filter_pass_only: bool = True
    ) -> TruthComparisonResult:
        """
        Perform comparison using bcftools isec.
        
        This is the internal method that performs the actual bcftools isec
        comparison. It is called by compare() when som.py metrics are not
        available.
        
        Args:
            query_vcf: Path to the query VCF file to compare
            filter_pass_only: If True, filter query VCF to include only PASS
                            variants before comparison
        
        Returns:
            TruthComparisonResult containing all comparison metrics
        
        Raises:
            subprocess.CalledProcessError: If bcftools isec fails
        """
        # Create temporary directory for intermediate files
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Step 1: Filter query VCF to PASS variants if requested
            if filter_pass_only:
                filtered_query = self._filter_to_pass_variants(query_vcf, tmpdir)
            else:
                filtered_query = query_vcf
            
            # Step 2: Restrict to high-confidence regions if BED provided
            if self.hc_bed:
                restricted_query = self._restrict_to_regions(filtered_query, tmpdir, "query")
                restricted_truth = self._restrict_to_regions(self.truth_vcf, tmpdir, "truth")
            else:
                restricted_query = filtered_query
                restricted_truth = self.truth_vcf
            
            # Step 3: Run bcftools isec
            tp, fp, fn = self._run_bcftools_isec(restricted_query, restricted_truth, tmpdir)
            
            # Step 4: Calculate metrics (handles F1 = 0 when precision + recall = 0)
            precision, recall, f1_score = self._calculate_metrics(tp, fp, fn)
            
            # Calculate totals
            query_total = tp + fp
            truth_total = tp + fn
            
            return TruthComparisonResult(
                tp=tp,
                fp=fp,
                fn=fn,
                precision=precision,
                recall=recall,
                f1_score=f1_score,
                query_total=query_total,
                truth_total=truth_total,
                source="bcftools_isec"
            )
    
    def compare_from_som_py_metrics(
        self,
        metrics_json_path: Path
    ) -> TruthComparisonResult:
        """
        Parse som.py metrics.json output file to extract comparison results.
        
        This method parses pre-computed metrics from som.py output files,
        which follow a standard JSON format containing TP, FP, FN counts
        and calculated metrics.
        
        Expected file format: {caller}.metrics.json
        Example paths:
            - DS.metrics.json (DeepSomatic)
            - M2.metrics.json (Mutect2)
            - S2.metrics.json (Strelka2)
        
        Expected JSON structure:
            {
                "type": "somatic",
                "tp": 1234,
                "fp": 56,
                "fn": 78,
                "precision": 0.956,
                "recall": 0.940,
                "f1": 0.948
            }
        
        Args:
            metrics_json_path: Path to the som.py metrics.json file
        
        Returns:
            TruthComparisonResult containing parsed metrics
        
        Raises:
            FileNotFoundError: If metrics_json_path does not exist
            json.JSONDecodeError: If file is not valid JSON
            KeyError: If required fields are missing from JSON
        
        Example:
            >>> comparator = TruthSetComparator()
            >>> result = comparator.compare_from_som_py_metrics(
            ...     Path("/path/to/DS.metrics.json")
            ... )
            >>> print(f"F1: {result.f1_score:.3f}")
        """
        metrics_json_path = Path(metrics_json_path)
        if not metrics_json_path.exists():
            raise FileNotFoundError(f"Metrics JSON not found: {metrics_json_path}")
        
        with open(metrics_json_path, "r") as f:
            metrics = json.load(f)
        
        # Extract required fields
        tp = int(metrics["tp"])
        fp = int(metrics["fp"])
        fn = int(metrics["fn"])
        
        # Use provided metrics or calculate if not present
        precision = float(metrics.get("precision", 0.0))
        recall = float(metrics.get("recall", 0.0))
        f1_score = float(metrics.get("f1", metrics.get("f1_score", 0.0)))
        
        # If metrics not provided, calculate them
        if precision == 0.0 and recall == 0.0 and (tp > 0 or fp > 0 or fn > 0):
            precision, recall, f1_score = self._calculate_metrics(tp, fp, fn)
        
        # Calculate totals
        query_total = tp + fp
        truth_total = tp + fn
        
        return TruthComparisonResult(
            tp=tp,
            fp=fp,
            fn=fn,
            precision=precision,
            recall=recall,
            f1_score=f1_score,
            query_total=query_total,
            truth_total=truth_total,
            source="som_py_metrics"
        )
    
    def parse_som_py_stats_csv(self, stats_csv_path: Path) -> pd.DataFrame:
        """
        Parse som.py stats.csv output file for detailed statistics.
        
        This method parses the detailed statistics CSV file produced by som.py,
        which contains per-variant-type breakdowns and additional metrics.
        
        Expected file format: {caller}.stats.csv
        
        Args:
            stats_csv_path: Path to the som.py stats.csv file
        
        Returns:
            pandas DataFrame with parsed statistics
        
        Raises:
            FileNotFoundError: If stats_csv_path does not exist
            pd.errors.ParserError: If file cannot be parsed as CSV
        
        Example:
            >>> comparator = TruthSetComparator()
            >>> stats_df = comparator.parse_som_py_stats_csv(
            ...     Path("/path/to/DS.stats.csv")
            ... )
            >>> print(stats_df.columns.tolist())
        """
        stats_csv_path = Path(stats_csv_path)
        if not stats_csv_path.exists():
            raise FileNotFoundError(f"Stats CSV not found: {stats_csv_path}")
        
        return pd.read_csv(stats_csv_path)
    
    def _check_tools_available(self) -> bool:
        """
        Check if bcftools is installed and accessible.
        
        Returns:
            True if bcftools is available, False otherwise
        """
        return shutil.which("bcftools") is not None
    
    def _filter_to_pass_variants(
        self,
        vcf_path: Path,
        output_dir: Path
    ) -> Path:
        """
        Filter VCF to include only PASS variants.
        
        Uses bcftools view to filter the VCF to include only variants
        with FILTER=PASS.
        
        Args:
            vcf_path: Path to input VCF file
            output_dir: Directory for output file
        
        Returns:
            Path to filtered VCF file
        
        Raises:
            subprocess.CalledProcessError: If bcftools view fails
        """
        output_path = output_dir / f"{vcf_path.stem}.pass.vcf.gz"
        
        cmd = [
            "bcftools", "view",
            "-f", "PASS",
            "-Oz",
            "-o", str(output_path),
            str(vcf_path)
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
        
        # Index the output
        subprocess.run(
            ["bcftools", "index", "-t", str(output_path)],
            check=True,
            capture_output=True
        )
        
        return output_path
    
    def _restrict_to_regions(
        self,
        vcf_path: Path,
        output_dir: Path,
        prefix: str
    ) -> Path:
        """
        Restrict VCF to high-confidence regions using BED file.
        
        Uses bcftools view with -R option to restrict variants to
        regions defined in the high-confidence BED file.
        
        Args:
            vcf_path: Path to input VCF file
            output_dir: Directory for output file
            prefix: Prefix for output filename
        
        Returns:
            Path to restricted VCF file
        
        Raises:
            subprocess.CalledProcessError: If bcftools view fails
        """
        output_path = output_dir / f"{prefix}.restricted.vcf.gz"
        
        cmd = [
            "bcftools", "view",
            "-R", str(self.hc_bed),
            "-Oz",
            "-o", str(output_path),
            str(vcf_path)
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
        
        # Index the output
        subprocess.run(
            ["bcftools", "index", "-t", str(output_path)],
            check=True,
            capture_output=True
        )
        
        return output_path
    
    def _run_bcftools_isec(
        self,
        query_vcf: Path,
        truth_vcf: Path,
        output_dir: Path
    ) -> Tuple[int, int, int]:
        """
        Run bcftools isec and return (TP, FP, FN) counts.
        
        Uses bcftools isec to compute the intersection between query and
        truth VCFs. The output directories contain:
        - 0000.vcf: Variants unique to query (FP)
        - 0001.vcf: Variants unique to truth (FN)
        - 0002.vcf: Variants in both, from query (TP)
        - 0003.vcf: Variants in both, from truth (TP verification)
        
        Args:
            query_vcf: Path to query VCF file
            truth_vcf: Path to truth VCF file
            output_dir: Directory for bcftools isec output
        
        Returns:
            Tuple of (TP, FP, FN) counts
        
        Raises:
            subprocess.CalledProcessError: If bcftools isec fails
        """
        isec_dir = output_dir / "isec_output"
        
        cmd = [
            "bcftools", "isec",
            "-p", str(isec_dir),
            str(query_vcf),
            str(truth_vcf)
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
        
        # Count variants in each output file
        fp = self._count_variants(isec_dir / "0000.vcf")  # Query only
        fn = self._count_variants(isec_dir / "0001.vcf")  # Truth only
        tp = self._count_variants(isec_dir / "0002.vcf")  # Both (from query)
        
        return tp, fp, fn
    
    def _count_variants(self, vcf_path: Path) -> int:
        """
        Count the number of variants in a VCF file.
        
        Uses bcftools view to count non-header lines.
        
        Args:
            vcf_path: Path to VCF file
        
        Returns:
            Number of variants in the VCF
        """
        if not vcf_path.exists():
            return 0
        
        result = subprocess.run(
            ["bcftools", "view", "-H", str(vcf_path)],
            capture_output=True,
            text=True
        )
        
        # Count non-empty lines
        lines = [line for line in result.stdout.strip().split("\n") if line]
        return len(lines)
    
    def _calculate_metrics(
        self,
        tp: int,
        fp: int,
        fn: int
    ) -> Tuple[float, float, float]:
        """
        Calculate Precision, Recall, F1 from TP/FP/FN.
        
        Handles edge cases where denominators are zero:
        - If TP + FP = 0: Precision = 0.0
        - If TP + FN = 0: Recall = 0.0
        - If Precision + Recall = 0: F1 = 0.0
        
        Args:
            tp: True Positives count
            fp: False Positives count
            fn: False Negatives count
        
        Returns:
            Tuple of (precision, recall, f1_score)
        
        Example:
            >>> comparator = TruthSetComparator()
            >>> p, r, f1 = comparator._calculate_metrics(100, 10, 5)
            >>> print(f"P={p:.3f}, R={r:.3f}, F1={f1:.3f}")
            P=0.909, R=0.952, F1=0.930
        """
        # Calculate Precision: TP / (TP + FP)
        if tp + fp > 0:
            precision = tp / (tp + fp)
        else:
            precision = 0.0
        
        # Calculate Recall: TP / (TP + FN)
        if tp + fn > 0:
            recall = tp / (tp + fn)
        else:
            recall = 0.0
        
        # Calculate F1 Score: 2 * (P * R) / (P + R)
        if precision + recall > 0:
            f1_score = 2 * (precision * recall) / (precision + recall)
        else:
            f1_score = 0.0
        
        return precision, recall, f1_score


# Module-level convenience functions

def compare_vcfs(
    query_vcf: Path,
    truth_vcf: Path,
    high_confidence_bed: Optional[Path] = None,
    filter_pass_only: bool = True
) -> TruthComparisonResult:
    """
    Convenience function to compare a query VCF against a truth set.
    
    Args:
        query_vcf: Path to the query VCF file
        truth_vcf: Path to the truth VCF file
        high_confidence_bed: Optional path to high-confidence regions BED
        filter_pass_only: If True, filter query to PASS variants only
    
    Returns:
        TruthComparisonResult containing comparison metrics
    
    Example:
        >>> result = compare_vcfs(
        ...     query_vcf=Path("/path/to/query.vcf.gz"),
        ...     truth_vcf=Path("/path/to/truth.vcf.gz"),
        ...     high_confidence_bed=Path("/path/to/hc.bed")
        ... )
        >>> print(f"F1: {result.f1_score:.3f}")
    """
    comparator = TruthSetComparator(
        truth_vcf=truth_vcf,
        high_confidence_bed=high_confidence_bed
    )
    return comparator.compare(query_vcf, filter_pass_only=filter_pass_only)


def parse_som_py_metrics(metrics_json_path: Path) -> TruthComparisonResult:
    """
    Convenience function to parse som.py metrics JSON file.
    
    Args:
        metrics_json_path: Path to the som.py metrics.json file
    
    Returns:
        TruthComparisonResult containing parsed metrics
    
    Example:
        >>> result = parse_som_py_metrics(Path("/path/to/DS.metrics.json"))
        >>> print(f"TP: {result.tp}, FP: {result.fp}, FN: {result.fn}")
    """
    comparator = TruthSetComparator()
    return comparator.compare_from_som_py_metrics(metrics_json_path)


print("✓ Truth Set Comparison module imported successfully")
