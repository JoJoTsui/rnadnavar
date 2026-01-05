#!/usr/bin/env python3
"""
Chromosome Filtering Utilities

This module provides utilities for filtering variants to canonical chromosomes
(1-22, X, Y, M/MT) and managing chromosome naming conventions. It handles
different chromosome naming schemes (with/without 'chr' prefix) and provides
functions for both variant-level and VCF header-level filtering.

Key Features:
- Canonical chromosome detection
- Chromosome name normalization
- VCF header contig filtering
- Integration with variant reading pipelines
"""

from .classification_config import CANONICAL_SET
from .io_utils import normalize_chromosome


def is_canonical_chromosome(chrom):
    """
    Check if chromosome is canonical (1-22, X, Y, M/MT).

    Args:
        chrom (str): Chromosome name (e.g., 'chr1', '1', 'chrX', 'chrM')

    Returns:
        bool: True if chromosome is canonical

    Example:
        >>> is_canonical_chromosome('chr1')
        True
        >>> is_canonical_chromosome('1')
        True
        >>> is_canonical_chromosome('chrX')
        True
        >>> is_canonical_chromosome('chrM')
        True
        >>> is_canonical_chromosome('chrUn_gl000220')
        False
    """
    chrom_norm = normalize_chromosome(chrom)

    # Handle MT/M equivalence
    if chrom_norm == "MT":
        chrom_norm = "M"

    return chrom_norm in CANONICAL_SET


def get_canonical_contig_names(header):
    """
    Extract canonical contig names from VCF header.

    Args:
        header: pysam.VariantHeader or cyvcf2.VCF header

    Returns:
        list: List of canonical contig names as they appear in the header
    """
    canonical_contigs = []

    # Handle pysam header
    if hasattr(header, "contigs"):
        for contig in header.contigs:
            if is_canonical_chromosome(contig.name):
                canonical_contigs.append(contig.name)
    # Handle cyvcf2 header
    elif hasattr(header, "seqnames"):
        for contig in header.seqnames:
            if is_canonical_chromosome(contig):
                canonical_contigs.append(contig)

    return canonical_contigs


def filter_vcf_header_contigs(input_header, output_header, keep_non_canonical=False):
    """
    Filter VCF header to only include canonical chromosomes.

    This function modifies the output header to only include contigs for
    canonical chromosomes when keep_non_canonical is False.

    Args:
        input_header: Input VCF header (pysam.VariantHeader)
        output_header: Output VCF header to modify (pysam.VariantHeader)
        keep_non_canonical (bool): If True, keep all contigs

    Returns:
        output_header: Modified header with filtered contigs
    """
    if keep_non_canonical:
        return output_header

    # Clear existing contigs
    for contig in list(output_header.contigs):
        output_header.contigs.remove_header(contig)

    # Add only canonical contigs from input
    for contig in input_header.contigs:
        if is_canonical_chromosome(contig.name):
            # Copy contig info
            if hasattr(contig, "length") and contig.length:
                output_header.contigs.add(contig.name, length=contig.length)
            else:
                output_header.contigs.add(contig.name)

    return output_header


def should_skip_variant(variant, include_non_canonical=False, use_cyvcf2=False):
    """
    Check if a variant should be skipped based on chromosome filtering.

    Args:
        variant: Variant object (cyvcf2.Variant or pysam.VariantRecord)
        include_non_canonical (bool): If True, include all chromosomes
        use_cyvcf2 (bool): If True, variant is cyvcf2.Variant

    Returns:
        bool: True if variant should be skipped
    """
    if include_non_canonical:
        return False

    if use_cyvcf2:
        chrom = variant.CHROM
    else:
        chrom = variant.chrom

    return not is_canonical_chromosome(chrom)


def get_canonical_chromosome_list(include_chr_prefix=False):
    """
    Get list of all canonical chromosome names.

    Args:
        include_chr_prefix (bool): If True, include 'chr' prefix

    Returns:
        list: List of canonical chromosome names
    """
    chromosomes = []

    # Autosomes 1-22
    for i in range(1, 23):
        chrom = f"chr{i}" if include_chr_prefix else str(i)
        chromosomes.append(chrom)

    # Sex chromosomes
    for sex in ["X", "Y"]:
        chrom = f"chr{sex}" if include_chr_prefix else sex
        chromosomes.append(chrom)

    # Mitochondrial
    mt = "chrM" if include_chr_prefix else "M"
    chromosomes.append(mt)

    return chromosomes


def standardize_chromosome_name(chrom, add_chr_prefix=False):
    """
    Standardize chromosome name to consistent format.

    Args:
        chrom (str): Input chromosome name
        add_chr_prefix (bool): If True, add 'chr' prefix to output

    Returns:
        str: Standardized chromosome name

    Example:
        >>> standardize_chromosome_name('chr1', add_chr_prefix=False)
        '1'
        >>> standardize_chromosome_name('1', add_chr_prefix=True)
        'chr1'
        >>> standardize_chromosome_name('chrMT', add_chr_prefix=False)
        'M'
    """
    # First normalize (remove chr prefix)
    norm = normalize_chromosome(chrom)

    # Handle MT -> M conversion
    if norm == "MT":
        norm = "M"

    # Add prefix if requested
    if add_chr_prefix and norm in CANONICAL_SET:
        return f"chr{norm}"
    else:
        return norm


def get_chromosome_sort_key(chrom):
    """
    Get sort key for chromosome to enable proper ordering.

    Args:
        chrom (str): Chromosome name

    Returns:
        tuple: Sort key (type, value) where type 0=autosome, 1=sex, 2=mito
    """
    norm = normalize_chromosome(chrom)

    # Try to convert to integer (autosomes)
    try:
        return (0, int(norm))
    except ValueError:
        pass

    # Sex chromosomes
    if norm == "X":
        return (1, 0)
    elif norm == "Y":
        return (1, 1)
    # Mitochondrial
    elif norm in ["M", "MT"]:
        return (2, 0)
    # Non-canonical (sort last)
    else:
        return (3, norm)
