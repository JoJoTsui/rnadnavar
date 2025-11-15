"""
Filtering functions for VCF variant processing.

This module provides reusable functions for normalizing and categorizing
variant filter values from different variant callers. It handles the diverse
filter naming conventions used by different callers (Mutect2, Strelka,
DeepSomatic) and maps them to standardized filter names and categories.

Constants:
    FILTER_MAP: Dictionary mapping caller-specific filter names to normalized names
    FILTER_CATEGORIES: Dictionary grouping normalized filters into major categories

Functions:
    normalize_filter: Normalize filter strings using the mapping
    categorize_filter: Categorize normalized filters into major groups

Example:
    >>> from vcf_utils.filters import normalize_filter, categorize_filter
    >>> 
    >>> # Normalize different caller filter formats
    >>> mutect2_filter = "germline;base_qual"
    >>> normalized = normalize_filter(mutect2_filter)
    >>> print(normalized)
    Germline;LowBaseQuality
    >>> 
    >>> # Categorize the normalized filter
    >>> category = categorize_filter(normalized)
    >>> print(category)
    germline;technical
    >>> 
    >>> # Handle Strelka filters
    >>> strelka_filter = "LowDepth;LowEVS"
    >>> print(normalize_filter(strelka_filter))
    LowDepth;LowEvidenceScore
"""
import re


# Filter normalization mapping - based on actual caller outputs
FILTER_MAP = {
    # DeepSomatic filters
    'GERMLINE': 'Germline',
    'RefCall': 'ReferenceCall',
    
    # Strelka filters
    'LowDepth': 'LowDepth',
    'LowEVS': 'LowEvidenceScore',  # EVS = Empirical Variant Score
    'LowDepth;LowEVS': 'LowDepth;LowEvidenceScore',
    
    # Mutect2 common filters (add as encountered)
    'germline': 'Germline',
    'panel_of_normals': 'PanelOfNormals',
    'base_qual': 'LowBaseQuality',
    'strand_bias': 'StrandBias',
    'weak_evidence': 'WeakEvidence',
    'clustered_events': 'ClusteredEvents',
    'multiallelic': 'Multiallelic',
    'fragment': 'FragmentEvidence',
    'haplotype': 'Haplotype',
    'slippage': 'Slippage',
    'orientation': 'OrientationBias',
    't_lod': 'LowTumorLOD',
    'normal_artifact': 'NormalArtifact',
    'contamination': 'Contamination',
    'duplicate': 'Duplicate',
    'position': 'PositionFilter',
    'strict_strand': 'StrictStrandBias',
    'n_ratio': 'NormalReadRatio',
    'map_qual': 'LowMappingQuality',
    'germline_risk': 'GermlineRisk',
    
    # Generic mappings
    'LowQual': 'LowQuality',
    'low_qual': 'LowQuality',
    'LQ': 'LowQuality',
    'QualityFilter': 'LowQuality',
    'low_depth': 'LowDepth',
    'DP': 'LowDepth',
    'MinDepth': 'LowDepth',
    'min_depth': 'LowDepth',
    'StrandBias': 'StrandBias',
    'SB': 'StrandBias',
    'strand_bias': 'StrandBias',
    'STRAND_BIAS': 'StrandBias',
    'LowMQ': 'LowMappingQuality',
    'MQ': 'LowMappingQuality',
    'low_mq': 'LowMappingQuality',
    'LowAF': 'LowAlleleFrequency',
    'low_af': 'LowAlleleFrequency',
    'MinAF': 'LowAlleleFrequency',
    'clustered_events': 'ClusteredEvents',
    'ClusteredEvents': 'ClusteredEvents',
    'LowBaseQual': 'LowBaseQuality',
    'base_qual': 'LowBaseQuality',
}

# Filter categories for unified classification
FILTER_CATEGORIES = {
    'quality': ['LowQuality', 'LowEvidenceScore', 'LowTumorLOD', 'WeakEvidence'],
    'depth': ['LowDepth', 'LowAlleleFrequency'],
    'bias': ['StrandBias', 'OrientationBias', 'StrictStrandBias'],
    'germline': ['Germline', 'GermlineRisk'],
    'artifact': ['NormalArtifact', 'PanelOfNormals', 'Contamination'],
    'technical': ['LowMappingQuality', 'LowBaseQuality', 'Duplicate', 'FragmentEvidence'],
    'reference': ['ReferenceCall'],
}


def normalize_filter(filter_str):
    """
    Normalize filter string using mapping.
    
    This function standardizes filter names from different variant callers
    into a consistent naming scheme. It handles PASS, missing values, and
    compound filters (semicolon, comma, or pipe-separated).
    
    Args:
        filter_str (str or None): Filter string from VCF FILTER field. Can be:
            - 'PASS': Variant passed all filters
            - '.': Missing filter value
            - None or empty: Missing filter value
            - Single filter: e.g., 'germline', 'LowDepth'
            - Compound filter: e.g., 'LowDepth;LowEVS', 'germline,base_qual'
    
    Returns:
        str: Normalized filter string. Returns 'PASS' if the variant passed
            all filters or if no filters are specified. Otherwise returns
            semicolon-separated normalized filter names.
    
    Notes:
        - Unknown filters are kept but sanitized (spaces replaced with underscores)
        - Compound filters are split on semicolon, comma, or pipe characters
        - Filter names are mapped using the FILTER_MAP constant
        - Empty filter lists after normalization return 'PASS'
    
    Example:
        >>> # Mutect2 filters
        >>> normalize_filter('germline')
        'Germline'
        >>> normalize_filter('germline;base_qual;strand_bias')
        'Germline;LowBaseQuality;StrandBias'
        >>> 
        >>> # Strelka filters
        >>> normalize_filter('LowDepth;LowEVS')
        'LowDepth;LowEvidenceScore'
        >>> 
        >>> # DeepSomatic filters
        >>> normalize_filter('GERMLINE')
        'Germline'
        >>> normalize_filter('RefCall')
        'ReferenceCall'
        >>> 
        >>> # PASS and missing values
        >>> normalize_filter('PASS')
        'PASS'
        >>> normalize_filter('.')
        'PASS'
        >>> normalize_filter(None)
        'PASS'
    """
    # Handle missing, empty, or PASS
    if not filter_str or filter_str == 'PASS' or filter_str == '.':
        return 'PASS'
    
    # Split compound filters (e.g., "LowDepth;LowEVS")
    filters = re.split('[;,|]', filter_str)
    normalized = []
    
    for f in filters:
        f = f.strip()
        if not f or f == 'PASS' or f == '.':
            continue
        
        # Apply mapping if exists, otherwise keep original
        if f in FILTER_MAP:
            normalized.append(FILTER_MAP[f])
        else:
            # Keep unknown filter but sanitize it
            normalized.append(f.replace(' ', '_'))
    
    return ';'.join(normalized) if normalized else 'PASS'


def categorize_filter(normalized_filter):
    """
    Categorize normalized filter into major categories.
    
    This function groups normalized filter names into high-level categories
    for easier analysis and reporting. Categories include quality, depth,
    bias, germline, artifact, technical, and reference.
    
    Args:
        normalized_filter (str): Normalized filter string from normalize_filter().
            Can be 'PASS' or semicolon-separated normalized filter names.
    
    Returns:
        str: Category string. Possible values:
            - 'PASS': Variant passed all filters
            - 'quality': Low quality or evidence score filters
            - 'depth': Low depth or allele frequency filters
            - 'bias': Strand bias or orientation bias filters
            - 'germline': Germline variant filters
            - 'artifact': Artifact or normal contamination filters
            - 'technical': Technical quality filters (mapping, base quality, etc.)
            - 'reference': Reference call filters
            - 'Other': Filter doesn't match any known category
            - Multiple categories separated by semicolon if compound filter
              spans multiple categories
    
    Notes:
        - Uses FILTER_CATEGORIES constant for category membership
        - Returns 'PASS' unchanged
        - Returns 'Other' for unknown filters
        - For compound filters with multiple categories, returns semicolon-separated
          sorted category names
    
    Example:
        >>> # Single category filters
        >>> categorize_filter('PASS')
        'PASS'
        >>> categorize_filter('Germline')
        'germline'
        >>> categorize_filter('LowDepth')
        'depth'
        >>> categorize_filter('StrandBias')
        'bias'
        >>> 
        >>> # Compound filters
        >>> categorize_filter('Germline;LowBaseQuality')
        'germline;technical'
        >>> categorize_filter('LowDepth;LowEvidenceScore')
        'depth;quality'
        >>> 
        >>> # Unknown filters
        >>> categorize_filter('UnknownFilter')
        'Other'
    """
    if normalized_filter == 'PASS':
        return 'PASS'
    
    filters = normalized_filter.split(';')
    categories = set()
    
    for filt in filters:
        for cat, members in FILTER_CATEGORIES.items():
            if filt in members:
                categories.add(cat)
                break
    
    if not categories:
        return 'Other'
    elif len(categories) == 1:
        return list(categories)[0]
    else:
        # Multiple categories - return concatenated
        return ';'.join(sorted(categories))
