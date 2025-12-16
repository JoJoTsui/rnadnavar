"""
I/O utility functions for VCF file operations.

This module provides reusable functions for reading and writing VCF files,
creating VCF headers with custom INFO fields, and handling variant identifiers.
It supports both consensus and rescue workflows with appropriate metadata fields.

Functions:
    get_caller_name: Extract caller name from VCF filename
    normalize_chromosome: Normalize chromosome names for consistent sorting
    variant_key: Create unique variant identifier
    is_snv: Determine if variant is a single nucleotide variant
    create_output_header: Create VCF header with consensus/rescue INFO fields
    write_union_vcf: Write aggregated variants to VCF file

Example:
    >>> from vcf_utils.io_utils import create_output_header, write_union_vcf
    >>> from cyvcf2 import VCF
    >>> 
    >>> # Read template header
    >>> template_vcf = VCF('input.vcf.gz')
    >>> template_header = template_vcf
    >>> 
    >>> # Create output header for consensus
    >>> header = create_output_header(template_header, 'sample123', include_rescue_fields=False)
    >>> 
    >>> # Write aggregated variants
    >>> write_union_vcf(aggregated_data, template_header, 'sample123', 
    ...                 'output.vcf.gz', 'vcf.gz', ['mutect2', 'strelka'])
"""
from pathlib import Path


def get_caller_name(filename):
    """
    Extract caller name from VCF filename.
    
    This function parses VCF filenames to extract the variant caller name,
    handling common naming patterns used in variant calling pipelines.
    
    Args:
        filename (str): VCF filename (can include path)
    
    Returns:
        str: Caller name extracted from filename
    
    Notes:
        - Handles '.variants.' pattern: sample.strelka.variants.vcf.gz -> 'strelka'
        - Handles '.consensus.' pattern: sample.consensus.vcf.gz -> 'consensus'
        - Fallback: takes second dot-separated component
    
    Example:
        >>> get_caller_name('sample123.mutect2.variants.vcf.gz')
        'mutect2'
        >>> get_caller_name('sample123.strelka.variants.vcf.gz')
        'strelka'
        >>> get_caller_name('sample123.consensus.vcf.gz')
        'consensus'
    """
    name = Path(filename).name
    # Handle .variants. pattern (e.g., sample.strelka.variants.vcf.gz)
    if '.variants.' in name:
        parts = name.split('.')
        for i, part in enumerate(parts):
            if part == 'variants' and i > 0:
                return parts[i-1]
    # Handle .consensus. pattern
    if '.consensus.' in name:
        return 'consensus'
    # Fallback: take second component
    parts = name.split('.')
    return parts[1] if len(parts) > 1 else parts[0]


def normalize_chromosome(chrom):
    """
    Normalize chromosome names for consistent sorting.
    
    This function removes 'chr' or 'Chr' prefixes from chromosome names
    to enable consistent sorting and comparison across different reference
    genome naming conventions.
    
    Args:
        chrom (str or int): Chromosome name (e.g., 'chr1', 'Chr1', '1', 'chrX')
    
    Returns:
        str: Normalized chromosome name without prefix (e.g., '1', 'X', 'MT')
    
    Example:
        >>> normalize_chromosome('chr1')
        '1'
        >>> normalize_chromosome('Chr22')
        '22'
        >>> normalize_chromosome('chrX')
        'X'
        >>> normalize_chromosome('1')
        '1'
    """
    chrom = str(chrom).replace('chr', '').replace('Chr', '')
    return chrom


def variant_key(variant, use_cyvcf2=True):
    """
    Create unique key for variant including all ALT alleles.
    
    This function generates a unique identifier for a variant based on its
    genomic position and alleles. It works with both cyvcf2.Variant objects
    and dictionary representations.
    
    Args:
        variant (cyvcf2.Variant or dict): Variant object or dictionary with
            CHROM, POS, REF, and ALT fields
        use_cyvcf2 (bool, optional): If True, treat variant as cyvcf2.Variant
            object. If False, treat as dictionary. Default: True
    
    Returns:
        str: Unique variant key in format "chrom:pos:ref:alt" where:
            - chrom is normalized (no 'chr' prefix)
            - pos is 1-based position
            - ref is reference allele
            - alt is comma-separated alternate alleles
    
    Example:
        >>> from cyvcf2 import VCF
        >>> vcf = VCF('variants.vcf.gz')
        >>> for variant in vcf:
        ...     key = variant_key(variant, use_cyvcf2=True)
        ...     print(key)
        ...     break
        1:12345:A:G
        >>> 
        >>> # With dictionary
        >>> var_dict = {'CHROM': 'chr2', 'POS': 67890, 'REF': 'C', 'ALT': 'T,G'}
        >>> variant_key(var_dict, use_cyvcf2=False)
        '2:67890:C:T,G'
    """
    if use_cyvcf2:
        # cyvcf2.Variant object
        chrom = normalize_chromosome(variant.CHROM)
        pos = variant.POS
        ref = variant.REF
        # Handle multiple ALT alleles
        alts = ','.join(variant.ALT) if variant.ALT else '.'
    else:
        # dict representation
        chrom = normalize_chromosome(variant['CHROM'])
        pos = variant['POS']
        ref = variant['REF']
        alts = variant['ALT']
    
    return f"{chrom}:{pos}:{ref}:{alts}"


def is_snv(ref, alt_list):
    """
    Determine if variant is a single nucleotide variant (SNV).
    
    A variant is classified as SNV if both the reference and all alternate
    alleles are single nucleotides. Multi-nucleotide variants, insertions,
    deletions, and complex variants are not SNVs.
    
    Args:
        ref (str): Reference allele
        alt_list (list): List of alternate alleles
    
    Returns:
        bool: True if variant is SNV, False otherwise
    
    Notes:
        - Returns False if ref is not a single base
        - Returns False if any alt allele is not a single base
        - Returns False if any alt allele is '*' (deletion indicator)
    
    Example:
        >>> is_snv('A', ['G'])
        True
        >>> is_snv('A', ['G', 'T'])
        True
        >>> is_snv('AT', ['A'])
        False
        >>> is_snv('A', ['AT'])
        False
        >>> is_snv('A', ['*'])
        False
    """
    if len(ref) != 1:
        return False
    for alt in alt_list:
        if len(alt) != 1 or alt == '*':  # * indicates deletion in some formats
            return False
    return True



def create_output_header(template_header, sample_name, include_rescue_fields=False):
    """
    Create output VCF header with all INFO fields.
    
    This function creates a pysam VariantHeader with all necessary INFO and
    FILTER fields for consensus or rescue VCF output. It preserves metadata
    from the template header and adds custom fields for variant aggregation.
    
    Args:
        template_header (cyvcf2.VCF): cyvcf2 VCF header object to use as template.
            Metadata lines (contigs, references, etc.) are preserved.
        sample_name (str): Sample name for the output VCF. If None or empty,
            defaults to 'SAMPLE'.
        include_rescue_fields (bool, optional): If True, add rescue-specific
            INFO fields for modality tracking (MODALITIES, CALLERS_BY_MODALITY,
            DNA_SUPPORT, RNA_SUPPORT, CROSS_MODALITY, DP_DNA_MEAN, DP_RNA_MEAN,
            VAF_DNA_MEAN, VAF_RNA_MEAN). Default: False
    
    Returns:
        pysam.VariantHeader: VCF header object with all necessary INFO and
            FILTER fields defined. Ready for use with pysam.VariantFile.
    
    INFO Fields Added (Consensus):
        - N_CALLERS: Total number of aggregated callers
        - CALLERS: List of all aggregated callers (pipe-separated)
        - N_SUPPORT_CALLERS: Number of callers that detected this variant
        - CALLERS_SUPPORT: Callers that detected this variant (pipe-separated)
        - FILTERS_ORIGINAL: Original filter values from each caller
        - FILTERS_NORMALIZED: Normalized filter categories
        - FILTERS_CATEGORY: Filter categories
        - UNIFIED_FILTER: Unified filter status based on majority
        - PASSES_CONSENSUS: Whether variant passes consensus threshold (YES/NO)
        - RESCUED: Variant included via cross-modality consensus (YES/NO)
        - QUAL_MEAN/MIN/MAX: Quality score statistics
        - CONSENSUS_GT: Consensus genotype across callers
        - GT_BY_CALLER: Genotypes from each caller
        - DP_MEAN/MIN/MAX: Depth statistics
        - DP_BY_CALLER: Depth values from each caller
        - VAF_MEAN/MIN/MAX: VAF statistics
        - VAF_BY_CALLER: VAF values from each caller
    
    INFO Fields Added (Rescue, when include_rescue_fields=True):
        - MODALITIES: Modalities where variant was detected
        - CALLERS_BY_MODALITY: Callers grouped by modality
        - DNA_SUPPORT: Number of DNA callers supporting this variant
        - RNA_SUPPORT: Number of RNA callers supporting this variant
        - CROSS_MODALITY: Whether variant has cross-modality support (YES/NO)
        - DP_DNA_MEAN: Mean depth across DNA callers
        - DP_RNA_MEAN: Mean depth across RNA callers
        - VAF_DNA_MEAN: Mean VAF across DNA callers
        - VAF_RNA_MEAN: Mean VAF across RNA callers
    
    FILTER Fields Added:
        - LowQuality, LowDepth, StrandBias, Germline, Artifact, NoConsensus,
          LowEvidenceScore, ReferenceCall
    
    Example:
        >>> from cyvcf2 import VCF
        >>> import pysam
        >>> 
        >>> # Create consensus header
        >>> template = VCF('input.vcf.gz')
        >>> header = create_output_header(template, 'sample123', include_rescue_fields=False)
        >>> 
        >>> # Create rescue header
        >>> rescue_header = create_output_header(template, 'sample123', include_rescue_fields=True)
        >>> 
        >>> # Use with pysam
        >>> vcf_out = pysam.VariantFile('output.vcf.gz', 'wz', header=header)
    """
    import pysam
    
    # Convert cyvcf2 header to pysam header
    header_str = str(template_header.raw_header)
    new_header = pysam.VariantHeader()
    
    # Add lines from template (excluding sample line)
    for line in header_str.strip().split('\n'):
        if line.startswith('#CHROM'):
            continue
        if line.startswith('#'):
            new_header.add_line(line)
    
    # Helpers to avoid duplicate header entries when input template already defines them
    def add_info_safe(h, ident, number, typ, desc):
        if ident not in h.info:
            h.info.add(ident, number, typ, desc)
    
    def add_filter_safe(h, ident, number, typ, desc):
        if ident not in h.filters:
            h.filters.add(ident, number, typ, desc)

    # Add custom INFO fields for caller aggregation
    add_info_safe(new_header, 'N_CALLERS', '1', 'Integer', 'Total number of aggregated variant callers (excludes consensus)')
    add_info_safe(new_header, 'CALLERS', '.', 'String', 'List of all aggregated variant callers with modality prefix (pipe-separated)')
    add_info_safe(new_header, 'N_SUPPORT_CALLERS', '1', 'Integer', 'Number of variant callers that detected this variant (excludes consensus)')
    add_info_safe(new_header, 'CALLERS_SUPPORT', '.', 'String', 'Variant callers that detected this variant with modality prefix (pipe-separated, excludes consensus)')
    add_info_safe(new_header, 'N_CONSENSUS_SUPPORT', '1', 'Integer', 'Number of consensus callers that detected this variant')
    add_info_safe(new_header, 'CONSENSUS_SUPPORT', '.', 'String', 'Consensus callers that detected this variant with modality prefix (pipe-separated)')
    add_info_safe(new_header, 'N_DNA_CALLERS', '1', 'Integer', 'Total number of DNA variant callers used in analysis')
    add_info_safe(new_header, 'N_RNA_CALLERS', '1', 'Integer', 'Total number of RNA variant callers used in analysis')
    add_info_safe(new_header, 'N_DNA_CALLERS_SUPPORT', '1', 'Integer', 'Number of DNA callers that detected this variant')
    add_info_safe(new_header, 'N_RNA_CALLERS_SUPPORT', '1', 'Integer', 'Number of RNA callers that detected this variant')
    add_info_safe(new_header, 'FILTERS_ORIGINAL', '.', 'String', 'Original filter values from each caller with modality prefix (format: MODALITY_caller:filter|..., excludes consensus)')
    add_info_safe(new_header, 'FILTERS_NORMALIZED', '.', 'String', 'Normalized filter categories with modality prefix (format: MODALITY_caller:filter|..., excludes consensus)')
    add_info_safe(new_header, 'FILTERS_CATEGORY', '.', 'String', 'Filter categories with modality prefix (format: MODALITY_caller:category|..., excludes consensus)')
    add_info_safe(new_header, 'UNIFIED_FILTER', '1', 'String', 'Unified filter status based on majority of all callers')
    add_info_safe(new_header, 'UNIFIED_FILTER_DNA', '1', 'String', 'Unified filter status for DNA callers')
    add_info_safe(new_header, 'UNIFIED_FILTER_RNA', '1', 'String', 'Unified filter status for RNA callers')
    
    # Consensus flags
    add_info_safe(new_header, 'PASSES_CONSENSUS', '1', 'String', 'Whether variant passes consensus threshold overall (YES/NO)')
    add_info_safe(new_header, 'PASSES_CONSENSUS_DNA', '1', 'String', 'Whether variant passes DNA consensus threshold (YES/NO)')
    add_info_safe(new_header, 'PASSES_CONSENSUS_RNA', '1', 'String', 'Whether variant passes RNA consensus threshold (YES/NO)')
    
    # Quality aggregation
    add_info_safe(new_header, 'QUAL_MEAN', '1', 'Float', 'Mean QUAL score across callers')
    add_info_safe(new_header, 'QUAL_MIN', '1', 'Float', 'Minimum QUAL score across callers')
    add_info_safe(new_header, 'QUAL_MAX', '1', 'Float', 'Maximum QUAL score across callers')
    
    # Genotype aggregation
    add_info_safe(new_header, 'CONSENSUS_GT', '1', 'String', 'Consensus genotype across callers')
    add_info_safe(new_header, 'GT_BY_CALLER', '.', 'String', 'Genotypes from each caller with modality prefix (format: MODALITY_caller:GT|...)')
    
    # Depth aggregation
    add_info_safe(new_header, 'DP_MEAN', '1', 'Float', 'Mean depth across callers')
    add_info_safe(new_header, 'DP_MIN', '1', 'Integer', 'Minimum depth across callers')
    add_info_safe(new_header, 'DP_MAX', '1', 'Integer', 'Maximum depth across callers')
    add_info_safe(new_header, 'DP_BY_CALLER', '.', 'String', 'Depth values from each caller with modality prefix (format: MODALITY_caller:DP|...)')
    
    # VAF aggregation
    add_info_safe(new_header, 'VAF_MEAN', '1', 'Float', 'Mean variant allele frequency across callers')
    add_info_safe(new_header, 'VAF_MIN', '1', 'Float', 'Minimum VAF across callers')
    add_info_safe(new_header, 'VAF_MAX', '1', 'Float', 'Maximum VAF across callers')
    add_info_safe(new_header, 'VAF_BY_CALLER', '.', 'String', 'VAF values from each caller with modality prefix (format: MODALITY_caller:VAF|...)')

    # Rescue indicator
    add_info_safe(new_header, 'RESCUED', '1', 'String', 'Variant included via cross-modality consensus (YES/NO)')
    
    # Add rescue-specific modality tracking fields if requested
    if include_rescue_fields:
        add_info_safe(new_header, 'MODALITIES', '.', 'String', 'Modalities where variant was detected (pipe-separated)')
        add_info_safe(new_header, 'CALLERS_BY_MODALITY', '.', 'String', 'Callers grouped by modality (format: modality:caller1,caller2|...)')
        add_info_safe(new_header, 'DNA_SUPPORT', '1', 'Integer', 'Number of DNA callers supporting this variant')
        add_info_safe(new_header, 'RNA_SUPPORT', '1', 'Integer', 'Number of RNA callers supporting this variant')
        add_info_safe(new_header, 'CROSS_MODALITY', '1', 'String', 'Whether variant has cross-modality support (YES/NO)')
        
        # Modality-specific statistics
        add_info_safe(new_header, 'DP_DNA_MEAN', '1', 'Float', 'Mean depth across DNA callers')
        add_info_safe(new_header, 'DP_RNA_MEAN', '1', 'Float', 'Mean depth across RNA callers')
        add_info_safe(new_header, 'VAF_DNA_MEAN', '1', 'Float', 'Mean VAF across DNA callers')
        add_info_safe(new_header, 'VAF_RNA_MEAN', '1', 'Float', 'Mean VAF across RNA callers')
    
    # Add classification INFO fields
    add_info_safe(new_header, 'VC', '1', 'String', 'Variant Classification: Somatic, Germline, Reference, or Artifact')
    add_info_safe(new_header, 'VC_CALLERS', '.', 'String', 'Classification by each caller (format: caller:classification|...)')
    add_info_safe(new_header, 'VC_CONSENSUS', '1', 'String', 'Consensus biological classification across callers')
    
    # Add sample
    if (sample_name if sample_name else 'SAMPLE') not in new_header.samples:
        new_header.add_sample(sample_name if sample_name else 'SAMPLE')
    
    # Add biological category FILTER values (standardized classification)
    # These are the canonical FILTER values based on variant classification
    add_filter_safe(new_header, 'Somatic', None, None, 'High-confidence somatic variant specific to tumor')
    add_filter_safe(new_header, 'Germline', None, None, 'Germline variant detected in normal sample')
    add_filter_safe(new_header, 'Reference', None, None, 'Reference call - no variant detected')
    add_filter_safe(new_header, 'Artifact', None, None, 'Low quality variant, technical artifact, or inconsistent classifications across callers')
    add_filter_safe(new_header, 'NoConsensus', None, None, 'Does not meet consensus threshold or does not fit other classification categories')
    
    return new_header



def write_union_vcf(variant_data, template_header, sample_name, out_file, output_format, 
                    all_callers, modality_map=None, snv_threshold=2, indel_threshold=2):
    """
    Write union VCF with all variants and aggregated information using pysam.
    
    This function writes aggregated variant data to a VCF file, including all
    consensus/rescue INFO fields and properly formatted variant records. Variants
    are sorted by genomic position before writing.
    
    Args:
        variant_data (dict): Dictionary of variant data keyed by variant_key.
            Each value should be an aggregated variant dict from aggregate_variants().
        template_header (cyvcf2.VCF): cyvcf2 VCF header object to use as template
            for creating the output header.
        sample_name (str): Sample name for the output VCF.
        out_file (str): Output file path (e.g., 'output.vcf.gz').
        output_format (str): Output format - 'vcf', 'vcf.gz', or 'bcf'.
        all_callers (list): List of all caller names in the analysis. Used to
            populate N_CALLERS and CALLERS INFO fields. Should exclude consensus callers.
        modality_map (dict, optional): Dictionary mapping caller names to modality
            ('DNA' or 'RNA'). If provided, modality-specific INFO fields will be
            written (MODALITIES, CALLERS_BY_MODALITY, DNA_SUPPORT, RNA_SUPPORT,
            CROSS_MODALITY, DP_DNA_MEAN, DP_RNA_MEAN, VAF_DNA_MEAN, VAF_RNA_MEAN).
            Caller names will be prefixed with modality (e.g., DNA_mutect2).
            Default: None
        snv_threshold (int, optional): SNV consensus threshold for classification logic.
            Default: 2
        indel_threshold (int, optional): Indel consensus threshold for classification logic.
            Default: 2
    
    Returns:
        int: Number of variants written to the output file.
    
    Notes:
        - Variants are sorted by chromosome and position before writing
        - Chromosome sorting uses contig order from template header
        - Special chromosome handling: X=23, Y=24, M/MT=25
        - FILTER field is set based on unified filter computation
        - NoConsensus filter is added for variants not passing consensus threshold
        - Progress is printed every 10,000 variants
        - When modality_map is provided, caller names are prefixed with modality
    
    Example:
        >>> from cyvcf2 import VCF
        >>> 
        >>> # Write consensus VCF
        >>> template = VCF('input.vcf.gz')
        >>> n_written = write_union_vcf(
        ...     aggregated_data, template, 'sample123', 'consensus.vcf.gz',
        ...     'vcf.gz', ['mutect2', 'strelka', 'deepsomatic']
        ... )
        >>> print(f"Wrote {n_written} variants")
        >>> 
        >>> # Write rescue VCF with modality information
        >>> modality_map = {'mutect2': 'DNA', 'strelka': 'DNA', 'deepsomatic': 'RNA'}
        >>> n_written = write_union_vcf(
        ...     rescue_data, template, 'sample123', 'rescued.vcf.gz',
        ...     'vcf.gz', ['mutect2', 'strelka', 'deepsomatic'], modality_map
        ... )
    """
    import pysam
    from statistics import mean
    from collections import Counter
    
    print(f"- Writing union VCF to {out_file}")
    
    # Helper function to add modality prefix to caller name
    def prefix_caller(caller, modality_map):
        """Add modality prefix to caller name if modality_map is provided."""
        if not modality_map:
            return caller
        
        # Check if caller already has modality prefix (DNA_ or RNA_)
        if caller.startswith('DNA_') or caller.startswith('RNA_'):
            return caller
        
        # Add prefix if caller is in modality_map
        if caller in modality_map:
            modality = modality_map[caller]
            return f"{modality}_{caller}"
        
        return caller
    
    # Helper function to safely add modality prefix (prevents double-prefixing)
    def safe_prefix_caller(caller_name, modality):
        """Safely add modality prefix to caller name, preventing double-prefixing."""
        if caller_name.startswith(f'{modality}_'):
            return caller_name  # Already prefixed
        return f'{modality}_{caller_name}'
    
    # Helper function to check if caller is a consensus caller
    def is_consensus_caller(caller):
        """Check if caller is a consensus caller."""
        # Standardized naming: DNA_consensus, RNA_consensus
        return caller in ['DNA_consensus', 'RNA_consensus'] or (
            caller.endswith('_consensus') and 
            any(caller.startswith(prefix) for prefix in ['DNA_', 'RNA_'])
        )
    
    # Create output header with rescue fields if modality_map is provided
    include_rescue_fields = modality_map is not None
    output_header = create_output_header(template_header, sample_name, include_rescue_fields)
    
    # Determine write mode
    mode = 'w'
    if output_format == 'vcf.gz':
        mode = 'wz'
    elif output_format == 'bcf':
        mode = 'wb'
    
    # Open output VCF
    vcf_out = pysam.VariantFile(out_file, mode, header=output_header)
    
    # Get chromosome order for sorting
    chrom_order = {}
    for idx, line in enumerate(str(template_header.raw_header).split('\n')):
        if line.startswith('##contig=<ID='):
            chrom_name = line.split('ID=')[1].split(',')[0].split('>')[0]
            chrom_order[normalize_chromosome(chrom_name)] = idx
    
    # Sort variants by position
    def sort_key(item):
        vkey, data = item
        chrom = normalize_chromosome(data['CHROM'])
        chrom_idx = chrom_order.get(chrom, 999999)
        # Try to convert to int for numeric sorting, otherwise use string
        try:
            if chrom.isdigit():
                chrom_idx = int(chrom)
            elif chrom == 'X':
                chrom_idx = 23
            elif chrom == 'Y':
                chrom_idx = 24
            elif chrom == 'M' or chrom == 'MT':
                chrom_idx = 25
        except:
            pass
        return (chrom_idx, data['POS'])
    
    sorted_variants = sorted(variant_data.items(), key=sort_key)
    
    # Write each variant
    written_count = 0
    for vkey, data in sorted_variants:
        # Create new record
        alleles = tuple([data['REF']] + data['ALT'].split(','))
        
        record = vcf_out.new_record(
            contig=data['CHROM'],
            start=data['POS'] - 1,  # pysam uses 0-based
            stop=data['POS'] - 1 + len(data['REF']),
            alleles=alleles,
            id=';'.join(data['ids']) if data['ids'] else None,
            qual=mean(data['qualities']) if data['qualities'] else None
        )
        
        # Add aggregated INFO fields
        # Separate consensus callers from actual variant callers
        actual_callers_in_variant = [c for c in data['callers'] if not is_consensus_caller(c)]
        consensus_callers_in_variant = [c for c in data['callers'] if is_consensus_caller(c)]
        
        # Prefix caller names with modality if modality_map is provided
        prefixed_all_callers = [prefix_caller(c, modality_map) for c in all_callers]
        prefixed_support_callers = [prefix_caller(c, modality_map) for c in actual_callers_in_variant]
        prefixed_consensus_callers = [prefix_caller(c, modality_map) for c in consensus_callers_in_variant]
        
        record.info['N_CALLERS'] = len(all_callers)
        record.info['CALLERS'] = '|'.join(prefixed_all_callers)
        record.info['N_SUPPORT_CALLERS'] = len(set(actual_callers_in_variant))
        record.info['CALLERS_SUPPORT'] = '|'.join(prefixed_support_callers)
        
        # Add consensus support tracking
        if consensus_callers_in_variant:
            record.info['N_CONSENSUS_SUPPORT'] = len(set(consensus_callers_in_variant))
            record.info['CONSENSUS_SUPPORT'] = '|'.join(prefixed_consensus_callers)
        
        # Add DNA/RNA caller counts if modality_map is provided
        if modality_map:
            # Total number of DNA/RNA callers used in analysis (fixed)
            total_dna_callers = sum(1 for c in all_callers if modality_map.get(c) == 'DNA')
            total_rna_callers = sum(1 for c in all_callers if modality_map.get(c) == 'RNA')
            record.info['N_DNA_CALLERS'] = total_dna_callers
            record.info['N_RNA_CALLERS'] = total_rna_callers
            
            # Number of DNA/RNA callers that detected THIS variant
            dna_callers_support = sum(1 for c in actual_callers_in_variant if modality_map.get(c) == 'DNA')
            rna_callers_support = sum(1 for c in actual_callers_in_variant if modality_map.get(c) == 'RNA')
            record.info['N_DNA_CALLERS_SUPPORT'] = dna_callers_support
            record.info['N_RNA_CALLERS_SUPPORT'] = rna_callers_support
        
        # Prefix filter fields with caller names (with modality prefix) - EXCLUDE consensus
        # Note: Use proper VCF escaping for special characters instead of replacement
        prefixed_filters_original = []
        prefixed_filters_normalized = []
        prefixed_filters_category = []
        for i, caller in enumerate(data['callers']):
            if not is_consensus_caller(caller):  # Skip consensus callers
                prefixed_caller = prefix_caller(caller, modality_map)
                if i < len(data['filters_original']):
                    # Use proper VCF escaping for FILTER values (critical for Strelka compound filters like "LowDepth;LowEVS")
                    from vcf_utils.vcf_escaping import escape_filter_for_info
                    filter_val = escape_filter_for_info(data['filters_original'][i])
                    prefixed_filters_original.append(f"{prefixed_caller}:{filter_val}")
                if i < len(data['filters_normalized']):
                    from vcf_utils.vcf_escaping import escape_filter_for_info
                    filter_val = escape_filter_for_info(data['filters_normalized'][i])
                    prefixed_filters_normalized.append(f"{prefixed_caller}:{filter_val}")
                if i < len(data['filters_category']):
                    from vcf_utils.vcf_escaping import escape_filter_for_info
                    filter_val = escape_filter_for_info(data['filters_category'][i])
                    prefixed_filters_category.append(f"{prefixed_caller}:{filter_val}")
        
        record.info['FILTERS_ORIGINAL'] = '|'.join(prefixed_filters_original) if prefixed_filters_original else '.'
        record.info['FILTERS_NORMALIZED'] = '|'.join(prefixed_filters_normalized) if prefixed_filters_normalized else '.'
        record.info['FILTERS_CATEGORY'] = '|'.join(prefixed_filters_category) if prefixed_filters_category else '.'
        
        # Determine unified biological classification using new classification functions
        from vcf_utils.classification import compute_unified_classification_consensus, compute_unified_classification_rescue
        
        if modality_map:
            # Rescue mode: use new rescue classification logic
            unified_classification = compute_unified_classification_rescue(data, modality_map)
        else:
            # Consensus mode: use new consensus classification logic
            unified_classification = compute_unified_classification_consensus(data, snv_threshold, indel_threshold)
        
        record.info['UNIFIED_FILTER'] = unified_classification
        record.filter.clear()
        record.filter.add(unified_classification)
        
        # Compute modality-specific unified biological classifications if modality_map provided
        if modality_map:
            # DNA unified classification - use DNA consensus label if available
            dna_label = None
            for i, c in enumerate(data['callers']):
                if is_consensus_caller(c) and 'DNA' in c.upper():
                    dna_label = data['filters_normalized'][i]
                    break
            if dna_label:
                record.info['UNIFIED_FILTER_DNA'] = dna_label
            else:
                # Fallback: derive from individual DNA callers via majority vote
                dna_filters = [data['filters_normalized'][i] for i, c in enumerate(data['callers']) 
                               if not is_consensus_caller(c) and modality_map.get(c) == 'DNA']
                if dna_filters:
                    dna_counts = Counter(dna_filters)
                    dna_max_count = max(dna_counts.values())
                    dna_most_common = [cls for cls, count in dna_counts.items() if count == dna_max_count]
                    priority = ['Somatic', 'Germline', 'Reference', 'Artifact']
                    dna_unified = dna_most_common[0]
                    for cls in priority:
                        if cls in dna_most_common:
                            dna_unified = cls
                            break
                    record.info['UNIFIED_FILTER_DNA'] = dna_unified
            
            # RNA unified classification - use RNA consensus label if available
            rna_label = None
            for i, c in enumerate(data['callers']):
                if is_consensus_caller(c) and 'RNA' in c.upper():
                    rna_label = data['filters_normalized'][i]
                    break
            if rna_label:
                record.info['UNIFIED_FILTER_RNA'] = rna_label
            else:
                # Fallback: derive from individual RNA callers via majority vote
                rna_filters = [data['filters_normalized'][i] for i, c in enumerate(data['callers']) 
                               if not is_consensus_caller(c) and modality_map.get(c) == 'RNA']
                if rna_filters:
                    rna_counts = Counter(rna_filters)
                    rna_max_count = max(rna_counts.values())
                    rna_most_common = [cls for cls, count in rna_counts.items() if count == rna_max_count]
                    priority = ['Somatic', 'Germline', 'Reference', 'Artifact']
                    rna_unified = rna_most_common[0]
                    for cls in priority:
                        if cls in rna_most_common:
                            rna_unified = cls
                            break
                    record.info['UNIFIED_FILTER_RNA'] = rna_unified
        
        # Add consensus flag for informational purposes (but don't override FILTER)
        # The FILTER field is now set by the classification functions above
        record.info['PASSES_CONSENSUS'] = 'YES' if data['passes_consensus'] else 'NO'
        
        # Note: We no longer override FILTER here because the classification functions
        # already handle NoConsensus assignment based on more sophisticated logic
        
        # Add modality-specific consensus flags if modality_map provided
        if modality_map:
            # Check if DNA consensus or RNA consensus is in the callers (fixed naming)
            has_dna_consensus = any(c == 'DNA_consensus' for c in data['callers'])
            has_rna_consensus = any(c == 'RNA_consensus' for c in data['callers'])
            
            # Check if there are any DNA or RNA callers in this variant
            has_dna_callers = any(not is_consensus_caller(c) and modality_map.get(c) == 'DNA' for c in data['callers'])
            has_rna_callers = any(not is_consensus_caller(c) and modality_map.get(c) == 'RNA' for c in data['callers'])
            
            if has_dna_consensus or has_dna_callers:
                record.info['PASSES_CONSENSUS_DNA'] = 'YES' if has_dna_consensus else 'NO'
            if has_rna_consensus or has_rna_callers:
                record.info['PASSES_CONSENSUS_RNA'] = 'YES' if has_rna_consensus else 'NO'

        # Rescue indicator - use actual rescued status from variant data
        record.info['RESCUED'] = 'YES' if data.get('rescued', False) else 'NO'
        
        # Add modality-specific fields if modality_map is provided
        if modality_map:
            # Group callers by modality - EXCLUDE consensus callers
            modalities = set()
            dna_callers = []
            rna_callers = []
            
            for caller in data['callers']:
                if not is_consensus_caller(caller):  # Skip consensus callers
                    modality = modality_map.get(caller, 'UNKNOWN')
                    modalities.add(modality)
                    if modality == 'DNA':
                        dna_callers.append(caller)
                    elif modality == 'RNA':
                        rna_callers.append(caller)
            
            # Write modality fields
            record.info['MODALITIES'] = '|'.join(sorted(modalities)) if modalities else '.'
            
            # Format: DNA:caller1,caller2|RNA:caller3,caller4 - EXCLUDE consensus
            callers_by_mod = []
            if dna_callers:
                callers_by_mod.append(f"DNA:{','.join(dna_callers)}")
            if rna_callers:
                callers_by_mod.append(f"RNA:{','.join(rna_callers)}")
            if callers_by_mod:
                record.info['CALLERS_BY_MODALITY'] = '|'.join(callers_by_mod)
            
            record.info['DNA_SUPPORT'] = len(set(dna_callers))
            record.info['RNA_SUPPORT'] = len(set(rna_callers))
            record.info['CROSS_MODALITY'] = 'YES' if len(modalities) > 1 else 'NO'
            
            # Calculate modality-specific statistics
            agg = data['gt_aggregated']
            
            # DNA statistics
            dna_dp_values = []
            dna_vaf_values = []
            for i, caller in enumerate(data['callers']):
                if modality_map.get(caller) == 'DNA':
                    if i < len(agg['dp_by_caller']) and agg['dp_by_caller'][i] is not None:
                        dna_dp_values.append(agg['dp_by_caller'][i])
                    if i < len(agg['vaf_by_caller']) and agg['vaf_by_caller'][i] is not None:
                        dna_vaf_values.append(agg['vaf_by_caller'][i])
            
            if dna_dp_values:
                record.info['DP_DNA_MEAN'] = round(mean(dna_dp_values), 1)
            if dna_vaf_values:
                record.info['VAF_DNA_MEAN'] = round(mean(dna_vaf_values), 4)
            
            # RNA statistics
            rna_dp_values = []
            rna_vaf_values = []
            for i, caller in enumerate(data['callers']):
                if modality_map.get(caller) == 'RNA':
                    if i < len(agg['dp_by_caller']) and agg['dp_by_caller'][i] is not None:
                        rna_dp_values.append(agg['dp_by_caller'][i])
                    if i < len(agg['vaf_by_caller']) and agg['vaf_by_caller'][i] is not None:
                        rna_vaf_values.append(agg['vaf_by_caller'][i])
            
            if rna_dp_values:
                record.info['DP_RNA_MEAN'] = round(mean(rna_dp_values), 1)
            if rna_vaf_values:
                record.info['VAF_RNA_MEAN'] = round(mean(rna_vaf_values), 4)
        
        # Add quality statistics
        if data['qualities']:
            record.info['QUAL_MEAN'] = round(mean(data['qualities']), 2)
            record.info['QUAL_MIN'] = round(min(data['qualities']), 2)
            record.info['QUAL_MAX'] = round(max(data['qualities']), 2)
        
        # Add genotype aggregation
        agg = data['gt_aggregated']
        
        if agg['consensus_gt']:
            record.info['CONSENSUS_GT'] = agg['consensus_gt']
        
        # Prefix GT_BY_CALLER with modality - EXCLUDE consensus
        if agg['gt_by_caller']:
            prefixed_gt_by_caller = []
            for i, caller in enumerate(data['callers']):
                if not is_consensus_caller(caller):  # Skip consensus callers
                    prefixed_caller = prefix_caller(caller, modality_map)
                    if i < len(agg['gt_by_caller']):
                        prefixed_gt_by_caller.append(f"{prefixed_caller}:{agg['gt_by_caller'][i]}")
            if prefixed_gt_by_caller:
                record.info['GT_BY_CALLER'] = '|'.join(prefixed_gt_by_caller)
        
        # Add depth statistics with modality prefix - EXCLUDE consensus
        if agg['dp_values']:
            record.info['DP_MEAN'] = round(agg['dp_mean'], 1)
            record.info['DP_MIN'] = agg['dp_min']
            record.info['DP_MAX'] = agg['dp_max']
            
            prefixed_dp_by_caller = []
            for i, caller in enumerate(data['callers']):
                if not is_consensus_caller(caller):  # Skip consensus callers
                    prefixed_caller = prefix_caller(caller, modality_map)
                    if i < len(agg['dp_by_caller']):
                        dp_val = '' if agg['dp_by_caller'][i] is None else str(agg['dp_by_caller'][i])
                        prefixed_dp_by_caller.append(f"{prefixed_caller}:{dp_val}")
            if prefixed_dp_by_caller:
                record.info['DP_BY_CALLER'] = '|'.join(prefixed_dp_by_caller)
        
        # Add VAF statistics with modality prefix - EXCLUDE consensus
        if agg['vaf_values']:
            record.info['VAF_MEAN'] = round(agg['vaf_mean'], 4)
            record.info['VAF_MIN'] = round(agg['vaf_min'], 4)
            record.info['VAF_MAX'] = round(agg['vaf_max'], 4)
            
            prefixed_vaf_by_caller = []
            for i, caller in enumerate(data['callers']):
                if not is_consensus_caller(caller):  # Skip consensus callers
                    prefixed_caller = prefix_caller(caller, modality_map)
                    if i < len(agg['vaf_by_caller']):
                        vaf_val = '' if agg['vaf_by_caller'][i] is None else f"{agg['vaf_by_caller'][i]:.4f}"
                        prefixed_vaf_by_caller.append(f"{prefixed_caller}:{vaf_val}")
            if prefixed_vaf_by_caller:
                record.info['VAF_BY_CALLER'] = '|'.join(prefixed_vaf_by_caller)
        
        # Write record
        vcf_out.write(record)
        written_count += 1
        
        if written_count % 10000 == 0:
            print(f"  - Written {written_count:,} variants...")
    
    vcf_out.close()
    print(f"- Successfully wrote {written_count:,} variants to {out_file}")
    
    return written_count
