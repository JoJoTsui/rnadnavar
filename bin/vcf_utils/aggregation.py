"""
Variant aggregation functions for consensus and rescue workflows.

This module provides reusable functions for aggregating variants
from multiple callers and across modalities. It handles reading VCF files,
extracting genotype information, and combining variant data with support
for both within-modality consensus and cross-modality rescue operations.

Functions:
    extract_genotype_info: Extract genotype, depth, and VAF from a variant
    aggregate_genotypes: Aggregate genotype information across callers
    read_variants_from_vcf: Read variants from a single VCF file
    aggregate_variants: Aggregate variants from multiple collections

Example:
    Reading and aggregating variants from multiple callers:
    
    >>> from vcf_utils.aggregation import read_variants_from_vcf, aggregate_variants
    >>> 
    >>> # Read variants from each caller
    >>> mutect2_vars = read_variants_from_vcf('mutect2.vcf.gz', 'mutect2')
    >>> strelka_vars = read_variants_from_vcf('strelka.vcf.gz', 'strelka')
    >>> deepsomatic_vars = read_variants_from_vcf('deepsomatic.vcf.gz', 'deepsomatic',
    ...                                           exclude_refcall=True, exclude_germline=True)
    >>> 
    >>> # Prepare collections for aggregation
    >>> collections = [
    ...     ('mutect2', mutect2_vars, None),
    ...     ('strelka', strelka_vars, None),
    ...     ('deepsomatic', deepsomatic_vars, None)
    ... ]
    >>> 
    >>> # Aggregate with consensus thresholds
    >>> aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)
    >>> 
    >>> # Access aggregated data
    >>> for vkey, data in aggregated.items():
    ...     print(f"Variant: {vkey}")
    ...     print(f"  Supported by: {data['callers']}")
    ...     print(f"  Passes consensus: {data['passes_consensus']}")
    ...     print(f"  Mean VAF: {data['gt_aggregated']['vaf_mean']}")
"""


def extract_genotype_info(variant, caller):
    """
    Extract genotype, depth, and VAF information from variant.
    
    This function extracts genotype information from a cyvcf2 Variant object,
    with special handling for caller-specific formats (particularly Strelka).
    It attempts to extract GT, DP, AD, VAF, and GQ fields using multiple
    strategies to handle different VCF format variations.
    
    Args:
        variant (cyvcf2.Variant): Variant object from cyvcf2 VCF reader
        caller (str): Caller name for caller-specific parsing. Special handling
            is provided for 'strelka' which uses non-standard format fields
            (TAR, TIR, AU, CU, GU, TU, SGT)
    
    Returns:
        dict: Genotype information dictionary with keys:
            - GT (str or None): Genotype string (e.g., '0/1', '1|1', '0|1')
                where '/' indicates unphased and '|' indicates phased
            - DP (int or None): Total read depth at the variant position
            - AD (str or None): Allele depths as comma-separated string (e.g., '50,25')
            - VAF (float or None): Variant allele frequency (0.0 to 1.0)
            - GQ (int or None): Genotype quality score
    
    Notes:
        - For Strelka variants, the function uses TAR/TIR (indels) or AU/CU/GU/TU
          (SNVs) to compute AD and VAF when standard fields are missing
        - VAF is calculated from AD if not directly available
        - DP is calculated from AD if not directly available
        - Returns None for fields that cannot be extracted
        - Prints warnings to stdout if extraction errors occur
    
    Example:
        >>> from cyvcf2 import VCF
        >>> vcf = VCF('variants.vcf.gz')
        >>> for variant in vcf:
        ...     gt_info = extract_genotype_info(variant, 'mutect2')
        ...     print(f"GT: {gt_info['GT']}, DP: {gt_info['DP']}, VAF: {gt_info['VAF']}")
        ...     break
        GT: 0/1, DP: 150, VAF: 0.25
    """
    info = {
        'GT': None,
        'DP': None,
        'AD': None,
        'VAF': None,
        'GQ': None,
    }
    
    try:
        def fmt(key):
            try:
                return variant.format(key)
            except Exception:
                return None

        # Extract genotype
        try:
            if variant.genotypes and len(variant.genotypes) > 0:
                gt = variant.genotypes[0]
                if len(gt) >= 2:
                    a1 = '.' if gt[0] == -1 else str(gt[0])
                    a2 = '.' if gt[1] == -1 else str(gt[1])
                    phased = len(gt) > 2 and bool(gt[2])
                    sep = '|' if phased else '/'
                    info['GT'] = f"{a1}{sep}{a2}"
        except Exception:
            pass

        # Depth
        dp = fmt('DP')
        if dp is not None and len(dp) > 0:
            try:
                info['DP'] = int(dp[0]) if dp[0] is not None else None
            except Exception:
                pass

        # Allele depth
        ad = fmt('AD')
        if ad is not None and len(ad) > 0 and ad[0] is not None:
            try:
                info['AD'] = ','.join(map(str, ad[0]))
            except Exception:
                pass

        # VAF/AF
        if info['VAF'] is None:
            for vaf_field in ['VAF', 'AF', 'FREQ', 'FA']:
                vaf = fmt(vaf_field)
                if vaf is not None and len(vaf) > 0 and vaf[0] is not None:
                    val = vaf[0]
                    if isinstance(val, str) and '%' in val:
                        val = float(val.replace('%', '')) / 100.0
                    try:
                        info['VAF'] = float(val)
                    except Exception:
                        pass
                    break

        # Strelka-specific parsing when AD/AF are missing
        if caller.lower() == 'strelka':
            def choose_row():
                try:
                    tar = fmt('TAR')
                    tir = fmt('TIR')
                    if tar is not None and tir is not None and len(tar) > 0 and len(tir) > 0:
                        scores = []
                        for i in range(max(len(tar), len(tir))):
                            r = tar[i] if i < len(tar) else None
                            a = tir[i] if i < len(tir) else None
                            r1 = r[0] if hasattr(r, '__len__') and len(r) > 0 else r
                            a1 = a[0] if hasattr(a, '__len__') and len(a) > 0 else a
                            s = 0
                            try:
                                s += int(r1) if r1 is not None else 0
                            except Exception:
                                pass
                            try:
                                s += int(a1) if a1 is not None else 0
                            except Exception:
                                pass
                            scores.append(s)
                        if scores:
                            return max(range(len(scores)), key=lambda i: scores[i])
                    au, cu, gu, tu = fmt('AU'), fmt('CU'), fmt('GU'), fmt('TU')
                    bases = [au, cu, gu, tu]
                    if any(b is not None for b in bases):
                        nrows = 0
                        for b in bases:
                            if b is not None:
                                nrows = max(nrows, len(b))
                        scores = []
                        for i in range(nrows):
                            s = 0
                            for b in bases:
                                if b is None or i >= len(b):
                                    continue
                                v = b[i]
                                v1 = v[0] if hasattr(v, '__len__') and len(v) > 0 else v
                                try:
                                    s += int(v1) if v1 is not None else 0
                                except Exception:
                                    pass
                            scores.append(s)
                        if scores:
                            return max(range(len(scores)), key=lambda i: scores[i])
                except Exception:
                    pass
                return 0

            row_idx = choose_row()
            if info['AD'] is None:
                tir = fmt('TIR')
                tar = fmt('TAR')
                if tir is not None and tar is not None and len(tir) > 0 and len(tar) > 0:
                    try:
                        t_alt = tir[row_idx]
                        t_ref = tar[row_idx]
                        a1 = t_alt[0] if hasattr(t_alt, '__len__') and len(t_alt) > 0 else t_alt
                        r1 = t_ref[0] if hasattr(t_ref, '__len__') and len(t_ref) > 0 else t_ref
                        alt_c = int(a1) if a1 is not None else None
                        ref_c = int(r1) if r1 is not None else None
                        if ref_c is not None and alt_c is not None:
                            info['AD'] = f"{ref_c},{alt_c}"
                            if info['DP'] is None:
                                info['DP'] = ref_c + alt_c
                            if info['VAF'] is None and (ref_c + alt_c) > 0:
                                info['VAF'] = alt_c / (ref_c + alt_c)
                    except Exception:
                        pass
                else:
                    au = fmt('AU')
                    cu = fmt('CU')
                    gu = fmt('GU')
                    tu = fmt('TU')
                    if au is not None and cu is not None and gu is not None and tu is not None:
                        base_map = {'A': au, 'C': cu, 'G': gu, 'T': tu}
                        ref_base = variant.REF if hasattr(variant, 'REF') else None
                        alt_base = None
                        try:
                            if hasattr(variant, 'ALT') and variant.ALT and len(variant.ALT) > 0:
                                alt_base = variant.ALT[0]
                        except Exception:
                            pass
                        try:
                            ref_pair = base_map.get(str(ref_base))
                            alt_pair = base_map.get(str(alt_base))
                            if ref_pair is not None and alt_pair is not None and len(ref_pair) > 0 and len(alt_pair) > 0:
                                ref_counts = ref_pair[row_idx]
                                alt_counts = alt_pair[row_idx]
                                r1 = ref_counts[0] if hasattr(ref_counts, '__len__') and len(ref_counts) > 0 else ref_counts
                                a1 = alt_counts[0] if hasattr(alt_counts, '__len__') and len(alt_counts) > 0 else alt_counts
                                ref_c = int(r1) if r1 is not None else None
                                alt_c = int(a1) if a1 is not None else None
                                info['AD'] = f"{ref_c},{alt_c}"
                                if info['DP'] is None:
                                    info['DP'] = ref_c + alt_c
                                if info['VAF'] is None and (ref_c + alt_c) > 0:
                                    info['VAF'] = alt_c / (ref_c + alt_c)
                        except Exception:
                            pass

            # Use SGT when GT is absent
            if info['GT'] is None:
                sgt = fmt('SGT')
                if sgt is not None and len(sgt) > 0 and sgt[0] is not None:
                    try:
                        info['GT'] = str(sgt[0])
                    except Exception:
                        pass

        # Calculate VAF from AD if still not available
        if info['VAF'] is None and info['AD'] is not None:
            try:
                ad_values = [int(x) for x in info['AD'].split(',')]
                if len(ad_values) >= 2 and sum(ad_values) > 0:
                    info['VAF'] = ad_values[1] / sum(ad_values)
            except Exception:
                pass

        if info['DP'] is None and info['AD'] is not None:
            try:
                ad_values = [int(x) for x in info['AD'].split(',')]
                info['DP'] = sum(ad_values)
            except Exception:
                pass

        # Genotype quality
        gq = fmt('GQ')
        if gq is not None and len(gq) > 0:
            try:
                info['GQ'] = int(gq[0]) if gq[0] is not None else None
            except Exception:
                pass

    except Exception as e:
        print(f"Warning: Error extracting genotype info from {caller}: {e}")
    
    return info



def aggregate_genotypes(genotypes_by_caller, callers_order):
    """
    Aggregate genotype information across callers.
    
    This function combines genotype information from multiple variant callers,
    computing consensus genotypes and summary statistics for depth and VAF.
    The consensus genotype is determined by majority vote across callers.
    
    Args:
        genotypes_by_caller (dict): Dictionary mapping caller names to genotype
            info dictionaries (as returned by extract_genotype_info). Each
            genotype info dict should contain GT, DP, VAF, AD, and GQ fields.
        callers_order (list): Ordered list of caller names. This order is used
            to maintain consistent ordering in the output lists (gt_by_caller,
            dp_by_caller, vaf_by_caller).
    
    Returns:
        dict: Aggregated genotype statistics with keys:
            - consensus_gt (str or None): Most common genotype across callers
            - gt_list (list): List of all non-None genotypes
            - dp_values (list): List of all non-None depth values
            - dp_mean (float or None): Mean depth across callers
            - dp_min (int or None): Minimum depth across callers
            - dp_max (int or None): Maximum depth across callers
            - vaf_values (list): List of all non-None VAF values
            - vaf_mean (float or None): Mean VAF across callers
            - vaf_min (float or None): Minimum VAF across callers
            - vaf_max (float or None): Maximum VAF across callers
            - ad_values (list): List of all non-None allele depth strings
            - gq_values (list): List of all non-None genotype quality values
            - gt_by_caller (list): Ordered list of genotypes (or '.' if None) by caller
            - dp_by_caller (list): Ordered list of depths (or None) by caller
            - vaf_by_caller (list): Ordered list of VAFs (or None) by caller
    
    Example:
        >>> genotypes = {
        ...     'mutect2': {'GT': '0/1', 'DP': 100, 'VAF': 0.25, 'AD': '75,25', 'GQ': 99},
        ...     'strelka': {'GT': '0/1', 'DP': 120, 'VAF': 0.27, 'AD': '88,32', 'GQ': 95},
        ...     'deepsomatic': {'GT': '0/1', 'DP': 110, 'VAF': 0.26, 'AD': '81,29', 'GQ': 98}
        ... }
        >>> callers = ['mutect2', 'strelka', 'deepsomatic']
        >>> agg = aggregate_genotypes(genotypes, callers)
        >>> print(f"Consensus GT: {agg['consensus_gt']}")
        Consensus GT: 0/1
        >>> print(f"Mean DP: {agg['dp_mean']:.1f}, Mean VAF: {agg['vaf_mean']:.3f}")
        Mean DP: 110.0, Mean VAF: 0.260
    """
    from collections import defaultdict
    from statistics import mean
    
    agg = {
        'consensus_gt': None,
        'gt_list': [],
        'dp_values': [],
        'dp_mean': None,
        'dp_min': None,
        'dp_max': None,
        'vaf_values': [],
        'vaf_mean': None,
        'vaf_min': None,
        'vaf_max': None,
        'ad_values': [],
        'gq_values': [],
        'gt_by_caller': [],
        'dp_by_caller': [],
        'vaf_by_caller': [],
    }
    
    gt_counts = defaultdict(int)
    
    for caller in callers_order:
        info = genotypes_by_caller.get(caller, {})
        # Collect GTs
        if info and 'GT' in info and info['GT'] is not None:
            agg['gt_list'].append(info['GT'])
            gt_counts[info['GT']] += 1
            agg['gt_by_caller'].append(info['GT'])
        else:
            agg['gt_by_caller'].append('.')
        
        # Collect DPs
        if info and 'DP' in info and info['DP'] is not None:
            agg['dp_values'].append(info['DP'])
            agg['dp_by_caller'].append(info['DP'])
        else:
            agg['dp_by_caller'].append(None)
        
        # Collect VAFs
        if info and 'VAF' in info and info['VAF'] is not None:
            agg['vaf_values'].append(info['VAF'])
            agg['vaf_by_caller'].append(info['VAF'])
        else:
            agg['vaf_by_caller'].append(None)
        
        # Collect ADs
        if info and 'AD' in info and info['AD'] is not None:
            agg['ad_values'].append(info['AD'])
        
        # Collect GQs
        if info and 'GQ' in info and info['GQ'] is not None:
            agg['gq_values'].append(info['GQ'])
    
    # Determine consensus genotype (most common)
    if gt_counts:
        agg['consensus_gt'] = max(gt_counts, key=gt_counts.get)
    
    # Calculate DP statistics
    if agg['dp_values']:
        agg['dp_mean'] = mean(agg['dp_values'])
        agg['dp_min'] = min(agg['dp_values'])
        agg['dp_max'] = max(agg['dp_values'])
    
    # Calculate VAF statistics
    if agg['vaf_values']:
        agg['vaf_mean'] = mean(agg['vaf_values'])
        agg['vaf_min'] = min(agg['vaf_values'])
        agg['vaf_max'] = max(agg['vaf_values'])
    
    return agg



def read_variants_from_vcf(vcf_path, caller_name, modality=None, 
                          exclude_refcall=False, exclude_germline=False,
                          classify_variants=True):
    """
    Read variants from a single VCF file with biological classification.
    
    This function reads all variants from a VCF file and extracts relevant
    information including position, alleles, filters, quality, and genotype
    data. It supports optional filtering to exclude reference calls and
    germline variants, and can tag variants with modality information.
    Additionally, it classifies variants into biological categories.
    
    Args:
        vcf_path (str): Path to VCF file (can be .vcf, .vcf.gz, or .bcf)
        caller_name (str): Name of the variant caller (e.g., 'mutect2', 'strelka')
        modality (str, optional): Modality tag to add to variants ('DNA' or 'RNA').
            If None, no modality tag is added. Default: None
        exclude_refcall (bool, optional): If True, exclude variants with 'RefCall'
            in the FILTER field. Useful for DeepSomatic output. Default: False
        exclude_germline (bool, optional): If True, exclude variants with 'GERMLINE'
            in the FILTER field. Default: False
        classify_variants (bool, optional): If True, classify variants into
            biological categories (Somatic, Germline, Reference, Artifact).
            Default: True
    
    Returns:
        dict: Dictionary mapping variant_key to variant_data. Each variant_data
            dict contains:
            - CHROM (str): Chromosome name
            - POS (int): 1-based position
            - REF (str): Reference allele
            - ALT (str): Comma-separated alternate allele(s)
            - is_snv (bool): True if variant is a single nucleotide variant
            - caller (str): Caller name
            - modality (str, optional): Modality tag (only if modality arg provided)
            - filter_original (str): Original FILTER field value
            - filter_normalized (str): Normalized filter string (unified)
            - filter_category (str): Filter category
            - classification (str): Biological classification (if classify_variants=True)
            - quality (float or None): QUAL score
            - genotype (dict): Genotype information from extract_genotype_info()
            - id (str or None): Variant ID from ID field
    
    Notes:
        - Variant keys are generated using variant_key() function which creates
          unique identifiers in format: "chrom:pos:ref:alt"
        - Filter normalization and categorization use functions from vcf_utils.filters
        - Genotype extraction uses extract_genotype_info() with caller-specific handling
        - Classification uses vcf_utils.classification module for unified categorization
        - Filter values are unified: PASS, GERMLINE, RefCall, LowQuality
    
    Example:
        >>> # Read DNA variants from Mutect2 with classification
        >>> dna_vars = read_variants_from_vcf('mutect2.vcf.gz', 'mutect2', 
        ...                                    modality='DNA', exclude_germline=True)
        >>> print(f"Read {len(dna_vars)} DNA variants from Mutect2")
        >>> 
        >>> # Read RNA variants from Strelka
        >>> rna_vars = read_variants_from_vcf('strelka.vcf.gz', 'strelka', modality='RNA')
        >>> 
        >>> # Access variant data with classification
        >>> for vkey, data in list(dna_vars.items())[:3]:
        ...     print(f"{vkey}: {data['classification']} -> {data['filter_normalized']}")
    """
    from cyvcf2 import VCF
    from vcf_utils.io_utils import variant_key, is_snv
    from vcf_utils.filters import normalize_filter, categorize_filter
    
    # Import classification functions if needed
    if classify_variants:
        from vcf_utils.classification import (
            classify_variant_from_record,
            get_sample_indices,
            normalize_filter_value
        )
    
    variants = {}
    
    vcf = VCF(vcf_path)
    
    # Get sample indices for classification (needed for Strelka)
    sample_indices = None
    if classify_variants:
        sample_indices = get_sample_indices(vcf, caller_name)
    
    for variant in vcf:
        # Get filter and check exclusions
        filter_str = variant.FILTER if variant.FILTER else 'PASS'
        
        if exclude_refcall and 'RefCall' in filter_str:
            continue
        if exclude_germline and 'GERMLINE' in filter_str:
            continue
        
        vkey = variant_key(variant, use_cyvcf2=True)
        
        # Classify variant before creating data dict
        classification = None
        if classify_variants:
            try:
                classification = classify_variant_from_record(variant, caller_name, sample_indices)
            except Exception as e:
                classification = "Artifact"
        
        # Create variant data
        data = {
            'CHROM': variant.CHROM,
            'POS': variant.POS,
            'REF': variant.REF,
            'ALT': ','.join(variant.ALT) if variant.ALT else '.',
            'is_snv': is_snv(variant.REF, variant.ALT),
            'caller': caller_name,
            'filter_original': filter_str,
            'filter_normalized': normalize_filter_value(classification) if classification else normalize_filter(filter_str),
            'quality': float(variant.QUAL) if variant.QUAL is not None else None,
            'genotype': extract_genotype_info(variant, caller_name),
            'id': variant.ID if variant.ID else None,
        }
        
        # Add classification
        if classification:
            data['classification'] = classification
        
        # Add filter category based on normalized (unified) filter
        data['filter_category'] = categorize_filter(data['filter_normalized'])
        
        # Add modality if provided
        if modality:
            data['modality'] = modality
        
        variants[vkey] = data
    
    return variants



def aggregate_variants(variant_collections, snv_threshold=2, indel_threshold=2):
    """
    Aggregate variants from multiple collections.
    
    This function combines variants from multiple callers and/or modalities,
    merging information for variants at the same genomic position. It applies
    consensus thresholds to determine which variants have sufficient support,
    and aggregates genotype information across all supporting callers.
    
    Args:
        variant_collections (list): List of tuples, each containing:
            - caller_name (str): Name of the caller
            - variant_dict (dict): Dictionary of variants from read_variants_from_vcf()
            - modality (str or None): Modality tag ('DNA', 'RNA', or None)
        snv_threshold (int, optional): Minimum number of callers required for
            a SNV to pass consensus. Default: 2
        indel_threshold (int, optional): Minimum number of callers required for
            an indel to pass consensus. Default: 2
    
    Returns:
        dict: Dictionary mapping variant_key to aggregated variant_data. Each
            aggregated variant_data dict contains:
            - CHROM (str): Chromosome name
            - POS (int): 1-based position
            - REF (str): Reference allele
            - ALT (str): Comma-separated alternate allele(s)
            - is_snv (bool): True if variant is SNV
            - callers (list): List of all caller names that detected this variant
            - modalities (list): List of modalities (if provided in collections)
            - caller_modality_map (dict): Maps caller name to modality
            - filters_original (list): Original filter strings from each caller
            - filters_normalized (list): Normalized filter strings
            - filters_category (list): Filter categories
            - qualities (list): Quality scores from each caller
            - genotypes (dict): Maps caller name to genotype info dict
            - ids (list): Variant IDs from each caller
            - support_callers (set): Set of unique callers supporting this variant
            - passes_consensus (bool): True if variant meets consensus threshold
            - gt_aggregated (dict): Aggregated genotype statistics from
                aggregate_genotypes()
    
    Notes:
        - Variants are matched by genomic position (chrom:pos:ref:alt)
        - The consensus threshold is applied based on variant type (SNV vs indel)
        - Genotype aggregation is performed automatically for all variants
        - Modality information is preserved when provided in variant_collections
    
    Example:
        >>> # Aggregate within-modality consensus (DNA only)
        >>> mutect2_vars = read_variants_from_vcf('mutect2.vcf.gz', 'mutect2')
        >>> strelka_vars = read_variants_from_vcf('strelka.vcf.gz', 'strelka')
        >>> collections = [
        ...     ('mutect2', mutect2_vars, None),
        ...     ('strelka', strelka_vars, None)
        ... ]
        >>> aggregated = aggregate_variants(collections, snv_threshold=2, indel_threshold=2)
        >>> 
        >>> # Count consensus variants
        >>> consensus_count = sum(1 for v in aggregated.values() if v['passes_consensus'])
        >>> print(f"Variants passing consensus: {consensus_count}")
        >>> 
        >>> # Aggregate cross-modality rescue
        >>> dna_vars = read_variants_from_vcf('dna_consensus.vcf.gz', 'consensus', 'DNA')
        >>> rna_vars = read_variants_from_vcf('rna_consensus.vcf.gz', 'consensus', 'RNA')
        >>> collections = [
        ...     ('consensus', dna_vars, 'DNA'),
        ...     ('consensus', rna_vars, 'RNA')
        ... ]
        >>> rescue_aggregated = aggregate_variants(collections, snv_threshold=1, indel_threshold=1)
    """
    from collections import defaultdict
    
    aggregated = defaultdict(lambda: {
        'CHROM': None,
        'POS': None,
        'REF': None,
        'ALT': None,
        'is_snv': None,
        'callers': [],
        'modalities': [],
        'caller_modality_map': {},
        'filters_original': [],
        'filters_normalized': [],
        'filters_category': [],
        'qualities': [],
        'genotypes': {},
        'ids': [],
        'support_callers': set(),
    })
    
    # Process each variant collection
    for caller_name, variant_dict, modality in variant_collections:
        for vkey, variant_data in variant_dict.items():
            data = aggregated[vkey]
            
            # Store variant location info (first time)
            if data['CHROM'] is None:
                data['CHROM'] = variant_data['CHROM']
                data['POS'] = variant_data['POS']
                data['REF'] = variant_data['REF']
                data['ALT'] = variant_data['ALT']
                data['is_snv'] = variant_data['is_snv']
            
            # Store caller-specific information
            data['callers'].append(caller_name)
            data['support_callers'].add(caller_name)
            
            # Store modality information
            if modality:
                data['modalities'].append(modality)
                data['caller_modality_map'][caller_name] = modality
            
            # Store filters
            data['filters_original'].append(variant_data['filter_original'])
            data['filters_normalized'].append(variant_data['filter_normalized'])
            data['filters_category'].append(variant_data['filter_category'])
            
            # Store quality
            if variant_data['quality'] is not None:
                data['qualities'].append(variant_data['quality'])
            
            # Store ID
            if variant_data['id']:
                data['ids'].append(variant_data['id'])
            
            # Store genotype information
            data['genotypes'][caller_name] = variant_data['genotype']
    
    # Apply consensus thresholds and aggregate genotypes
    for vkey, data in aggregated.items():
        n_support = len(data['support_callers'])
        
        # Determine if variant passes consensus threshold
        if data['is_snv']:
            data['passes_consensus'] = n_support >= snv_threshold
        else:
            data['passes_consensus'] = n_support >= indel_threshold
        
        # Aggregate genotype information
        data['gt_aggregated'] = aggregate_genotypes(data['genotypes'], data['callers'])
    
    return dict(aggregated)
