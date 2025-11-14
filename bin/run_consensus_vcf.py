#!/usr/bin/env python3
"""
VCF Union/Consensus Script
Union mode: keeps ALL variants from all callers with aggregated information
Uses cyvcf2 for fast reading and pysam for writing
"""
import argparse
import sys
from pathlib import Path
from collections import defaultdict, Counter
from statistics import mean, median
import re
from cyvcf2 import VCF
import pysam


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


def argparser():
    parser = argparse.ArgumentParser(
        description="Union all variants from multiple VCF files with caller aggregation"
    )
    parser.add_argument("--input_dir", required=True, help="Directory containing input VCF files")
    parser.add_argument("--out_prefix", required=True, help="Prefix for output files")
    parser.add_argument("--snv_thr", type=int, default=2, help="Number of callers required for SNV consensus (default: 2)")
    parser.add_argument("--indel_thr", type=int, default=2, help="Number of callers required for indel consensus (default: 2)")
    parser.add_argument("--output_format", choices=['vcf', 'vcf.gz', 'bcf'], default='vcf.gz', 
                       help="Output format (default: vcf.gz)")
    parser.add_argument("--sample_name", help="Sample name (if not specified, use first sample from VCF)")
    parser.add_argument("--exclude_refcall", action='store_true', 
                       help="Exclude variants marked as RefCall (DeepSomatic reference calls)")
    parser.add_argument("--exclude_germline", action='store_true', 
                       help="Exclude variants marked as GERMLINE")
    return parser.parse_args()


def get_caller_name(filename):
    """Extract caller name from filename"""
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
    """Normalize chromosome names for consistent sorting"""
    chrom = str(chrom).replace('chr', '').replace('Chr', '')
    return chrom


def variant_key(variant, use_cyvcf2=True):
    """
    Create unique key for variant including all ALT alleles
    Works with both cyvcf2.Variant and dict representations
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
    """Determine if variant is SNV (all alleles are single base)"""
    if len(ref) != 1:
        return False
    for alt in alt_list:
        if len(alt) != 1 or alt == '*':  # * indicates deletion in some formats
            return False
    return True


def normalize_filter(filter_str):
    """
    Normalize filter string using mapping
    Handles: PASS, ., missing, and compound filters (semicolon-separated)
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
    Categorize normalized filter into major categories
    Returns: primary category or 'Other'
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


def extract_genotype_info(variant, caller):
    """Extract genotype, depth, and VAF information from variant"""
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


def should_exclude_variant(filter_str, exclude_refcall, exclude_germline):
    """Determine if variant should be excluded based on filters"""
    if exclude_refcall and 'RefCall' in filter_str:
        return True
    if exclude_germline and 'GERMLINE' in filter_str:
        return True
    return False


def read_all_variants(vcf_files, snv_threshold, indel_threshold, exclude_refcall, exclude_germline):
    """
    Read all variants from all VCF files using cyvcf2 (single pass)
    Returns: (variant_data dict, template_header, sample_name)
    """
    print("- Reading and aggregating all variants from all VCF files")
    
    variant_data = defaultdict(lambda: {
        'CHROM': None,
        'POS': None,
        'REF': None,
        'ALT': None,
        'is_snv': None,
        'callers': [],
        'filters_original': [],
        'filters_normalized': [],
        'filters_category': [],
        'qualities': [],
        'genotypes': {},  # caller -> GT info dict
        'ids': [],
    })
    
    template_header = None
    sample_name = None
    all_callers = list(vcf_files.keys())
    
    # Statistics
    stats = {
        'total_read': 0,
        'excluded_refcall': 0,
        'excluded_germline': 0,
    }
    
    # Read each VCF file with cyvcf2
    for caller, vcf_path in vcf_files.items():
        print(f"  - Reading {caller}: {vcf_path}")
        vcf = VCF(vcf_path)
        
        # Get template header and sample name from first VCF
        if template_header is None:
            template_header = vcf
            if vcf.samples:
                sample_name = vcf.samples[0]
        
        variant_count = 0
        excluded_count = 0
        
        for variant in vcf:
            stats['total_read'] += 1
            
            # Get filter and check exclusions
            filter_str = variant.FILTER if variant.FILTER else 'PASS'
            
            if should_exclude_variant(filter_str, exclude_refcall, exclude_germline):
                excluded_count += 1
                if 'RefCall' in filter_str:
                    stats['excluded_refcall'] += 1
                if 'GERMLINE' in filter_str:
                    stats['excluded_germline'] += 1
                continue
            
            vkey = variant_key(variant, use_cyvcf2=True)
            data = variant_data[vkey]
            
            # Store variant location info (first time)
            if data['CHROM'] is None:
                data['CHROM'] = variant.CHROM
                data['POS'] = variant.POS
                data['REF'] = variant.REF
                data['ALT'] = ','.join(variant.ALT) if variant.ALT else '.'
                data['is_snv'] = is_snv(variant.REF, variant.ALT)
            
            # Store caller-specific information
            data['callers'].append(caller)
            if 'support_callers' not in data:
                data['support_callers'] = set()
            data['support_callers'].add(caller)
            
            # Store and normalize filters
            data['filters_original'].append(filter_str)
            normalized = normalize_filter(filter_str)
            data['filters_normalized'].append(normalized)
            data['filters_category'].append(categorize_filter(normalized))
            
            # Get quality score
            if variant.QUAL is not None:
                data['qualities'].append(float(variant.QUAL))
            
            # Get ID
            if variant.ID:
                data['ids'].append(variant.ID)
            
            # Extract genotype information
            gt_info = extract_genotype_info(variant, caller)
            data['genotypes'][caller] = gt_info
            
            variant_count += 1
        
        print(f"    - Processed {variant_count:,} variants from {caller} (excluded {excluded_count:,})")
    
    print(f"\n- Total variants read: {stats['total_read']:,}")
    if exclude_refcall:
        print(f"- Excluded RefCall variants: {stats['excluded_refcall']:,}")
    if exclude_germline:
        print(f"- Excluded GERMLINE variants: {stats['excluded_germline']:,}")
    print(f"- Total unique variants across all callers: {len(variant_data):,}")
    
    return variant_data, template_header, sample_name, all_callers


def aggregate_variant_info(variant_data, snv_threshold, indel_threshold):
    """
    Process variant data and compute aggregate statistics
    """
    print("- Computing aggregate statistics for all variants")
    
    stats = {
        'total_variants': len(variant_data),
        'snvs': 0,
        'indels': 0,
        'snvs_consensus': 0,
        'indels_consensus': 0,
        'single_caller': 0,
        'multi_caller': 0,
    }
    
    for vkey, data in variant_data.items():
        n_callers = len(data['callers'])
        
        # Count variant types
        if data['is_snv']:
            stats['snvs'] += 1
            if len(set(data['callers'])) >= snv_threshold:
                stats['snvs_consensus'] += 1
                data['passes_consensus'] = True
            else:
                data['passes_consensus'] = False
        else:
            stats['indels'] += 1
            if len(set(data['callers'])) >= indel_threshold:
                stats['indels_consensus'] += 1
                data['passes_consensus'] = True
            else:
                data['passes_consensus'] = False
        
        # Count caller distribution
        if len(set(data['callers'])) == 1:
            stats['single_caller'] += 1
        else:
            stats['multi_caller'] += 1
        
        # Aggregate genotype information (stable caller order)
        data['gt_aggregated'] = aggregate_genotypes(data['genotypes'], data['callers'])
    
    print(f"\n- Statistics:")
    print(f"  - Total variants: {stats['total_variants']:,}")
    print(f"  - SNVs: {stats['snvs']:,} (consensus: {stats['snvs_consensus']:,})")
    print(f"  - Indels: {stats['indels']:,} (consensus: {stats['indels_consensus']:,})")
    print(f"  - Single caller: {stats['single_caller']:,}")
    print(f"  - Multiple callers: {stats['multi_caller']:,}")
    
    return variant_data


def aggregate_genotypes(genotypes_by_caller, callers_order):
    """
    Aggregate genotype information across callers
    Returns dict with aggregated stats
    """
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


def create_output_header(template_header, sample_name):
    """Create output VCF header with all INFO fields"""
    
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
    add_info_safe(new_header, 'N_CALLERS', '1', 'Integer', 'Total number of aggregated callers')
    add_info_safe(new_header, 'CALLERS', '.', 'String', 'List of all aggregated callers (pipe-separated)')
    add_info_safe(new_header, 'N_SUPPORT_CALLERS', '1', 'Integer', 'Number of callers that detected this variant')
    add_info_safe(new_header, 'CALLERS_SUPPORT', '.', 'String', 'Callers that detected this variant (pipe-separated)')
    add_info_safe(new_header, 'FILTERS_ORIGINAL', '.', 'String', 'Original filter values from each caller (pipe-separated)')
    add_info_safe(new_header, 'FILTERS_NORMALIZED', '.', 'String', 'Normalized filter categories (pipe-separated)')
    add_info_safe(new_header, 'FILTERS_CATEGORY', '.', 'String', 'Filter categories (pipe-separated)')
    add_info_safe(new_header, 'UNIFIED_FILTER', '1', 'String', 'Unified filter status based on majority')
    
    # Consensus flags
    add_info_safe(new_header, 'PASSES_CONSENSUS', '1', 'String', 'Whether variant passes consensus threshold (YES/NO)')
    
    # Quality aggregation
    add_info_safe(new_header, 'QUAL_MEAN', '1', 'Float', 'Mean QUAL score across callers')
    add_info_safe(new_header, 'QUAL_MIN', '1', 'Float', 'Minimum QUAL score across callers')
    add_info_safe(new_header, 'QUAL_MAX', '1', 'Float', 'Maximum QUAL score across callers')
    
    # Genotype aggregation
    add_info_safe(new_header, 'CONSENSUS_GT', '1', 'String', 'Consensus genotype across callers')
    add_info_safe(new_header, 'GT_BY_CALLER', '.', 'String', 'Genotypes from each caller (pipe-separated)')
    
    # Depth aggregation
    add_info_safe(new_header, 'DP_MEAN', '1', 'Float', 'Mean depth across callers')
    add_info_safe(new_header, 'DP_MIN', '1', 'Integer', 'Minimum depth across callers')
    add_info_safe(new_header, 'DP_MAX', '1', 'Integer', 'Maximum depth across callers')
    add_info_safe(new_header, 'DP_BY_CALLER', '.', 'String', 'Depth values from each caller (pipe-separated)')
    
    # VAF aggregation
    add_info_safe(new_header, 'VAF_MEAN', '1', 'Float', 'Mean variant allele frequency across callers')
    add_info_safe(new_header, 'VAF_MIN', '1', 'Float', 'Minimum VAF across callers')
    add_info_safe(new_header, 'VAF_MAX', '1', 'Float', 'Maximum VAF across callers')
    add_info_safe(new_header, 'VAF_BY_CALLER', '.', 'String', 'VAF values from each caller (pipe-separated)')

    # Rescue indicator
    add_info_safe(new_header, 'RESCUED', '1', 'String', 'Variant included via cross-modality consensus (YES/NO)')
    
    # Add sample
    if (sample_name if sample_name else 'SAMPLE') not in new_header.samples:
        new_header.add_sample(sample_name if sample_name else 'SAMPLE')
    
    # Add unified filter categories
    add_filter_safe(new_header, 'LowQuality', None, None, 'Low quality or evidence score')
    add_filter_safe(new_header, 'LowDepth', None, None, 'Low sequencing depth')
    add_filter_safe(new_header, 'StrandBias', None, None, 'Strand bias detected')
    add_filter_safe(new_header, 'Germline', None, None, 'Likely germline variant')
    add_filter_safe(new_header, 'Artifact', None, None, 'Likely artifact or present in normal')
    add_filter_safe(new_header, 'NoConsensus', None, None, 'Does not meet consensus threshold')
    add_filter_safe(new_header, 'LowEvidenceScore', None, None, 'Low empirical variant score (Strelka EVS)')
    add_filter_safe(new_header, 'ReferenceCall', None, None, 'Called as reference (DeepSomatic RefCall)')
    
    return new_header


def write_union_vcf(variant_data, template_header, sample_name, out_file, output_format, all_callers):
    """
    Write union VCF with all variants and aggregated information using pysam
    """
    print(f"- Writing union VCF to {out_file}")
    
    # Create output header
    output_header = create_output_header(template_header, sample_name)
    
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
        record.info['N_CALLERS'] = len(all_callers)
        record.info['CALLERS'] = '|'.join(all_callers)
        record.info['N_SUPPORT_CALLERS'] = len(set(data['callers']))
        record.info['CALLERS_SUPPORT'] = '|'.join(data['callers'])
        record.info['FILTERS_ORIGINAL'] = '|'.join(data['filters_original'])
        record.info['FILTERS_NORMALIZED'] = '|'.join(data['filters_normalized'])
        record.info['FILTERS_CATEGORY'] = '|'.join(data['filters_category'])
        
        # Determine unified filter (majority vote on categories)
        pass_count = sum(1 for f in data['filters_normalized'] if f == 'PASS')
        if pass_count >= len(set(data['callers'])) / 2:
            record.info['UNIFIED_FILTER'] = 'PASS'
        else:
            record.info['UNIFIED_FILTER'] = 'FAIL'
            # Find most common category
            non_pass_cats = [c for c in data['filters_category'] if c != 'PASS']
            if non_pass_cats:
                most_common_cat = Counter(non_pass_cats).most_common(1)[0][0]
                # Map category to filter
                if 'quality' in most_common_cat:
                    record.filter.add('LowQuality')
                elif 'depth' in most_common_cat:
                    record.filter.add('LowDepth')
                elif 'bias' in most_common_cat:
                    record.filter.add('StrandBias')
                elif 'germline' in most_common_cat:
                    record.filter.add('Germline')
                elif 'artifact' in most_common_cat:
                    record.filter.add('Artifact')
        
        # Add consensus flag
        record.info['PASSES_CONSENSUS'] = 'YES' if data['passes_consensus'] else 'NO'
        if not data['passes_consensus']:
            record.filter.add('NoConsensus')

        # Rescue indicator
        record.info['RESCUED'] = 'YES' if 'consensus' in set(data['callers']) else 'NO'
        
        # Add quality statistics
        if data['qualities']:
            record.info['QUAL_MEAN'] = round(mean(data['qualities']), 2)
            record.info['QUAL_MIN'] = round(min(data['qualities']), 2)
            record.info['QUAL_MAX'] = round(max(data['qualities']), 2)
        
        # Add genotype aggregation
        agg = data['gt_aggregated']
        
        if agg['consensus_gt']:
            record.info['CONSENSUS_GT'] = agg['consensus_gt']
        
        if agg['gt_by_caller']:
            record.info['GT_BY_CALLER'] = '|'.join(agg['gt_by_caller'])
        
        # Add depth statistics
        if agg['dp_values']:
            record.info['DP_MEAN'] = round(agg['dp_mean'], 1)
            record.info['DP_MIN'] = agg['dp_min']
            record.info['DP_MAX'] = agg['dp_max']
            record.info['DP_BY_CALLER'] = '|'.join('' if v is None else str(v) for v in agg['dp_by_caller'])
        
        # Add VAF statistics
        if agg['vaf_values']:
            record.info['VAF_MEAN'] = round(agg['vaf_mean'], 4)
            record.info['VAF_MIN'] = round(agg['vaf_min'], 4)
            record.info['VAF_MAX'] = round(agg['vaf_max'], 4)
            record.info['VAF_BY_CALLER'] = '|'.join('' if v is None else f"{v:.4f}" for v in agg['vaf_by_caller'])
        
        # Write record
        vcf_out.write(record)
        written_count += 1
        
        if written_count % 10000 == 0:
            print(f"  - Written {written_count:,} variants...")
    
    vcf_out.close()
    print(f"- Successfully wrote {written_count:,} variants to {out_file}")


def main():
    args = argparser()
    
    # Find VCF files
    input_dir = Path(args.input_dir)
    vcf_files = {}
    
    print(f"Searching for VCF files in: {input_dir}")
    for vcf_path in sorted(input_dir.glob("*.vcf*")):
        if vcf_path.suffix in ['.vcf', '.gz'] or vcf_path.name.endswith('.vcf.gz'):
            caller = get_caller_name(str(vcf_path))
            vcf_files[caller] = str(vcf_path)
            print(f"  - Found: {caller} -> {vcf_path.name}")
    
    if not vcf_files:
        print(f"ERROR: No VCF files found in {input_dir}")
        sys.exit(1)
    
    print(f"\nProcessing {len(vcf_files)} VCF files with callers: {', '.join(vcf_files.keys())}")
    print(f"SNV consensus threshold: {args.snv_thr}")
    print(f"Indel consensus threshold: {args.indel_thr}")
    if args.exclude_refcall:
        print("Excluding RefCall variants")
    if args.exclude_germline:
        print("Excluding GERMLINE variants")
    print()
    
    # Read all variants from all VCFs (single pass with cyvcf2)
    variant_data, template_header, sample_name, all_callers = read_all_variants(
        vcf_files, args.snv_thr, args.indel_thr, args.exclude_refcall, args.exclude_germline
    )
    
    # Override sample name if provided
    if args.sample_name:
        sample_name = args.sample_name
    
    # Aggregate statistics
    variant_data = aggregate_variant_info(variant_data, args.snv_thr, args.indel_thr)
    
    # Write output VCF (with pysam)
    out_file = f"{args.out_prefix}.{args.output_format}"
    write_union_vcf(variant_data, template_header, sample_name, out_file, args.output_format, all_callers)
    
    print(f"\n{'='*60}")
    print("DONE! Union VCF created successfully.")
    print(f"Output: {out_file}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
