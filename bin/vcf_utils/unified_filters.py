"""
Unified filtering logic for VCF files (consensus and rescue).

This module provides common filtering functions that work with both
consensus VCFs and rescue VCFs, applying consistent RaVeX filtering rules.
"""

import pysam


def is_multiallelic(ref, alts):
    """
    Check if variant is multiallelic (multiple ALT alleles).
    
    Args:
        ref: Reference allele
        alts: List or tuple of alternate alleles
    
    Returns:
        bool: True if multiple ALT alleles present
    """
    if not alts:
        return False
    # Filter out None and empty values
    valid_alts = [a for a in alts if a and a != '.']
    return len(valid_alts) > 1


def filter_homopolymer(ref_context, alt, hp_length=6):
    """Check if the variant is in a homopolymer context"""
    if len(alt) > 1:
        alt = alt[1:]
    elif alt == "-":
        return False
    elif ref_context is None:
        return None
    
    try:
        length_to_consider = int((len(ref_context) - 1) / 2)
    except TypeError:
        return None
    
    ref_context = list(ref_context)
    ref_context[length_to_consider] = alt
    ref_context = "".join(ref_context)
    context_to_consider = ref_context[length_to_consider - hp_length + 1 : length_to_consider + hp_length]
    
    for idx, base in enumerate(context_to_consider):
        context_window = context_to_consider[idx : idx + hp_length]
        if len(context_window) < hp_length:
            break
        elif len(set(context_window)) == 1:
            return True
    return False


def add_context(chrom, pos, ref, genome, flank=10):
    """Extract genomic context around variant"""
    if pos < 10:
        flank = pos - 1
    try:
        context = genome.fetch(chrom, pos - 1 - flank, pos + flank).upper()
    except (ValueError, KeyError):
        return None
    
    if ref != "-":
        try:
            assert ref[0] == context[flank]
        except (AssertionError, IndexError):
            return None
    return context


def is_noncoding(csq_string):
    """Check if variant is in noncoding region based on VEP consequence"""
    noncoding_list = [
        "intron_variant", "intergenic_variant", "non_coding_transcript_variant",
        "non_coding_transcript_exon_variant", "mature_miRNA_variant",
        "regulatory_region_variant", "IGR", "INTRON", "RNA"
    ]
    if not csq_string:
        return False
    
    # Parse first consequence from CSQ field
    consequences = csq_string.split("&")[0].split(",")[0]
    return consequences in noncoding_list


def is_ig_pseudo(biotype):
    """Check if gene is IG or pseudogene"""
    if not biotype:
        return False
    ig_pseudo_patterns = [
        "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
        "TR_C_gene", "TR_J_gene", "TR_V_gene", "pseudogene"
    ]
    return any(pattern in biotype for pattern in ig_pseudo_patterns)


# Common gnomAD AF INFO field names, matched case-insensitively.
# 'GNOMAD_AF' (uppercase) is what gnomad_annotator.py writes via
# bcftools annotate -c INFO/GNOMAD_AF:=INFO/AF; the rest are historical aliases.
GNOMAD_AF_FIELDS = ['MAX_AF', 'gnomAD_AF', 'AF_gnomad', 'gnomad_AF', 'GNOMAD_AF']


def get_gnomad_af(variant, use_cyvcf2=False):
    """
    Extract maximum gnomAD allele frequency from INFO field.

    Field names are matched case-insensitively against GNOMAD_AF_FIELDS.
    Single-element tuple/list values (pysam Number=A fields) are unwrapped.

    Args:
        variant: Variant object (cyvcf2 or pysam)
        use_cyvcf2: If True, variant is cyvcf2.Variant, else pysam.VariantRecord

    Returns:
        float: Maximum gnomAD AF, or 0.0 if not found
    """
    info = variant.INFO if use_cyvcf2 else variant.info
    try:
        # cyvcf2 INFO iterates as (key, value) pairs; pysam info iterates keys
        names = [key for key, _ in info] if use_cyvcf2 else list(info)
    except TypeError:
        return 0.0
    names_by_lower = {name.lower(): name for name in names}

    for field in GNOMAD_AF_FIELDS:
        actual_name = names_by_lower.get(field.lower())
        if actual_name is None:
            continue
        try:
            value = info.get(actual_name)
        except (ValueError, TypeError, KeyError):
            continue
        if isinstance(value, (tuple, list)):
            value = value[0] if value else None
        if value is not None:
            try:
                return float(value)
            except (ValueError, TypeError):
                pass
    return 0.0


# Biological-class FILTER values written by the consensus/rescue stages.
# These are biological labels (the pipeline's label contract), not caller
# rejection filters, so they must not trip the vc_filter check.
BIOLOGICAL_CLASS_FILTERS = frozenset(
    ["Somatic", "Germline", "Reference", "Artifact", "NoConsensus", "RNAedit"]
)


def get_tumor_alt_count(variant, use_cyvcf2=False):
    """
    Best available tumor alt-read count for a variant.

    Consensus/rescue VCFs written by write_union_vcf have no FORMAT column;
    they carry per-caller tumor alt counts in INFO (ALT_COUNT_BY_CALLER and
    ALT_COUNT_MAX), extracted from the tumor sample of each caller VCF. The
    min_alt_reads rule uses the MAXIMUM tumor alt count across callers: a
    record has alt support if at least one caller saw enough tumor alt reads
    (mirrors the consensus semantics where one strong caller suffices).
    FORMAT/AD of the first sample is only a fallback for VCFs that genuinely
    have it (e.g. raw caller VCFs).

    Args:
        variant: Variant object (cyvcf2 or pysam)
        use_cyvcf2: If True, variant is cyvcf2.Variant, else pysam.VariantRecord

    Returns:
        int: Tumor alt-read count, or 0 if no evidence is available
    """
    info = variant.INFO if use_cyvcf2 else variant.info

    def info_get(key):
        try:
            return info.get(key)
        except (KeyError, ValueError, TypeError):
            return None

    # 1. Consensus aggregate INFO: ALT_COUNT_MAX
    alt_max = info_get("ALT_COUNT_MAX")
    if isinstance(alt_max, (tuple, list)):
        alt_max = alt_max[0] if alt_max else None
    if alt_max is not None:
        try:
            return int(alt_max)
        except (ValueError, TypeError):
            pass

    # 2. Consensus per-caller INFO: ALT_COUNT_BY_CALLER ("caller:count|...")
    by_caller = info_get("ALT_COUNT_BY_CALLER")
    if isinstance(by_caller, (tuple, list)):
        by_caller = by_caller[0] if by_caller else None
    if by_caller:
        counts = []
        for entry in str(by_caller).split("|"):
            if ":" not in entry:
                continue
            _, value = entry.rsplit(":", 1)
            try:
                counts.append(int(value))
            except (ValueError, TypeError):
                continue
        if counts:
            return max(counts)

    # 3. Fallback: FORMAT/AD of the first sample (raw caller VCFs)
    try:
        if use_cyvcf2:
            ad = variant.format("AD")
            if ad is not None and len(ad) > 0 and len(ad[0]) > 1:
                return int(ad[0][1])
        else:
            if len(variant.samples) > 0:
                sample = variant.samples[0]
                if "AD" in sample:
                    ad = sample["AD"]
                    if ad and len(ad) > 1:
                        return int(ad[1])
    except (KeyError, IndexError, TypeError):
        pass

    return 0


def get_csq_field(variant, field_name, use_cyvcf2=False):
    """
    Extract field from VEP CSQ annotation.
    
    Args:
        variant: Variant object (cyvcf2 or pysam)
        field_name: Name of CSQ field to extract
        use_cyvcf2: If True, variant is cyvcf2.Variant, else pysam.VariantRecord
    
    Returns:
        str or None: Field value or None if not found
    """
    try:
        if use_cyvcf2:
            csq = variant.INFO.get('CSQ')
        else:
            # pysam raises ValueError("Invalid header") if INFO/CSQ is undefined
            csq = variant.info.get('CSQ')
    except (KeyError, AttributeError, ValueError):
        return None
    
    if not csq:
        return None
    
    # Parse first annotation
    annotations = csq.split(',')[0].split('|')
    
    # Common VEP field positions (may need adjustment based on VEP version)
    csq_fields = {
        'Consequence': 1,
        'BIOTYPE': 7,
        'SYMBOL': 3
    }
    
    if field_name in csq_fields:
        try:
            return annotations[csq_fields[field_name]]
        except IndexError:
            return None
    return None


def apply_ravex_filters(variant, args, genome=None, blacklist_regions=None, use_cyvcf2=False):
    """
    Apply RaVeX filtering rules to a variant.
    
    Args:
        variant: Variant object (cyvcf2 or pysam)
        args: Argparse namespace with filter parameters
        genome: pysam.FastaFile for reference sequence (optional)
        blacklist_regions: List of (chrom, start, end) tuples (optional)
        use_cyvcf2: If True, variant is cyvcf2.Variant, else pysam.VariantRecord
    
    Returns:
        list: List of filter flags that failed
    """
    filters = []
    
    # Get variant info
    if use_cyvcf2:
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        alt = variant.ALT[0] if variant.ALT else ""
    else:
        chrom = variant.chrom
        pos = variant.pos
        ref = variant.ref
        alt = variant.alts[0] if variant.alts else ""
    
    # Multiallelic filter
    if hasattr(args, 'filter_multiallelic') and args.filter_multiallelic:
        if use_cyvcf2:
            if is_multiallelic(ref, variant.ALT):
                filters.append("multiallelic")
        else:
            if is_multiallelic(ref, variant.alts):
                filters.append("multiallelic")
    
    # Alt read count filter. Alt support comes from the consensus per-caller
    # tumor alt-count INFO when present (consensus VCFs have no FORMAT
    # column); FORMAT/AD is only a fallback for raw caller VCFs.
    if hasattr(args, 'min_alt_reads'):
        alt_count = get_tumor_alt_count(variant, use_cyvcf2=use_cyvcf2)
        if alt_count < args.min_alt_reads:
            filters.append("min_alt_reads")
    
    # gnomAD filter
    if hasattr(args, 'gnomad_thr'):
        gnomad_af = get_gnomad_af(variant, use_cyvcf2=use_cyvcf2)
        if gnomad_af >= args.gnomad_thr:
            filters.append("gnomad")
    
    # Blacklist filter
    if blacklist_regions:
        for bl_chrom, bl_start, bl_end in blacklist_regions:
            if chrom == bl_chrom and bl_start <= pos <= bl_end:
                filters.append("blacklist")
                break
    
    # Noncoding filter
    csq = get_csq_field(variant, 'Consequence', use_cyvcf2=use_cyvcf2)
    if is_noncoding(csq):
        filters.append("noncoding")
    
    # IG/pseudogene filter
    biotype = get_csq_field(variant, 'BIOTYPE', use_cyvcf2=use_cyvcf2)
    if is_ig_pseudo(biotype):
        filters.append("ig_pseudo")
    
    # Homopolymer filter
    if genome:
        context = add_context(chrom, pos, ref, genome)
        if context and filter_homopolymer(context, alt):
            filters.append("homopolymer")
    
    return filters


def add_filter_info_fields(header):
    """
    Add RaVeX filter INFO field definitions to VCF header.
    
    Args:
        header: pysam.VariantHeader object
    
    Returns:
        pysam.VariantHeader: Updated header
    """
    # Add RaVeX_FILTER INFO field if not present
    if 'RaVeX_FILTER' not in header.info:
        header.info.add(
            'RaVeX_FILTER', '.', 'String',
            'RaVeX filter reasons: semicolon-separated list of filter flags'
        )
    
    # Add individual filter flag INFO fields
    filter_flags = {
        'min_alt_reads': 'Variant filtered due to insufficient alternate reads',
        'gnomad': 'Variant filtered due to high gnomAD allele frequency',
        'blacklist': 'Variant in blacklisted region',
        'noncoding': 'Variant in noncoding region',
        'ig_pseudo': 'Variant in immunoglobulin or pseudogene',
        'homopolymer': 'Variant in homopolymer region',
        'vc_filter': 'Variant failed variant caller filters',
        'not_consensus': 'Variant not in consensus',
        'multiallelic': 'Multiallelic site with multiple ALT alleles'
    }
    
    for flag, description in filter_flags.items():
        if flag not in header.info:
            header.info.add(flag, '0', 'Flag', description)
    
    return header


def add_classification_filters_to_header(header):
    """
    Add biological classification FILTER definitions to VCF header.
    
    Args:
        header: pysam.VariantHeader object
    
    Returns:
        pysam.VariantHeader: Updated header
    """
    classification_filters = {
        'Somatic': 'Somatic variant',
        'Germline': 'Germline variant',
        'Reference': 'Reference/wildtype',
        'Artifact': 'Artifact/technical error'
    }
    
    for filt, description in classification_filters.items():
        if filt not in header.filters:
            header.filters.add(filt, None, None, description)
    
    return header
