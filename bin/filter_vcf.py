#!/usr/bin/env python3
"""
Author: Adapted from filter_mutations.py by @RaqManzano
Script: Filter variants from a VCF file producing another VCF file with the new filters added.
Applies the same filtering logic as filter_mutations.py but works directly with VCF format.
"""
import argparse
import pysam


def argparser():
    parser = argparse.ArgumentParser(description="Filter VCF files with RaVeX filtering logic")
    parser.add_argument("-i", "--input", help="VCF file input", required=True)
    parser.add_argument("-o", "--output", help="VCF file output", default="filtered.vcf.gz")
    parser.add_argument("-g", "--gnomad_thr", help="Gnomad threshold for variants", default=0.0001, type=float)
    parser.add_argument("--whitelist", help="BED file with variants to keep (CHROM POS REF ALT)")
    parser.add_argument("--blacklist", help="BED file with regions to remove (CHROM START END)")
    parser.add_argument("--filters", help="Other filters to be considered as PASS", default=["PASS"], nargs="+")
    parser.add_argument("--ref", help="FASTA reference file to extract context")
    parser.add_argument("--min_alt_reads", help="Minimum alt reads", default=2, type=int)
    return parser.parse_args()


def filter_homopolymer(ref_context, alt, hp_length=6):
    """Check if the variant is in a homopolymer context"""
    homopolymer = False
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
            homopolymer = True
            break
    return homopolymer


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


def get_gnomad_af(info_dict):
    """Extract maximum gnomAD allele frequency from INFO field"""
    # Try common gnomAD field names
    for field in ['MAX_AF', 'gnomAD_AF', 'AF_gnomad', 'gnomad_AF']:
        if field in info_dict:
            try:
                return float(info_dict[field])
            except (ValueError, TypeError):
                pass
    return 0.0


def get_csq_field(info_dict, field_name):
    """Extract field from VEP CSQ annotation"""
    if 'CSQ' not in info_dict:
        return None
    
    csq = info_dict['CSQ']
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


def apply_filters(vcf_in, vcf_out, args, genome):
    """Apply RaVeX filtering logic to VCF"""
    
    # Read whitelist if provided
    whitelist_vars = set()
    if args.whitelist:
        with open(args.whitelist) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    var_id = f"{parts[0]}:{parts[1]}:{parts[2]}:{parts[3]}"
                    whitelist_vars.add(var_id)
    
    # Read blacklist if provided
    blacklist_regions = []
    if args.blacklist:
        with open(args.blacklist) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    blacklist_regions.append((parts[0], int(parts[1]), int(parts[2])))
    
    # Add RaVeX_FILTER to header
    vcf_out.header.filters.add("RaVeX_FILTER", None, None, "RaVeX filtering applied")
    vcf_out.header.info.add("RaVeX_FILTER", ".", "String", "RaVeX filter reasons")
    
    for record in vcf_in:
        filters = []
        
        # Get variant info
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        alt = str(record.alts[0]) if record.alts else ""
        var_id = f"{chrom}:{pos}:{ref}:{alt}"
        
        # Check whitelist first
        if whitelist_vars and var_id in whitelist_vars:
            record.info['RaVeX_FILTER'] = "PASS"
            vcf_out.write(record)
            continue
        
        # Get alt read count from FORMAT field
        alt_count = 0
        if record.samples:
            sample = list(record.samples.values())[0]
            if 'AD' in sample:
                ad = sample['AD']
                if len(ad) > 1:
                    alt_count = ad[1]
        
        # Apply filters
        if alt_count < args.min_alt_reads:
            filters.append("min_alt_reads")
        
        # gnomAD filter
        gnomad_af = get_gnomad_af(record.info)
        if gnomad_af >= args.gnomad_thr:
            filters.append("gnomad")
        
        # Blacklist filter
        for bl_chrom, bl_start, bl_end in blacklist_regions:
            if chrom == bl_chrom and bl_start <= pos <= bl_end:
                filters.append("blacklist")
                break
        
        # Noncoding filter
        csq = get_csq_field(record.info, 'Consequence')
        if is_noncoding(csq):
            filters.append("noncoding")
        
        # IG/pseudogene filter
        biotype = get_csq_field(record.info, 'BIOTYPE')
        if is_ig_pseudo(biotype):
            filters.append("ig_pseudo")
        
        # Homopolymer filter
        if genome:
            context = add_context(chrom, pos, ref, genome)
            if context and filter_homopolymer(context, alt):
                filters.append("homopolymer")
        
        # Variant caller filter
        if record.filter.keys():
            filter_vals = list(record.filter.keys())
            if not any(f in args.filters for f in filter_vals):
                filters.append("vc_filter")
        
        # Set filter
        if filters:
            record.info['RaVeX_FILTER'] = ";".join(filters)
            record.filter.add("RaVeX_FILTER")
        else:
            record.info['RaVeX_FILTER'] = "PASS"
        
        vcf_out.write(record)


def main():
    args = argparser()
    
    # Open genome if provided
    genome = None
    if args.ref:
        genome = pysam.FastaFile(args.ref)
    
    # Open input VCF
    vcf_in = pysam.VariantFile(args.input)
    
    # Create output VCF
    vcf_out = pysam.VariantFile(args.output, 'w', header=vcf_in.header)
    
    # Apply filters
    apply_filters(vcf_in, vcf_out, args, genome)
    
    # Close files
    vcf_in.close()
    vcf_out.close()
    
    if genome:
        genome.close()
    
    # Index output
    pysam.tabix_index(args.output, preset='vcf', force=True)
    
    print(f"Done! Filtered VCF written to '{args.output}'")


if __name__ == "__main__":
    main()
