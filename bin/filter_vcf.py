#!/usr/bin/env python3
"""
VCF Filtering Script
Applies the same filtering logic as filter_mutations.py but works natively with VCF files
Uses cyvcf2 for reading and pysam for writing
"""
import argparse
from cyvcf2 import VCF
import pysam
import sys
from pathlib import Path

# Add vcf_utils to path
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir))

from vcf_utils.filters import normalize_filter
from vcf_utils.io_utils import normalize_chromosome, variant_key
from vcf_utils.unified_filters import (
    apply_ravex_filters,
    add_filter_info_fields,
    add_classification_filters_to_header,
    is_multiallelic
)
from vcf_utils.stripped_writer import write_vcf_stripped


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
    parser.add_argument("--filter_multiallelic", help="Filter out multiallelic sites (multiple ALT alleles)", action="store_true")
    parser.add_argument("--output_stripped", help="Output stripped VCF file (multiallelic-filtered, FORMAT-stripped). If not provided, will be <output>.stripped.vcf.gz")
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


def get_gnomad_af_cyvcf2(variant):
    """Extract maximum gnomAD allele frequency from INFO field"""
    # Try common gnomAD field names
    for field in ['MAX_AF', 'gnomAD_AF', 'AF_gnomad', 'gnomad_AF']:
        try:
            value = variant.INFO.get(field)
            if value is not None:
                return float(value)
        except (ValueError, TypeError):
            pass
    return 0.0


def get_csq_field_cyvcf2(variant, field_name):
    """Extract field from VEP CSQ annotation"""
    csq = variant.INFO.get('CSQ')
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


def apply_filters(vcf_in, args, genome):
    """Apply RaVeX filtering logic to VCF using cyvcf2 for reading"""
    
    # Read whitelist if provided
    whitelist_vars = set()
    if args.whitelist:
        with open(args.whitelist) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    # Use normalized chromosome for consistency
                    chrom = normalize_chromosome(parts[0])
                    var_id = f"{chrom}:{parts[1]}:{parts[2]}:{parts[3]}"
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
                    # Normalize chromosome for consistency
                    chrom = normalize_chromosome(parts[0])
                    blacklist_regions.append((chrom, int(parts[1]), int(parts[2])))
    
    seen_variants = set()
    filtered_variants = []
    
    for variant in vcf_in:
        filters = []
        
        # Use shared variant_key function for consistency
        vkey = variant_key(variant, use_cyvcf2=True)
        
        # Skip duplicates
        if vkey in seen_variants:
            continue
        seen_variants.add(vkey)
        
        # Get variant info
        chrom = normalize_chromosome(variant.CHROM)
        pos = variant.POS
        ref = variant.REF
        alt = variant.ALT[0] if variant.ALT else ""
        
        # Check whitelist first
        if whitelist_vars and vkey in whitelist_vars:
            filtered_variants.append((variant, "PASS", []))
            continue
        
        # Get alt read count from FORMAT field
        alt_count = 0
        try:
            ad = variant.format('AD')
            if ad is not None and len(ad) > 0 and len(ad[0]) > 1:
                alt_count = ad[0][1]
        except (KeyError, IndexError, TypeError):
            pass
        
        # Apply filters
        if alt_count < args.min_alt_reads:
            filters.append("min_alt_reads")
        
        # gnomAD filter
        gnomad_af = get_gnomad_af_cyvcf2(variant)
        if gnomad_af >= args.gnomad_thr:
            filters.append("gnomad")
        
        # Blacklist filter
        for bl_chrom, bl_start, bl_end in blacklist_regions:
            if chrom == bl_chrom and bl_start <= pos <= bl_end:
                filters.append("blacklist")
                break
        
        # Noncoding filter
        csq = get_csq_field_cyvcf2(variant, 'Consequence')
        if is_noncoding(csq):
            filters.append("noncoding")
        
        # IG/pseudogene filter
        biotype = get_csq_field_cyvcf2(variant, 'BIOTYPE')
        if is_ig_pseudo(biotype):
            filters.append("ig_pseudo")
        
        # Homopolymer filter
        if genome:
            context = add_context(chrom, pos, ref, genome)
            if context and filter_homopolymer(context, alt):
                filters.append("homopolymer")
        
        # Variant caller filter - use normalize_filter for consistency
        if variant.FILTER:
            normalized = normalize_filter(variant.FILTER)
            if normalized != 'PASS' and variant.FILTER not in args.filters:
                filters.append("vc_filter")
        
            # Multiallelic filter (multiple ALT alleles)
            if args.filter_multiallelic and is_multiallelic(ref, variant.ALT):
                filters.append("multiallelic")
        
        # Store variant with filter info
        filter_status = "PASS" if not filters else "RaVeX_FILTER"
        filtered_variants.append((variant, filter_status, filters))
    
    return filtered_variants


def write_filtered_vcf(filtered_variants, input_vcf_path, output_path, strip_format=False):
    """Write filtered VCF using pysam for proper handling"""
    
    # Open input VCF with pysam to get header
    input_vcf = pysam.VariantFile(input_vcf_path)
    
    # Create new header
    new_header = input_vcf.header.copy()
    
    # Strip FORMAT fields if requested
    if strip_format:
        # Remove all FORMAT definitions. Keep sample names in header; we simply
        # won't write any FORMAT/sample data in records, which results in no
        # FORMAT column in output.
        for fmt_key in list(new_header.formats.keys()):
            new_header.formats.remove_header(fmt_key)
    
    # Ensure biological classification FILTER fields are defined
    classification_filters = {
        'Somatic': 'Somatic variant',
        'Germline': 'Germline variant',
        'Reference': 'Reference/wildtype',
        'Artifact': 'Artifact/technical error'
    }
    
    for filt, description in classification_filters.items():
        if filt not in new_header.filters:
            new_header.filters.add(filt, None, None, description)
    
    # Add RaVeX filter definitions (check if they already exist)
    if 'RaVeX_FILTER' not in new_header.filters:
        new_header.filters.add(
            'RaVeX_FILTER', None, None,
            'RaVeX filtering applied'
        )
    
    if 'RaVeX_FILTER' not in new_header.info:
        new_header.info.add(
            'RaVeX_FILTER', '.', 'String',
            'RaVeX filter reasons: semicolon-separated list of filter flags'
        )
    
    # Add individual filter flag INFO fields to prevent bcftools warnings
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
        if flag not in new_header.info:
            new_header.info.add(flag, '0', 'Flag', description)
    
    # Create output VCF
    output_vcf = pysam.VariantFile(output_path, 'w', header=new_header)
    
    # Re-open input with pysam for writing
    input_vcf_pysam = pysam.VariantFile(input_vcf_path)
    
    # Create a mapping of variants
    variant_map = {}
    for variant, filter_status, filter_list in filtered_variants:
        vkey = f"{variant.CHROM}:{variant.POS}:{variant.REF}:{variant.ALT[0]}"
        variant_map[vkey] = (filter_status, filter_list)
    
    # Get chromosome order for sorting
    chrom_order = {contig: idx for idx, contig in enumerate(new_header.contigs)}
    
    # Collect records to write
    records_to_write = []
    
    for record in input_vcf_pysam:
        vkey = f"{record.contig}:{record.pos}:{record.ref}:{record.alts[0]}"
        
        if vkey not in variant_map:
            continue
        
        filter_status, filter_list = variant_map[vkey]
        records_to_write.append((record, filter_status, filter_list))
    
    # Sort records
    def sort_key(item):
        record, _, _ = item
        chrom_idx = chrom_order.get(record.contig, 999999)
        return (chrom_idx, record.start)
    
    records_to_write.sort(key=sort_key)
    
    # Write sorted records
    for record, filter_status, filter_list in records_to_write:
        # Create new record with new header
        new_record = output_vcf.new_record(
            contig=record.contig,
            start=record.start,
            stop=record.stop,
            alleles=record.alleles,
            id=record.id,
            qual=record.qual
        )
        
        # Copy INFO fields (only if defined in header)
        for key in record.info:
            # Skip if not in output header
            if key not in new_header.info:
                continue
            try:
                value = record.info[key]
                # Skip None or empty values
                if value is not None and value != '':
                    new_record.info[key] = value
            except Exception:
                pass
        
        # Copy FORMAT and sample data (only if defined in header and not stripped)
        if not strip_format and record.format:
            for fmt_key in record.format.keys():
                # Skip if not in output header
                if fmt_key not in new_header.formats:
                    continue
                try:
                    for sample_idx, sample in enumerate(record.samples):
                        value = record.samples[sample_idx][fmt_key]
                        # Skip None values
                        if value is not None:
                            new_record.samples[sample_idx][fmt_key] = value
                except Exception:
                    pass
        
        # CRITICAL: Preserve original FILTER field (Somatic/Germline/Reference/Artifact)
        # Copy original FILTER values from input record
        if hasattr(record, 'filter') and record.filter:
            # pysam represents FILTER as VariantRecordFilter object
            for filt in record.filter:
                if filt != 'PASS' and filt in new_header.filters:
                    new_record.filter.add(filt)
        
        # Set RaVeX filter INFO tag (not FILTER field)
        if filter_status == "PASS":
            new_record.info['RaVeX_FILTER'] = "PASS"
        else:
            # Store filter reasons in INFO field only (use comma separator for multi-valued String fields)
            new_record.info['RaVeX_FILTER'] = ",".join(filter_list)
            # Also set individual filter flags
            for flag in filter_list:
                if flag in new_header.info:
                    new_record.info[flag] = True
        
        output_vcf.write(new_record)
    
    # Close files
    output_vcf.close()
    input_vcf_pysam.close()
    input_vcf.close()


def main():
    args = argparser()

    # Auto-generate stripped output path if not provided
    if not args.output_stripped:
        from pathlib import Path
        output_path = Path(args.output)
        args.output_stripped = str(output_path.parent / f"{output_path.stem}.stripped.vcf.gz")

    print("=" * 80)
    print("Consensus VCF Filtering with RaVeX Rules")
    print("=" * 80)
    print(f"Input:             {args.input}")
    print(f"Output (standard): {args.output}")
    print(f"Output (stripped): {args.output_stripped}")

    # Open genome if provided (optional)
    genome = None
    if args.ref:
        try:
            genome = pysam.FastaFile(args.ref)
        except Exception:
            print("Warning: Failed to open reference; continuing without it")
            genome = None

    # Read and filter variants using cyvcf2
    vcf_in = VCF(args.input)
    filtered_variants = apply_filters(vcf_in, args, genome)

    # Generate standard output (all variants, FORMAT preserved)
    print("\nGenerating standard output...")
    write_filtered_vcf(filtered_variants, args.input, args.output, strip_format=False)
    print(f"✓ Standard output written to: {args.output}")

    # Re-read for stripped output
    vcf_in = VCF(args.input)
    filtered_variants_stripped = apply_filters(vcf_in, args, genome)

    # Filter multiallelic variants for stripped output
    if args.filter_multiallelic:
        filtered_variants_stripped = [
            (v, s, f) for v, s, f in filtered_variants_stripped
            if not is_multiallelic(v.REF, v.ALT)
        ]

    # Generate stripped output (multiallelic-filtered, FORMAT removed)
    print("\nGenerating stripped output...")
    write_vcf_stripped(filtered_variants_stripped, args.input, args.output_stripped, use_cyvcf2=True)
    print(f"✓ Stripped output written to: {args.output_stripped}")

    # Index outputs
    print("\nIndexing outputs...")
    import subprocess
    for vcf_path in [args.output, args.output_stripped]:
        try:
            subprocess.run(['tabix', '-p', 'vcf', vcf_path], check=True, capture_output=True)
            print(f"  ✓ Indexed: {vcf_path}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            print(f"  ⚠ Warning: Failed to index {vcf_path}")

    print("\n" + "=" * 80)
    print("Filtering completed successfully!")
    print("=" * 80)


if __name__ == "__main__":
    main()
