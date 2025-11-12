#!/usr/bin/env python3
"""
VCF Consensus Script
Mirrors the logic of run_consensus.R but works natively with VCF files
Uses cyvcf2 for fast reading and pysam for writing
"""
import argparse
import sys
from pathlib import Path
from cyvcf2 import VCF
import pysam


def argparser():
    parser = argparse.ArgumentParser(description="Find consensus variants across multiple VCF files")
    parser.add_argument("--input_dir", required=True, help="Directory containing input VCF files")
    parser.add_argument("--out_prefix", required=True, help="Prefix for output files")
    parser.add_argument("--thr", type=int, default=2, help="Number of callers required for consensus")
    parser.add_argument("--cpu", type=int, default=1, help="Number of CPU cores")
    parser.add_argument("--id", help="Sample ID")
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
    # Fallback
    return name.split('.')[1] if '.' in name else 'unknown'


def variant_key(variant):
    """Create unique key for variant (works with both cyvcf2 and pysam)"""
    if hasattr(variant, 'CHROM'):  # cyvcf2
        return f"{variant.CHROM}:{variant.POS}:{variant.REF}:{variant.ALT[0]}"
    else:  # pysam
        return f"{variant.chrom}:{variant.pos}:{variant.ref}:{variant.alts[0]}"


def find_overlaps(vcf_files):
    """Find overlapping variants across VCF files using cyvcf2 for speed"""
    print("- Finding overlaps")
    
    # Store variants by caller
    variants_by_caller = {}
    variant_details = {}
    
    for caller, vcf_path in vcf_files.items():
        print(f"  - Processing {caller}")
        vcf = VCF(vcf_path)
        variants_by_caller[caller] = set()
        
        for variant in vcf:
            vkey = variant_key(variant)
            variants_by_caller[caller].add(vkey)
            
            # Store variant details (only first occurrence)
            if vkey not in variant_details:
                variant_details[vkey] = {
                    'CHROM': variant.CHROM,
                    'POS': variant.POS,
                    'REF': variant.REF,
                    'ALT': variant.ALT[0],
                    'callers': [],
                    'filters': [],
                    'source_caller': caller,  # Remember which caller to get full record from
                    'source_file': vcf_path
                }
            
            variant_details[vkey]['callers'].append(caller)
            variant_details[vkey]['filters'].append(variant.FILTER or 'PASS')
    
    return variants_by_caller, variant_details


def determine_consensus(variant_details, threshold, callers_list):
    """Determine which variants pass consensus threshold"""
    print(f"- Determining consensus (threshold={threshold})")
    
    snvs = 0
    indels = 0
    consensus_variants = {}
    
    for vkey, details in variant_details.items():
        n_callers = len(details['callers'])
        is_snv = len(details['REF']) == 1 and len(details['ALT']) == 1
        
        # SNVs require exact match across callers
        if is_snv and n_callers >= threshold:
            consensus_variants[vkey] = details
            snvs += 1
        # Indels just need overlap
        elif not is_snv and n_callers >= 1:
            consensus_variants[vkey] = details
            indels += 1
    
    print(f"- There are {snvs:,} SNVs that are consensus")
    print(f"- There are {indels:,} indels that are consensus")
    
    return consensus_variants


def write_consensus_vcf(consensus_variants, vcf_files, out_file, callers_list):
    """Write consensus VCF file using pysam for proper INFO field handling"""
    print("- Writing consensus VCF")
    
    # Open template VCF with pysam
    template_path = list(vcf_files.values())[0]
    template_vcf = pysam.VariantFile(template_path)
    
    # Create a new header based on template
    new_header = template_vcf.header.copy()
    
    # Add custom INFO fields to header
    new_header.info.add(
        'callers', '.', 'String',
        'Variant callers that called this mutation, separated by |'
    )
    new_header.info.add(
        'filters', '.', 'String',
        'Filters provided by each variant caller, separated by |'
    )
    new_header.info.add(
        'consensus_filter', '1', 'String',
        'PASS if 50% or more of the callers give PASS, otherwise FAIL'
    )
    
    # Add FAIL filter to header
    new_header.filters.add(
        'FAIL', None, None,
        'More than half the callers did not give a PASS'
    )
    
    # Create output VCF with new header
    output_vcf = pysam.VariantFile(out_file, 'w', header=new_header)
    
    # Close template
    template_vcf.close()
    
    # Create lookup of consensus variants
    consensus_keys = set(consensus_variants.keys())
    
    # Open all VCF files with pysam and collect consensus records
    pysam_vcfs = {caller: pysam.VariantFile(path) for caller, path in vcf_files.items()}
    
    # Collect all consensus records with their details
    records_to_write = []
    seen_variants = set()
    
    for caller, vcf in pysam_vcfs.items():
        for record in vcf:
            vkey = variant_key(record)
            
            # Skip if not in consensus or already collected
            if vkey not in consensus_keys or vkey in seen_variants:
                continue
            
            details = consensus_variants[vkey]
            records_to_write.append((record, details))
            seen_variants.add(vkey)
    
    # Close input VCFs
    for vcf in pysam_vcfs.values():
        vcf.close()
    
    # Sort records by chromosome and position
    # Get chromosome order from header
    chrom_order = {contig: idx for idx, contig in enumerate(new_header.contigs)}
    
    def sort_key(item):
        record, _ = item
        chrom_idx = chrom_order.get(record.contig, 999999)
        return (chrom_idx, record.start)
    
    records_to_write.sort(key=sort_key)
    
    # Write sorted records
    for record, details in records_to_write:
        # Create a new record with the new header
        new_record = output_vcf.new_record(
            contig=record.contig,
            start=record.start,
            stop=record.stop,
            alleles=record.alleles,
            id=record.id,
            qual=record.qual
        )
        
        # Copy INFO fields from original record
        for key in record.info:
            try:
                new_record.info[key] = record.info[key]
            except Exception:
                pass
        
        # Copy FORMAT and sample data
        if record.format:
            for fmt_key in record.format.keys():
                try:
                    new_record.samples[0][fmt_key] = record.samples[0][fmt_key]
                except Exception:
                    pass
        
        # Determine consensus filter
        pass_count = sum(1 for f in details['filters'] if f == 'PASS')
        if pass_count / len(details['filters']) >= 0.5:
            consensus_filter = 'PASS'
            # Leave filter empty for PASS
        else:
            consensus_filter = 'FAIL'
            new_record.filter.add('FAIL')
        
        # Add consensus INFO fields
        new_record.info['callers'] = '|'.join(details['callers'])
        new_record.info['filters'] = '|'.join(details['filters'])
        new_record.info['consensus_filter'] = consensus_filter
        
        # Write the new record
        output_vcf.write(new_record)
    
    # Close output VCF
    output_vcf.close()
    
    print(f"- Output in: {out_file}")


def main():
    args = argparser()
    
    # Find VCF files
    input_dir = Path(args.input_dir)
    vcf_files = {}
    
    for vcf_path in input_dir.glob("*.vcf*"):
        if vcf_path.suffix in ['.vcf', '.gz'] or vcf_path.name.endswith('.vcf.gz'):
            caller = get_caller_name(str(vcf_path))
            vcf_files[caller] = str(vcf_path)
    
    if not vcf_files:
        print("ERROR: No VCF files found")
        sys.exit(1)
    
    print(f"Found {len(vcf_files)} VCF files: {', '.join(vcf_files.keys())}")
    
    # Find overlaps (using cyvcf2 for speed)
    variants_by_caller, variant_details = find_overlaps(vcf_files)
    
    # Determine consensus
    consensus_variants = determine_consensus(variant_details, args.thr, list(vcf_files.keys()))
    
    # Write output (using pysam for proper VCF writing)
    out_file = f"{args.out_prefix}.vcf"
    write_consensus_vcf(consensus_variants, vcf_files, out_file, list(vcf_files.keys()))
    
    print(f"Done! Total consensus variants: {len(consensus_variants):,}")


if __name__ == "__main__":
    main()
