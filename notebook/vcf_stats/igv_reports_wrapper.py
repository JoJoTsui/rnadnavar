"""
IGV-Reports based visualization for tiered rescue variants.

This module provides wrappers around igv-reports to generate professional,
interactive HTML reports for variant-BAM visualization, organized by tier.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import subprocess
import tempfile
import shutil
import pandas as pd


def create_subset_vcf(
    variants_df: pd.DataFrame,
    original_vcf: Path,
    output_vcf: Path,
) -> Path:
    """
    Create a subset VCF containing only the sampled variants.
    
    Uses bcftools to efficiently extract variant regions from the original VCF.
    The subset VCF preserves all INFO fields and will be indexed.
    
    Args:
        variants_df: DataFrame with columns: chrom, pos, ref, alt
        original_vcf: Path to original rescue VCF (bgzipped and indexed)
        output_vcf: Path where subset VCF will be written
        
    Returns:
        Path to output VCF file (bgzipped)
        
    Raises:
        RuntimeError: If bcftools is not available or VCF operations fail
    """
    output_vcf = Path(output_vcf)
    
    # Check if bcftools is available
    try:
        subprocess.run(
            ["bcftools", "--version"],
            capture_output=True,
            check=True,
            text=True
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise RuntimeError("bcftools is required but not found. Install with: conda install -c bioconda bcftools")
    
    # Create a BED file with variant regions (±10bp for context)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as bed_file:
        for _, row in variants_df.iterrows():
            chrom = row['chrom']
            pos = int(row['pos'])
            start = max(0, pos - 10)
            end = pos + 10
            bed_file.write(f"{chrom}\t{start}\t{end}\n")
        bed_path = Path(bed_file.name)
    
    try:
        # Use bcftools to extract regions
        extract_cmd = [
            "bcftools",
            "view",
            "-R", str(bed_path),
            "-o", str(output_vcf),
            "-O", "z",  # Output bgzipped VCF
            str(original_vcf)
        ]
        
        result = subprocess.run(
            extract_cmd,
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"bcftools view failed: {result.stderr}")
        
        # Index the output VCF
        subprocess.run(
            ["bcftools", "index", "-t", str(output_vcf)],
            capture_output=True,
            check=False
        )
        
        return output_vcf
        
    finally:
        # Clean up temporary BED file
        if bed_path.exists():
            bed_path.unlink()


def check_igv_reports_available() -> bool:
    """Check if igv-reports is installed and available."""
    try:
        import igvreports
        return True
    except ImportError:
        return False


def get_alignment_index_path(alignment_path: Path) -> Path:
    """
    Get the expected index file path for a BAM or CRAM alignment file.
    
    Args:
        alignment_path: Path to BAM or CRAM file
        
    Returns:
        Path to expected index file (.bai for BAM, .crai for CRAM)
    """
    alignment_path = Path(alignment_path)
    if str(alignment_path).endswith(".cram"):
        return Path(str(alignment_path) + ".crai")
    else:
        return Path(str(alignment_path) + ".bai")


def generate_igv_report(
    subset_vcf: Path,
    bam_files: Dict[str, Path],
    ref_fasta: Path,
    output_dir: Path,
    report_title: str = "Variant IGV Report",
) -> Optional[Path]:
    """
    Generate an interactive IGV-reports HTML report.
    
    Prerequisites:
    - igv-reports must be installed: pip install igv-reports
    - Reference FASTA must be indexed (.fai file present)
    - BAM/CRAM files must be indexed (.bai for BAM, .crai for CRAM)
    
    Args:
        subset_vcf: Path to subset VCF with variants to visualize
        bam_files: Dictionary mapping sample name to BAM/CRAM path
        ref_fasta: Path to reference FASTA (must be indexed)
        output_dir: Directory where report HTML will be saved
        report_title: Title for the HTML report
        
    Returns:
        Path to generated HTML report, or None if generation failed
        
    Raises:
        RuntimeError: If igv-reports is not installed
    """
    try:
        from igvreports.igvreports import create_report
    except ImportError:
        raise RuntimeError(
            "igv-reports is not installed. Install manually with: pip install igv-reports"
        )
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Prepare tracks for IGV report
    tracks = []
    
    # Add BAM tracks in preferred order
    track_order = ["DNA_NORMAL", "DNA_TUMOR", "RNA_TUMOR"]
    for sample_name in track_order:
        if sample_name in bam_files:
            bam_path = bam_files[sample_name]
            if bam_path.exists():
                tracks.append({
                    "url": str(bam_path),
                    "name": sample_name,
                    "format": "bam" if str(bam_path).endswith(".bam") else "cram"
                })
    
    # Generate the report
    try:
        report_path = output_dir / "index.html"
        
        # igvreports create_report parameters
        create_report(
            vcf=str(subset_vcf),
            bam=tracks,
            fasta=str(ref_fasta),
            outputFile=str(report_path),
            title=report_title,
            flanking=100,  # bp of context on each side
        )
        
        print(f"✓ IGV report generated: {report_path}")
        return report_path
        
    except Exception as e:
        print(f"✗ Error generating IGV report: {e}")
        return None


def create_tier_summary_html(
    tier_reports: Dict[str, Path],
    tier_descriptions: Dict[str, str],
    output_file: Path,
) -> bool:
    """
    Create an HTML landing page that links to per-tier reports.
    
    Args:
        tier_reports: Dictionary mapping tier (e.g., 'T1') to report HTML path
        tier_descriptions: Dictionary mapping tier to description text
        output_file: Path where landing page HTML will be written
        
    Returns:
        True if successful, False otherwise
    """
    try:
        html_parts = [
            "<!DOCTYPE html>",
            "<html>",
            "<head>",
            "  <meta charset='utf-8'>",
            "  <meta name='viewport' content='width=device-width, initial-scale=1'>",
            "  <title>Variant Tiering & IGV Reports</title>",
            "  <style>",
            "    body { font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }",
            "    .container { max-width: 1000px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; }",
            "    h1 { color: #333; border-bottom: 3px solid #0066cc; padding-bottom: 10px; }",
            "    .tier-section { margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }",
            "    .tier-section h2 { color: #0066cc; margin-top: 0; }",
            "    .tier-id { display: inline-block; background: #0066cc; color: white; padding: 3px 8px; border-radius: 3px; margin-right: 10px; }",
            "    .tier-desc { color: #666; margin: 10px 0; }",
            "    .tier-link { display: inline-block; margin-top: 10px; padding: 10px 20px; background: #0066cc; color: white; text-decoration: none; border-radius: 4px; }",
            "    .tier-link:hover { background: #0052a3; }",
            "    .no-variants { color: #999; font-style: italic; }",
            "  </style>",
            "</head>",
            "<body>",
            "<div class='container'>",
            "  <h1>Variant Tiering & IGV-Reports</h1>",
            "  <p>Interactive BAM alignment visualization organized by variant tier confidence.</p>",
        ]
        
        # Add tier sections
        tier_order = [f"T{i}" for i in range(1, 9)]
        for tier in tier_order:
            if tier in tier_reports:
                report_path = tier_reports[tier]
                description = tier_descriptions.get(tier, "")
                
                # Calculate relative path for link
                rel_path = Path(report_path).relative_to(output_file.parent)
                
                html_parts.extend([
                    f"  <div class='tier-section'>",
                    f"    <h2><span class='tier-id'>{tier}</span></h2>",
                    f"    <p class='tier-desc'>{description}</p>",
                    f"    <a href='{rel_path}' class='tier-link'>View IGV Report</a>",
                    f"  </div>",
                ])
        
        html_parts.extend([
            "  <p style='color: #999; margin-top: 30px; padding-top: 20px; border-top: 1px solid #ddd;'>",
            "    Reports generated using <a href='https://github.com/igvteam/igv-reports' target='_blank'>igv-reports</a>",
            "  </p>",
            "</div>",
            "</body>",
            "</html>",
        ])
        
        html_content = "\n".join(html_parts)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(html_content)
        
        return True
        
    except Exception as e:
        print(f"✗ Error creating tier summary HTML: {e}")
        return False


def visualize_rescue_variants_with_igvreports(
    variants_df: pd.DataFrame,
    original_vcf: Path,
    bam_files: Dict[str, Path],
    ref_fasta: Path,
    output_dir: Path,
    k_per_tier: int = 3,
) -> Dict[str, Path]:
    """
    End-to-end function to generate IGV-reports for tiered rescue variants.
    
    Prerequisites:
    - igv-reports must be installed: pip install igv-reports
    - Reference FASTA must be indexed (.fai file present)
    - BAM/CRAM files must be indexed (.bai for BAM, .crai for CRAM)
    
    Workflow:
    1. Group variants by category and tier
    2. Sample k representatives per tier per category
    3. Create subset VCF
    4. Generate IGV report per tier per category
    5. Create landing page
    
    Args:
        variants_df: DataFrame with tiered variants (must have 'tier' and 'filter_category' columns)
        original_vcf: Path to original rescue VCF
        bam_files: Sample name -> BAM/CRAM path mapping
        ref_fasta: Path to indexed reference FASTA
        output_dir: Base output directory for all reports
        k_per_tier: Number of representative variants per tier per category
        
    Returns:
        Dictionary mapping tier to report HTML path
        
    Raises:
        RuntimeError: If igv-reports not installed or prerequisites not met
    """
    from .tiering import sample_tier_representatives
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    tier_reports = {}
    tier_descriptions = {
        "T1": "DNA consensus (≥2 callers) + RNA consensus (≥2 callers) - Strongest, both modalities agree",
        "T2": "DNA consensus (≥2) + RNA single caller - DNA strong with RNA corroboration",
        "T3": "DNA consensus (≥2) only - DNA consensus without RNA detection",
        "T4": "DNA single + RNA support - Single DNA caller with RNA corroboration",
        "T5": "DNA single only - Single DNA caller, no RNA support",
        "T6": "RNA consensus (≥2) only - RNA consensus without DNA detection",
        "T7": "RNA single only - Single RNA caller, no DNA support",
        "T8": "No caller support - Edge case, extremely rare",
    }
    
    # Group by category and tier
    categories = sorted(variants_df["filter_category"].unique())
    
    for category in categories:
        cat_variants = variants_df[variants_df["filter_category"] == category]
        print(f"\n  {category} ({len(cat_variants)} variants):")
        
        # Within this category, group by tier
        cat_tier_groups = cat_variants.groupby('tier')
        
        for tier in [f"T{i}" for i in range(1, 9)]:
            if tier not in cat_variants['tier'].values:
                continue
            
            # Get variants for this tier within this category
            tier_df = cat_tier_groups.get_group(tier)
            count = len(tier_df)
            pct = (count / len(cat_variants)) * 100
            print(f"    {tier}: {count:6d} variants ({pct:5.1f}%)")
            
            # Sample representatives
            if len(tier_df) <= k_per_tier:
                sampled = tier_df
            else:
                # Random sample
                sampled = tier_df.sample(n=k_per_tier, random_state=42)
            
            # Create category/tier output directory
            cat_tier_dir = output_dir / category / tier
            cat_tier_dir.mkdir(parents=True, exist_ok=True)
            
            # Create subset VCF
            subset_vcf = cat_tier_dir / f"{category}_{tier}_subset.vcf.gz"
            try:
                create_subset_vcf(sampled, original_vcf, subset_vcf)
            except Exception as e:
                print(f"      ✗ Failed to create subset VCF: {e}")
                continue
            
            # Generate IGV report
            try:
                report_path = generate_igv_report(
                    subset_vcf,
                    bam_files,
                    ref_fasta,
                    cat_tier_dir,
                    report_title=f"{category} - Tier {tier} Variants"
                )
                
                if report_path:
                    # Use category/tier as key for hierarchical organization
                    key = f"{category}/{tier}"
                    tier_reports[key] = report_path
            except RuntimeError as e:
                # igv-reports not available - re-raise
                raise e
            except Exception as e:
                print(f"      ✗ Error generating IGV report: {e}")
    
    # Create landing page
    landing_page = output_dir / "index.html"
    if create_tier_summary_html(tier_reports, tier_descriptions, landing_page):
        print(f"\n✓ Created landing page: {landing_page}")
    
    return tier_reports
