"""
IGV-Reports based visualization for category-tiered rescue variants.

This module provides a comprehensive interface for generating professional,
interactive HTML reports using the igv-reports command-line tool (create_report).

Key Design Decisions:
---------------------

1. **Command-Line API (NOT Python API)**
   - Uses subprocess to call 'create_report' command-line script
   - Reason: Generates self-contained HTML with all data embedded
   - No external server required - reports work offline
   - Output HTML is stable and reliable

2. **Category-Tiered Organization**
   - VCF files grouped hierarchically: Category → Tier
   - Each category/tier combo gets its own subset VCF and report
   - Hierarchical landing page with category sections
   - Example: Somatic/T1, Somatic/T2, Germline/T1, etc.

3. **Indexed Files Requirement**
   - VCF: Must be bgzip-compressed and indexed with .tbi (tabix)
   - Reference: Must have .fai index (samtools faidx)
   - BAM/CRAM: Must be indexed (.bai for BAM, .crai for CRAM)
   - Indexes enable efficient variant subsetting and visualization

Workflow:
---------

For each category-tier combination:
  1. Create BED file with variant regions
  2. Use bcftools to extract subset VCF (-R flag)
  3. Compress with bgzip
  4. Index with tabix
  5. Call create_report command with proper arguments
  6. Generate self-contained HTML report

Landing Page:
  - Hierarchical HTML with categories as sections
  - Shows variant counts per tier
  - Links to individual reports with relative paths
  - Professional styling with modern UI

Usage Example:
--------------

    from vcf_stats import organize_by_category_tier
    
    tier_reports = organize_by_category_tier(
        variants_df=tiered_variants,
        original_vcf=Path("rescue.vcf.gz"),
        bam_files={
            "DNA_NORMAL": Path("dna_normal.bam"),
            "DNA_TUMOR": Path("dna_tumor.bam"),
            "RNA_TUMOR": Path("rna_tumor.cram"),
        },
        ref_fasta=Path("reference.fasta"),
        output_dir=Path("igv_reports"),
        k_per_tier=None  # Use all variants
    )
    
    # Returns: {"Somatic/T1": (8, Path(...)), "Somatic/T2": (7, Path(...))}
    # Reports at: igv_reports/Somatic/T1/report.html, etc.
    # Landing page at: igv_reports/index.html

Prerequisites Installation:
----------------------------

    # Install igv-reports command-line tool
    pip install igv-reports
    
    # Install required bioinformatics tools
    conda install -c bioconda bcftools tabix samtools

    # Verify installation
    create_report --help
    bcftools --version
    tabix -h
    samtools --version
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
        subprocess.run(
            ["create_report", "--help"],
            capture_output=True,
            check=True,
            text=True
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
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


def generate_igv_report_subprocess(
    subset_vcf: Path,
    bam_files: Dict[str, Path],
    ref_fasta: Path,
    output_dir: Path,
    report_title: str = "Variant IGV Report",
    info_columns: Optional[List[str]] = None,
) -> Optional[Path]:
    """
    Generate an interactive IGV-reports HTML report using command-line create_report.
    
    This is the correct method: uses the igv-reports command-line script (create_report)
    which generates self-contained HTML with all data embedded as blobs.
    
    Prerequisites:
    - igv-reports must be installed: pip install igv-reports
    - Reference FASTA must be indexed (.fai file present)
    - VCF file must be indexed (.tbi file present)
    - BAM/CRAM files must be indexed (.bai for BAM, .crai for CRAM)
    
    Args:
        subset_vcf: Path to subset VCF with variants to visualize (must be bgzip+tabix indexed)
        bam_files: Dictionary mapping sample name to BAM/CRAM path
        ref_fasta: Path to reference FASTA (must be indexed with .fai)
        output_dir: Directory where report HTML will be saved
        report_title: Title for the HTML report
        info_columns: List of VCF INFO columns to display (optional)
        
    Returns:
        Path to generated HTML report, or None if generation failed
        
    Raises:
        RuntimeError: If create_report command not found
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if create_report is available
    try:
        subprocess.run(
            ["create_report", "--help"],
            capture_output=True,
            check=True,
            text=True
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise RuntimeError(
            "create_report command not found. Install igv-reports with: pip install igv-reports"
        )
    
    try:
        report_path = output_dir / "report.html"
        
        # Build command: create_report sites.vcf.gz --fasta ref.fa --tracks bam1 bam2 --output report.html --title "Title"
        cmd = [
            "create_report",
            str(subset_vcf),
            "--fasta", str(ref_fasta),
            "--output", str(report_path),
            "--title", report_title,
        ]
        
        # Add tracks in preferred order
        track_order = ["DNA_NORMAL", "DNA_TUMOR", "RNA_TUMOR"]
        tracks = []
        for sample_name in track_order:
            if sample_name in bam_files:
                bam_path = bam_files[sample_name]
                if Path(bam_path).exists():
                    tracks.append(str(bam_path))
        
        if tracks:
            cmd.extend(["--tracks"] + tracks)
        
        # Add info columns if specified
        if info_columns:
            cmd.extend(["--info-columns"] + info_columns)
        
        # Run create_report
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode != 0:
            print(f"✗ create_report error: {result.stderr}")
            return None
        
        if report_path.exists():
            print(f"✓ IGV report generated: {report_path}")
            return report_path
        else:
            print(f"✗ Report file not created at {report_path}")
            return None
            
    except Exception as e:
        print(f"✗ Error generating IGV report: {e}")
        return None



def create_navigation_page(
    tier_reports: Dict[str, Tuple[int, Path]],
    output_file: Path,
) -> bool:
    """
    Create hierarchical landing page for category-tiered reports.
    
    Args:
        tier_reports: Dictionary mapping "Category/Tier" to (variant_count, report_path)
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
            "  <title>Variant Review by Category and Tier</title>",
            "  <style>",
            "    body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 0; background: #f0f2f5; }",
            "    .container { max-width: 1200px; margin: 0 auto; background: white; min-height: 100vh; }",
            "    .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 40px 20px; }",
            "    .header h1 { margin: 0; font-size: 32px; }",
            "    .header p { margin: 10px 0 0 0; opacity: 0.95; }",
            "    .content { padding: 30px 20px; }",
            "    .category-section { margin-bottom: 40px; }",
            "    .category-title { display: flex; align-items: center; margin-bottom: 20px; padding-bottom: 15px; border-bottom: 2px solid #667eea; }",
            "    .category-name { font-size: 24px; font-weight: bold; color: #333; }",
            "    .category-stats { margin-left: auto; font-size: 14px; color: #666; }",
            "    .tier-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 20px; }",
            "    .tier-card { border: 1px solid #e0e0e0; border-radius: 8px; overflow: hidden; transition: all 0.3s ease; }",
            "    .tier-card:hover { box-shadow: 0 4px 12px rgba(102, 126, 234, 0.15); transform: translateY(-2px); }",
            "    .tier-header { background: #f8f9fa; padding: 15px; border-bottom: 2px solid #667eea; }",
            "    .tier-id { display: inline-block; background: #667eea; color: white; padding: 4px 10px; border-radius: 4px; font-weight: bold; }",
            "    .tier-count { float: right; color: #667eea; font-weight: bold; }",
            "    .tier-body { padding: 15px; }",
            "    .tier-desc { color: #666; font-size: 13px; line-height: 1.5; margin-bottom: 15px; }",
            "    .view-btn { display: inline-block; width: 100%; text-align: center; padding: 10px; background: #667eea; color: white; text-decoration: none; border-radius: 4px; transition: background 0.2s; }",
            "    .view-btn:hover { background: #764ba2; }",
            "    .empty-category { color: #999; font-style: italic; padding: 20px; text-align: center; }",
            "    .footer { background: #f8f9fa; padding: 20px; text-align: center; color: #666; border-top: 1px solid #e0e0e0; font-size: 13px; }",
            "  </style>",
            "</head>",
            "<body>",
            "<div class='container'>",
            "  <div class='header'>",
            "    <h1>Variant Review by Category and Tier</h1>",
            "    <p>Interactive BAM alignment visualization organized by variant biological classification and caller confidence</p>",
            "  </div>",
            "  <div class='content'>",
        ]
        
        # Group reports by category
        category_map = {}
        for key, (count, path) in tier_reports.items():
            parts = key.split('/')
            if len(parts) == 2:
                category, tier = parts
                if category not in category_map:
                    category_map[category] = {}
                category_map[category][tier] = (count, path)
        
        # Tier descriptions
        tier_descriptions = {
            "T1": "DNA consensus (≥2 callers) + RNA consensus (≥2 callers)",
            "T2": "DNA consensus (≥2) + RNA single caller",
            "T3": "DNA consensus (≥2) only",
            "T4": "DNA single + RNA support",
            "T5": "DNA single only",
            "T6": "RNA consensus (≥2) only",
            "T7": "RNA single only",
            "T8": "No caller support",
        }
        
        # Render by category
        for category in sorted(category_map.keys()):
            tier_map = category_map[category]
            total_variants = sum(count for count, _ in tier_map.values())
            
            html_parts.append("    <div class='category-section'>")
            html_parts.append(f"      <div class='category-title'>")
            html_parts.append(f"        <span class='category-name'>{category}</span>")
            html_parts.append(f"        <span class='category-stats'>{total_variants:,} variants ({len(tier_map)} tiers)</span>")
            html_parts.append(f"      </div>")
            html_parts.append(f"      <div class='tier-grid'>")
            
            for tier in [f"T{i}" for i in range(1, 9)]:
                if tier not in tier_map:
                    continue
                
                count, path = tier_map[tier]
                rel_path = Path(path).relative_to(output_file.parent)
                description = tier_descriptions.get(tier, "")
                pct = (count / total_variants) * 100
                
                html_parts.append(f"        <div class='tier-card'>")
                html_parts.append(f"          <div class='tier-header'>")
                html_parts.append(f"            <span class='tier-id'>{tier}</span>")
                html_parts.append(f"            <span class='tier-count'>{count} variants ({pct:.1f}%)</span>")
                html_parts.append(f"          </div>")
                html_parts.append(f"          <div class='tier-body'>")
                html_parts.append(f"            <p class='tier-desc'>{description}</p>")
                html_parts.append(f"            <a href='{rel_path}' class='view-btn'>View Report</a>")
                html_parts.append(f"          </div>")
                html_parts.append(f"        </div>")
            
            html_parts.append(f"      </div>")
            html_parts.append(f"    </div>")
        
        html_parts.extend([
            "  </div>",
            "  <div class='footer'>",
            "    <p>Reports generated using <a href='https://github.com/igvteam/igv-reports' target='_blank'>igv-reports</a> command-line interface</p>",
            "  </div>",
            "</div>",
            "</body>",
            "</html>",
        ])
        
        html_content = "\n".join(html_parts)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(html_content)
        
        return True
        
    except Exception as e:
        print(f"✗ Error creating navigation page: {e}")
        return False


def organize_by_category_tier(
    variants_df: pd.DataFrame,
    original_vcf: Path,
    bam_files: Dict[str, Path],
    ref_fasta: Path,
    output_dir: Path,
    k_per_tier: Optional[int] = None,
) -> Dict[str, Tuple[int, Path]]:
    """
    Main orchestration function: generate category-tiered IGV reports.
    
    This function creates a hierarchical directory structure (Category/Tier)
    with indexed VCF subsets and interactive IGV-reports for each combination.
    
    Prerequisites:
    - igv-reports must be installed: pip install igv-reports
    - create_report command must be available in PATH
    - Reference FASTA must be indexed (.fai file present)
    - All BAM/CRAM files must be indexed (.bai or .crai)
    - Original VCF must be indexed (.tbi)
    
    Args:
        variants_df: DataFrame with tiered variants (must have 'tier' and 'filter_category' columns)
        original_vcf: Path to original rescue VCF (bgzip+tabix indexed)
        bam_files: Sample name -> BAM/CRAM path mapping
        ref_fasta: Path to indexed reference FASTA
        output_dir: Base output directory for all reports
        k_per_tier: Optional max variants per tier (None = all variants)
        
    Returns:
        Dictionary mapping "Category/Tier" to (variant_count, report_path)
        
    Raises:
        RuntimeError: If create_report not found or prerequisites not met
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check prerequisites
    if not check_igv_reports_available():
        raise RuntimeError(
            "igv-reports is not installed. Install with: pip install igv-reports"
        )
    
    # Verify create_report is available
    try:
        subprocess.run(
            ["create_report", "--help"],
            capture_output=True,
            check=True,
            text=True,
            timeout=5
        )
    except:
        raise RuntimeError(
            "create_report command not found. Install igv-reports with: pip install igv-reports"
        )
    
    # Verify reference is indexed
    ref_fasta = Path(ref_fasta)
    ref_index = Path(str(ref_fasta) + ".fai")
    if not ref_index.exists():
        raise RuntimeError(f"Reference FASTA not indexed: {ref_index}")
    
    # Verify VCF is indexed
    vcf_index = Path(str(original_vcf) + ".tbi")
    if not vcf_index.exists():
        raise RuntimeError(f"VCF not indexed: {vcf_index}")
    
    tier_reports = {}
    
    # Group by category
    categories = sorted(variants_df["filter_category"].unique())
    
    for category in categories:
        cat_variants = variants_df[variants_df["filter_category"] == category]
        print(f"\nProcessing category: {category} ({len(cat_variants)} variants)")
        
        # Within this category, process each tier
        for tier in sorted(cat_variants["tier"].unique()):
            tier_df = cat_variants[cat_variants["tier"] == tier]
            count = len(tier_df)
            pct = (count / len(cat_variants)) * 100
            
            print(f"  {tier}: {count:6d} variants ({pct:5.1f}%)", end=" → ")
            
            # Sample if needed
            if k_per_tier and len(tier_df) > k_per_tier:
                sampled = tier_df.sample(n=k_per_tier, random_state=42)
                print(f"sampled to {k_per_tier}", end=" → ")
            else:
                sampled = tier_df
            
            # Create category/tier output directory
            cat_tier_dir = output_dir / category / tier
            cat_tier_dir.mkdir(parents=True, exist_ok=True)
            
            # Create subset VCF
            subset_vcf = cat_tier_dir / f"{category}_{tier}_subset.vcf.gz"
            try:
                print("creating VCF", end=" → ")
                create_subset_vcf(sampled, original_vcf, subset_vcf)
            except Exception as e:
                print(f"✗ Failed to create subset VCF: {e}")
                continue
            
            # Generate IGV report
            try:
                print("generating report", end=" → ")
                report_path = generate_igv_report_subprocess(
                    subset_vcf,
                    bam_files,
                    ref_fasta,
                    cat_tier_dir,
                    report_title=f"{category} - Tier {tier} Variants"
                )
                
                if report_path:
                    key = f"{category}/{tier}"
                    tier_reports[key] = (count, report_path)
                    print("✓")
                else:
                    print("✗ Report not generated")
            except Exception as e:
                print(f"✗ {e}")
    
    # Create landing page
    print(f"\nGenerating hierarchical landing page...")
    landing_page = output_dir / "index.html"
    if create_navigation_page(tier_reports, landing_page):
        print(f"✓ Created landing page: {landing_page}")
    
    return tier_reports


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
    
    This is a wrapper that calls organize_by_category_tier() and maintains
    backward compatibility with the previous API.
    
    Prerequisites:
    - igv-reports must be installed: pip install igv-reports
    - Reference FASTA must be indexed (.fai file present)
    - BAM/CRAM files must be indexed (.bai for BAM, .crai for CRAM)
    
    Args:
        variants_df: DataFrame with tiered variants (must have 'tier' and 'filter_category' columns)
        original_vcf: Path to original rescue VCF
        bam_files: Sample name -> BAM/CRAM path mapping
        ref_fasta: Path to indexed reference FASTA
        output_dir: Base output directory for all reports
        k_per_tier: Number of representative variants per tier per category (used for sampling)
        
    Returns:
        Dictionary mapping "Category/Tier" to report HTML path
        
    Raises:
        RuntimeError: If igv-reports not installed or prerequisites not met
    """
    # Call new implementation with category-tier organization
    tier_reports_with_counts = organize_by_category_tier(
        variants_df,
        original_vcf,
        bam_files,
        ref_fasta,
        output_dir,
        k_per_tier=k_per_tier
    )
    
    # Return in backward-compatible format (key -> path, dropping counts)
    return {key: path for key, (count, path) in tier_reports_with_counts.items()}
