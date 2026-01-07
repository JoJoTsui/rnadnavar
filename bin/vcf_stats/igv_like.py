#!/usr/bin/env python3
"""
IGV-like BAM Visualization for Variants

Generates static, IGV-like plots around variant positions for provided BAM/CRAM files.
Shows reference sequence once at the top, then ordered sample tracks (DNA_NORMAL, DNA_TUMOR, RNA_TUMOR).
Highlights mismatches against reference and exports PNG and HTML.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


class IGVLikePlotter:
    def __init__(self, bam_files: Dict[str, Path], ref_fasta: Path):
        self.bam_files = {k: Path(v) for k, v in bam_files.items()}
        self.ref_fasta = Path(ref_fasta)

        try:
            import matplotlib.pyplot as plt  # noqa: F401
            import pysam  # noqa: F401
        except ImportError as e:
            raise RuntimeError("matplotlib and pysam are required for IGV-like plotting") from e
        if not self.ref_fasta.exists():
            raise RuntimeError("Reference FASTA (indexed) is required for visualization")

    def plot_variant_region(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        flank: int = 60,
        samples_subset: Optional[List[str]] = None,
        out_png: Optional[Path] = None,
        out_html: Optional[Path] = None,
        title: Optional[str] = None,
    ) -> Tuple[Optional[Path], Optional[Path]]:
        import matplotlib.pyplot as plt
        import pysam

        start = max(1, pos - flank)
        end = pos + flank

        sample_names = list(self.bam_files.keys()) if samples_subset is None else [s for s in samples_subset if s in self.bam_files]
        if not sample_names:
            return None, None

        # Enforce desired order
        order = [s for s in ["DNA_NORMAL", "DNA_TUMOR", "RNA_TUMOR"] if s in sample_names]
        order += [s for s in sample_names if s not in order]
        sample_names = order

        n = len(sample_names)
        fig, axes = plt.subplots(n + 1, 1, figsize=(13, 3.0 * (n + 1)), sharex=True)
        # Ensure axes is a list: when n+1=1, axes is a scalar Axes; when n+1>1, axes is numpy array
        if hasattr(axes, '__len__'):
            axes = list(axes)  # numpy array or tuple -> list
        else:
            axes = [axes]  # scalar Axes -> list with one element

        # Reference track at top
        ref_ax = axes[0]
        ref_fa = pysam.FastaFile(str(self.ref_fasta))
        ref_seq = ref_fa.fetch(chrom, start - 1, end).upper()
        xs = list(range(start, end + 1))
        ref_ax.plot(xs, [0] * len(xs), alpha=0)
        for i, base in enumerate(ref_seq):
            x = start + i
            ref_ax.text(x, 1.0, base, fontsize=9, ha="center", va="top", color="#37474f")
        ref_ax.axvline(pos, color="#ef5350", linestyle="--", linewidth=1.0)
        ref_ax.set_ylim(0, 2)
        ref_ax.set_ylabel("REFERENCE", fontsize=10)
        ref_ax.grid(True, axis="y", linestyle=":", alpha=0.2)

        # Sample tracks
        for ax, sample in zip(axes[1:], sample_names):
            bam_path = str(self.bam_files[sample])
            try:
                bam = pysam.AlignmentFile(bam_path, "rb")
            except ValueError:
                bam = pysam.AlignmentFile(bam_path, "rc")

            # Draw alignments and highlight mismatches
            y = 1
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped:
                    continue
                coords = read.get_reference_positions(full_length=False)
                qpos = read.get_aligned_pairs(matches_only=False)
                if coords:
                    ax.hlines(y, min(coords) + 1, max(coords) + 1, color="#78909c", linewidth=2, alpha=0.9)
                for qidx, ref_idx in qpos:
                    if qidx is None or ref_idx is None:
                        continue
                    ref_pos = ref_idx + 1
                    if ref_pos < start or ref_pos > end:
                        continue
                    ref_base = ref_seq[ref_pos - start]
                    read_base = read.query_sequence[qidx].upper()
                    if read_base != ref_base:
                        ax.scatter(ref_pos, y, color="#d32f2f", s=16, marker="s")
                        ax.text(ref_pos, y + 0.15, read_base, fontsize=8, ha="center", color="#c62828", fontweight="bold")
                y += 1

            bam.close()

            ax.axvline(pos, color="#ef5350", linestyle="--", linewidth=1.0)
            ax.set_ylabel(sample, fontsize=9)
            ax.set_ylim(0, max(y, 3))
            ax.grid(True, axis="y", linestyle=":", alpha=0.25)

        axes[-1].set_xlabel(f"{chrom}:{start}-{end}", fontsize=10)
        ttl = title or f"{chrom}:{pos} {ref}>{alt}"
        fig.suptitle(ttl, fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.95])

        saved_png = None
        if out_png:
            out_png = Path(out_png)
            out_png.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_png, dpi=140)
            saved_png = out_png

        saved_html = None
        if out_html and saved_png:
            out_html = Path(out_html)
            out_html.parent.mkdir(parents=True, exist_ok=True)
            png_name = Path(saved_png).name
            html = """
<!doctype html>
<html>
  <head>
    <meta charset="utf-8" />
    <title>{title}</title>
    <style>
      body {{ font-family: system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 16px; }}
      h3 {{ margin: 8px 0 16px; }}
      .img-wrap {{ border: 1px solid #e0e0e0; box-shadow: 0 2px 6px rgba(0,0,0,0.08); padding: 8px; }}
    </style>
  </head>
  <body>
    <h3>{title}</h3>
    <div class="img-wrap">
      <img src="{png_name}" alt="{title}" style="max-width: 100%; height: auto;" />
    </div>
  </body>
</html>
""".format(title=ttl, png_name=png_name)
            out_html.write_text(html)

        import matplotlib.pyplot as plt
        plt.close(fig)
        return saved_png, saved_html


def visualize_variants_igv_like(
    variants: pd.DataFrame,
    bam_files: Dict[str, Path],
    output_dir: Path,
    ref_fasta: Path,
    flank: int = 60,
    samples_subset: Optional[List[str]] = None,
) -> List[Path]:
    plotter = IGVLikePlotter(bam_files, ref_fasta=ref_fasta)
    out_paths: List[Path] = []
    output_dir = Path(output_dir)

    for _, row in variants.iterrows():
        chrom, pos, ref, alt = row["chrom"], int(row["pos"]), row["ref"], row["alt"]
        cat = str(row.get("filter_category", "Unknown"))
        tier = str(row.get("tier", "T?"))
        vid = str(row.get("variant_id", f"{chrom}:{pos}:{ref}>{alt}"))

        out_png = output_dir / cat / tier / f"{vid}.png"
        out_html = output_dir / cat / tier / f"{vid}.html"
        title = f"{cat} {tier} | {chrom}:{pos} {ref}>{alt}"

        saved_png, saved_html = plotter.plot_variant_region(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            flank=flank,
            samples_subset=samples_subset,
            out_png=out_png,
            out_html=out_html,
            title=title,
        )
        if saved_png:
            out_paths.append(saved_png)

    return out_paths
