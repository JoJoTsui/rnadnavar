#!/usr/bin/env python3
"""
IGV-like BAM Visualization for Variants

Generates static, lightweight IGV-like plots around variant positions
for provided BAM/CRAM files. Uses pysam pileups and matplotlib.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


class IGVLikePlotter:
    def __init__(self, bam_files: Dict[str, Path], ref_fasta: Optional[Path] = None):
        """
        bam_files: mapping of sample_name -> BAM/CRAM path
        ref_fasta: optional reference FASTA (indexed) for context
        """
        self.bam_files = {k: Path(v) for k, v in bam_files.items()}
        self.ref_fasta = Path(ref_fasta) if ref_fasta else None

        try:
            import matplotlib.pyplot as plt  # noqa: F401
            import pysam  # noqa: F401
        except ImportError as e:
            raise RuntimeError("matplotlib and pysam are required for IGV-like plotting") from e

    @staticmethod
    def _base_counts_at(pileups) -> Dict[str, int]:
        counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "INS": 0, "DEL": 0}
        for pr in pileups:
            if pr.is_del:
                counts["DEL"] += 1
            elif pr.is_refskip:
                continue
            else:
                b = pr.alignment.query_sequence[pr.query_position]
                counts[b.upper() if b else "N"] = counts.get(b.upper() if b else "N", 0) + 1
        return counts

    def plot_variant_region(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        window: int = 100,
        samples_subset: Optional[List[str]] = None,
        out_path: Optional[Path] = None,
        title: Optional[str] = None,
    ) -> Optional[Path]:
        """
        Create a multi-panel plot with coverage and base counts at the variant site
        for each selected sample.
        """
        import matplotlib.pyplot as plt
        import pysam

        start = max(1, pos - window)
        end = pos + window

        sample_names = list(self.bam_files.keys()) if samples_subset is None else [s for s in samples_subset if s in self.bam_files]
        if not sample_names:
            return None

        n = len(sample_names)
        fig, axes = plt.subplots(n, 1, figsize=(10, 2.4 * n), sharex=True)
        if n == 1:
            axes = [axes]

        for ax, sample in zip(axes, sample_names):
            bam_path = str(self.bam_files[sample])
            try:
                bam = pysam.AlignmentFile(bam_path, "rb")
            except ValueError:
                # maybe CRAM
                bam = pysam.AlignmentFile(bam_path, "rc")

            # Coverage across window
            cov = [0] * (end - start + 1)
            for pileupcolumn in bam.pileup(chrom, start - 1, end, truncate=True, stepper="all"):
                idx = pileupcolumn.pos + 1 - start
                if 0 <= idx < len(cov):
                    cov[idx] = pileupcolumn.nsegments

            # Base counts exactly at variant site
            base_counts = {}
            for pileupcolumn in bam.pileup(chrom, pos - 1, pos, truncate=True, stepper="all"):
                if pileupcolumn.pos + 1 == pos:
                    base_counts = self._base_counts_at(pileupcolumn.pileups)
                    break

            bam.close()

            # Plot coverage
            xs = list(range(start, end + 1))
            ax.fill_between(xs, cov, step="mid", color="#cfd8dc")
            ax.plot(xs, cov, color="#90a4ae", linewidth=1)

            # Vertical line at variant pos
            ax.axvline(pos, color="#ef5350", linestyle="--", linewidth=1.2)

            # Annotate base counts at position
            if base_counts:
                lbl = ", ".join([f"{k}:{v}" for k, v in base_counts.items() if v > 0])
                ax.text(pos, max(cov) * 0.85 if cov else 1.0, lbl, fontsize=9, color="#37474f", ha="left", va="top")

            ax.set_ylabel(f"{sample}\ncoverage", fontsize=9)
            ax.grid(True, axis="y", linestyle=":", alpha=0.4)

        axes[-1].set_xlabel(f"{chrom}:{start}-{end}", fontsize=10)
        ttl = title or f"{chrom}:{pos} {ref}>{alt}"
        fig.suptitle(ttl, fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.96])

        if out_path:
            out_path = Path(out_path)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_path, dpi=140)
            plt.close(fig)
            return out_path

        plt.show()
        return None


def visualize_variants_igv_like(
    variants: pd.DataFrame,
    bam_files: Dict[str, Path],
    output_dir: Path,
    window: int = 100,
    samples_subset: Optional[List[str]] = None,
) -> List[Path]:
    """
    Convenience wrapper to render a list/DataFrame of variants.
    variants must contain columns: chrom, pos, ref, alt, filter_category, tier, variant_id
    """
    plotter = IGVLikePlotter(bam_files)
    out_paths: List[Path] = []
    output_dir = Path(output_dir)

    for _, row in variants.iterrows():
        chrom, pos, ref, alt = row["chrom"], int(row["pos"]), row["ref"], row["alt"]
        cat = str(row.get("filter_category", "Unknown"))
        tier = str(row.get("tier", "T?"))
        vid = str(row.get("variant_id", f"{chrom}:{pos}:{ref}>{alt}"))

        out_path = output_dir / cat / tier / f"{vid}.png"
        title = f"{cat} {tier} | {chrom}:{pos} {ref}>{alt}"

        saved = plotter.plot_variant_region(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            window=window,
            samples_subset=samples_subset,
            out_path=out_path,
            title=title,
        )
        if saved:
            out_paths.append(saved)

    return out_paths
