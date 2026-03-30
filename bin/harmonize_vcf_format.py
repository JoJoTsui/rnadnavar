#!/usr/bin/env python3
"""
VCF FORMAT Field Harmonizer

Normalizes FORMAT fields across DNA variant callers to Mutect2 conventions:
  - GT, AD, AF, DP, GQ as canonical fields
  - Caller-specific remapping rules applied per caller type
  - All non-target FORMAT fields preserved unchanged

Supported callers:
  - mutect2:      passthrough (fields already canonical)
  - strelka (SNV):  AU/CU/GU/TU base tier-0 counts → AD, AF, DP
  - strelka (Indel): TAR/TIR → AD, AF, DP
  - sage:         VF → AF; AD and DP already present
  - deepsomatic:  VAF → AF; AD and DP already present
"""

import argparse
import subprocess
import sys

from cyvcf2 import VCF, Writer


# ---------------------------------------------------------------------------
# Canonical FORMAT header definitions (Mutect2-compatible)
# ---------------------------------------------------------------------------

CANONICAL_FORMATS = {
    "AD": {
        "ID": "AD",
        "Number": "R",
        "Type": "Integer",
        "Description": "Allelic depths for the ref and alt alleles in the order listed",
    },
    "AF": {
        "ID": "AF",
        "Number": "A",
        "Type": "Float",
        "Description": "Allele fractions of alternate alleles in the tumor",
    },
    "DP": {
        "ID": "DP",
        "Number": "1",
        "Type": "Integer",
        "Description": "Approximate read depth",
    },
    "GQ": {
        "ID": "GQ",
        "Number": "1",
        "Type": "Integer",
        "Description": "Genotype Quality",
    },
}

TARGET_FIELDS = set(CANONICAL_FORMATS.keys())

# Base nucleotides for Strelka2 SNV tier-0 fields
STRELKA_SNV_FIELDS = {"A": "AU", "C": "CU", "G": "GU", "T": "TU"}
STRELKA_SNV_FIELD_NAMES = set(STRELKA_SNV_FIELDS.values())  # AU, CU, GU, TU
STRELKA_INDEL_FIELD_NAMES = {"TAR", "TIR"}


# ---------------------------------------------------------------------------
# Caller detection helpers
# ---------------------------------------------------------------------------


def detect_strelka_variant_type(vcf):
    """Return 'snv', 'indel', or None based on FORMAT fields in VCF header."""
    format_ids = {h["ID"] for h in vcf.header_iter() if h["HeaderType"] == "FORMAT"}
    if STRELKA_SNV_FIELD_NAMES & format_ids:
        return "snv"
    if STRELKA_INDEL_FIELD_NAMES & format_ids:
        return "indel"
    return None


# ---------------------------------------------------------------------------
# Per-sample FORMAT remapping functions
# ---------------------------------------------------------------------------


def _safe_int(val):
    """Return int or None if val is None/missing."""
    if val is None:
        return None
    try:
        return int(val)
    except (TypeError, ValueError):
        return None


def _safe_float(val):
    """Return float or None if val is None/missing."""
    if val is None:
        return None
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def remap_mutect2(variant, sample_idx):
    """Mutect2: passthrough — fields already canonical."""
    ad = variant.format("AD")
    af = variant.format("AF")
    dp = variant.format("DP")
    gq = variant.format("GQ")

    ad_val = list(ad[sample_idx]) if ad is not None else None
    af_val = list(af[sample_idx]) if af is not None else None
    dp_val = _safe_int(dp[sample_idx][0]) if dp is not None else None
    gq_val = _safe_int(gq[sample_idx][0]) if gq is not None else None

    return ad_val, af_val, dp_val, gq_val


def remap_strelka_snv(variant, sample_idx, ref_base, alt_base):
    """
    Strelka2 SNV: compute AD/AF/DP from base tier-0 counts.
    AU[0], CU[0], GU[0], TU[0] are the first-tier counts.
    """
    base_counts = {}
    for base, field in STRELKA_SNV_FIELDS.items():
        fmt = variant.format(field)
        if fmt is not None:
            # Each sample entry is a tuple (tier1, tier2); take tier1 (index 0)
            val = fmt[sample_idx]
            base_counts[base] = _safe_int(val[0]) if val is not None else None
        else:
            base_counts[base] = None

    ref_count = base_counts.get(ref_base.upper())
    alt_count = base_counts.get(alt_base.upper())

    # DP = sum of all base tier-0 counts
    all_counts = [v for v in base_counts.values() if v is not None]
    dp_val = sum(all_counts) if all_counts else None

    if ref_count is not None and alt_count is not None:
        ad_val = [ref_count, alt_count]
        total = ref_count + alt_count
        af_val = [alt_count / total] if total > 0 else [0.0]
    else:
        ad_val = None
        af_val = None

    # GQ passthrough
    gq_fmt = variant.format("GQ")
    gq_val = _safe_int(gq_fmt[sample_idx][0]) if gq_fmt is not None else None

    return ad_val, af_val, dp_val, gq_val


def remap_strelka_indel(variant, sample_idx):
    """
    Strelka2 Indel: compute AD/AF/DP from TAR/TIR.
    TAR[0] = ref tier-1 count, TIR[0] = alt tier-1 count.
    """
    tar_fmt = variant.format("TAR")
    tir_fmt = variant.format("TIR")

    tar = _safe_int(tar_fmt[sample_idx][0]) if tar_fmt is not None else None
    tir = _safe_int(tir_fmt[sample_idx][0]) if tir_fmt is not None else None

    if tar is not None and tir is not None:
        ad_val = [tar, tir]
        total = tar + tir
        af_val = [tir / total] if total > 0 else [0.0]
        dp_val = total
    else:
        ad_val = None
        af_val = None
        dp_val = None

    gq_fmt = variant.format("GQ")
    gq_val = _safe_int(gq_fmt[sample_idx][0]) if gq_fmt is not None else None

    return ad_val, af_val, dp_val, gq_val


def remap_sage(variant, sample_idx):
    """SAGE: remap VF → AF; AD and DP already present."""
    ad_fmt = variant.format("AD")
    dp_fmt = variant.format("DP")
    gq_fmt = variant.format("GQ")

    ad_val = list(ad_fmt[sample_idx]) if ad_fmt is not None else None
    dp_val = _safe_int(dp_fmt[sample_idx][0]) if dp_fmt is not None else None
    gq_val = _safe_int(gq_fmt[sample_idx][0]) if gq_fmt is not None else None

    # Use AF if present, otherwise fall back to VF
    af_fmt = variant.format("AF")
    if af_fmt is not None:
        af_val = [_safe_float(af_fmt[sample_idx][0])]
    else:
        vf_fmt = variant.format("VF")
        if vf_fmt is not None:
            af_val = [_safe_float(vf_fmt[sample_idx][0])]
        else:
            af_val = None

    return ad_val, af_val, dp_val, gq_val


def remap_deepsomatic(variant, sample_idx):
    """DeepSomatic: remap VAF → AF; AD and DP already present."""
    ad_fmt = variant.format("AD")
    dp_fmt = variant.format("DP")
    gq_fmt = variant.format("GQ")

    ad_val = list(ad_fmt[sample_idx]) if ad_fmt is not None else None
    dp_val = _safe_int(dp_fmt[sample_idx][0]) if dp_fmt is not None else None
    gq_val = _safe_int(gq_fmt[sample_idx][0]) if gq_fmt is not None else None

    # Use AF if present, otherwise fall back to VAF
    af_fmt = variant.format("AF")
    if af_fmt is not None:
        af_val = [_safe_float(af_fmt[sample_idx][0])]
    else:
        vaf_fmt = variant.format("VAF")
        if vaf_fmt is not None:
            af_val = [_safe_float(vaf_fmt[sample_idx][0])]
        else:
            af_val = None

    return ad_val, af_val, dp_val, gq_val


# ---------------------------------------------------------------------------
# Header rewriting
# ---------------------------------------------------------------------------


def rewrite_header(vcf):
    """
    Add/replace canonical FORMAT header lines in the VCF header.
    Returns a new VCF header string with canonical fields declared.
    """
    for field_id, meta in CANONICAL_FORMATS.items():
        # Remove existing declaration if present, then add canonical one
        vcf.add_format_to_header(
            {
                "ID": meta["ID"],
                "Number": meta["Number"],
                "Type": meta["Type"],
                "Description": meta["Description"],
            }
        )
    return vcf


# ---------------------------------------------------------------------------
# Main harmonization logic
# ---------------------------------------------------------------------------


def _compute_canonical_for_sample(variant, sample_idx, caller_lower, strelka_type, ref, alt_base):
    """Compute (ad, af, dp, gq) for a single sample using caller-specific rules."""
    if "mutect2" in caller_lower:
        return remap_mutect2(variant, sample_idx)

    elif "strelka" in caller_lower:
        original_fmt_keys = variant.FORMAT
        if strelka_type == "snv":
            return remap_strelka_snv(variant, sample_idx, ref[0], alt_base[0])
        elif strelka_type == "indel":
            return remap_strelka_indel(variant, sample_idx)
        else:
            # Per-variant fallback detection
            fmt_ids = set(original_fmt_keys)
            if STRELKA_SNV_FIELD_NAMES & fmt_ids:
                return remap_strelka_snv(variant, sample_idx, ref[0], alt_base[0])
            elif STRELKA_INDEL_FIELD_NAMES & fmt_ids:
                return remap_strelka_indel(variant, sample_idx)
            else:
                return None, None, None, None

    elif "sage" in caller_lower:
        return remap_sage(variant, sample_idx)

    elif "deepsomatic" in caller_lower:
        return remap_deepsomatic(variant, sample_idx)

    else:
        # Unknown caller: passthrough existing fields
        ad_fmt = variant.format("AD")
        af_fmt = variant.format("AF")
        dp_fmt = variant.format("DP")
        gq_fmt = variant.format("GQ")
        ad = list(ad_fmt[sample_idx]) if ad_fmt is not None else None
        af = [_safe_float(af_fmt[sample_idx][0])] if af_fmt is not None else None
        dp = _safe_int(dp_fmt[sample_idx][0]) if dp_fmt is not None else None
        gq = _safe_int(gq_fmt[sample_idx][0]) if gq_fmt is not None else None
        return ad, af, dp, gq


def harmonize(vcf_path, caller, output_path):
    """
    Read input VCF, apply caller-specific FORMAT remapping, write harmonized VCF.
    Output is written uncompressed to a temp path, then bgzip-compressed.
    """
    import numpy as np

    vcf = VCF(vcf_path)

    # Detect strelka variant type before rewriting header
    strelka_type = None
    if "strelka" in caller.lower():
        strelka_type = detect_strelka_variant_type(vcf)
        # Re-open after header scan
        vcf = VCF(vcf_path)

    # Rewrite header with canonical FORMAT declarations
    rewrite_header(vcf)

    # Determine output path (write uncompressed, pipe through bgzip)
    if output_path.endswith(".gz"):
        uncompressed_path = output_path[:-3]  # strip .gz
    else:
        uncompressed_path = output_path

    writer = Writer(uncompressed_path, vcf)

    n_samples = len(vcf.samples)
    caller_lower = caller.lower()
    INT_MISSING = -2147483648  # cyvcf2 sentinel for missing int

    for variant in vcf:
        ref = variant.REF
        alts = variant.ALT
        alt_base = alts[0] if alts else "N"

        # Collect per-sample canonical values
        all_ad = []   # list of [ref_count, alt_count] per sample (or None)
        all_af = []   # list of [freq] per sample (or None)
        all_dp = []   # list of int per sample (or None)
        all_gq = []   # list of int per sample (or None)

        for sample_idx in range(n_samples):
            ad, af, dp, gq = _compute_canonical_for_sample(
                variant, sample_idx, caller_lower, strelka_type, ref, alt_base
            )
            all_ad.append(ad)
            all_af.append(af)
            all_dp.append(dp)
            all_gq.append(gq)

        # Set AD (Number=R → 2 values: ref, alt)
        if any(v is not None for v in all_ad):
            ad_arr = np.full((n_samples, 2), INT_MISSING, dtype=np.int32)
            for i, v in enumerate(all_ad):
                if v is not None and len(v) >= 2:
                    ad_arr[i] = [int(v[0]) if v[0] is not None else INT_MISSING,
                                 int(v[1]) if v[1] is not None else INT_MISSING]
            try:
                variant.set_format("AD", ad_arr)
            except Exception:
                pass

        # Set AF (Number=A → 1 value: alt freq)
        if any(v is not None for v in all_af):
            af_arr = np.full((n_samples, 1), float("nan"), dtype=np.float32)
            for i, v in enumerate(all_af):
                if v is not None and len(v) >= 1 and v[0] is not None:
                    af_arr[i, 0] = float(v[0])
            try:
                variant.set_format("AF", af_arr)
            except Exception:
                pass

        # Set DP (Number=1 → scalar per sample)
        if any(v is not None for v in all_dp):
            dp_arr = np.full(n_samples, INT_MISSING, dtype=np.int32)
            for i, v in enumerate(all_dp):
                if v is not None:
                    dp_arr[i] = int(v)
            try:
                variant.set_format("DP", dp_arr)
            except Exception:
                pass

        # Set GQ (Number=1 → scalar per sample)
        if any(v is not None for v in all_gq):
            gq_arr = np.full(n_samples, INT_MISSING, dtype=np.int32)
            for i, v in enumerate(all_gq):
                if v is not None:
                    gq_arr[i] = int(v)
            try:
                variant.set_format("GQ", gq_arr)
            except Exception:
                pass

        writer.write_record(variant)

    writer.close()
    vcf.close()

    # bgzip compress if needed
    if output_path.endswith(".gz"):
        _bgzip(uncompressed_path, output_path)


def _bgzip(input_path, output_path):
    """Compress input_path to output_path using bgzip."""
    try:
        with open(output_path, "wb") as out_fh:
            subprocess.run(
                ["bgzip", "-c", input_path],
                stdout=out_fh,
                stderr=subprocess.PIPE,
                check=True,
            )
    except subprocess.CalledProcessError as e:
        print(
            f"ERROR: bgzip failed: {e.stderr.decode()}", file=sys.stderr
        )
        sys.exit(1)
    except FileNotFoundError:
        print("ERROR: bgzip not found on PATH", file=sys.stderr)
        sys.exit(1)

    # Remove uncompressed temp file
    import os
    try:
        os.remove(input_path)
    except OSError:
        pass


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def argparser():
    parser = argparse.ArgumentParser(
        description="Harmonize VCF FORMAT fields to Mutect2 conventions across callers",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--vcf", required=True, help="Input VCF file (may be bgzip-compressed)")
    parser.add_argument(
        "--caller",
        required=True,
        help="Variant caller name (mutect2, strelka, sage, deepsomatic, ...)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output VCF path (use .vcf.gz for bgzip-compressed output)",
    )
    return parser.parse_args()


def main():
    args = argparser()
    print(f"Harmonizing VCF: {args.vcf}")
    print(f"  Caller: {args.caller}")
    print(f"  Output: {args.output}")
    harmonize(args.vcf, args.caller, args.output)
    print("Done.")


if __name__ == "__main__":
    main()
