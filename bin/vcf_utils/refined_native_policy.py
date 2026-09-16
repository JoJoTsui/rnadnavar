"""Opt-in SEQC2 development policy; no truth, coordinates or sample names.

Inputs must be normalized biallelic records with caller-native measurements.
Missing measurements never become zero; caller read counts are not summed.
This module nominates candidates, leaving biological vetoes to classification.
"""
import math

POLICY = "seqc2_refined_v2"
SOFT = {"PASS", ".", "contamination", "weak_evidence"}


def number(value):
    if isinstance(value, (tuple, list)):
        if len(value) != 1:
            return None
        value = value[0]
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError, OverflowError):
        return None


def counts(value, length):
    if isinstance(value, str):
        value = value.split(",")
    if not isinstance(value, (list, tuple)) or len(value) != length:
        return None
    result = [number(x) for x in value]
    if any(x is None or x < 0 or not x.is_integer() for x in result):
        return None
    return result


def evaluate(data):
    """Return a branch name, or None; SNP rejection must not majority-fallback."""
    filters = dict(zip(data.get("callers", []), data.get("filters_original", [])))
    ds_filter, m2_filter = filters.get("deepsomatic"), filters.get("mutect2")
    soft = bool(m2_filter) and set(str(m2_filter).split(";")) <= SOFT
    evidence = data.get("native_evidence", {})
    ds, m2 = evidence.get("deepsomatic", {}), evidence.get("mutect2", {})
    quality = number(data.get("qualities_by_caller", {}).get("deepsomatic"))
    tlod, germq = number(m2.get("TLOD")), number(m2.get("GERMQ"))
    if data.get("is_multiallelic") or "," in str(data.get("ALT", "")):
        return None
    if data.get("is_snv"):
        if ds_filter in {"PASS", "."}:
            return "snv_native_pass"
        if (soft and quality is not None and quality > 0 and tlod is not None
                and tlod >= 12 and germq is not None and germq >= 60):
            return "snv_soft_corroborated"
        return None
    ref, alt = str(data.get("REF", "")), str(data.get("ALT", ""))
    if not ref or not alt or len(ref) == len(alt) or set(ref + alt) - set("ACGTN"):
        return None
    dt = counts(ds.get("tumor_AD"), 2)
    mt = counts(m2.get("tumor_AD"), 2)
    mn = counts(m2.get("normal_AD"), 2)
    nd = number(m2.get("normal_DP"))
    if (dt is None or mt is None or mn is None or nd is None or not nd.is_integer()
            or nd < 10 or mn[1] != 0 or dt[1] < 3
            or germq is None or germq < 20 or tlod is None or tlod <= 0
            or quality is None):
        return None
    if ds_filter in {"PASS", "."} and soft and mt[1] >= 2:
        if quality >= 20:
            return "indel_high_confidence"
        sb = counts(m2.get("tumor_SB"), 4)
        if quality >= 10 and sb is not None and min(sb[2:]) >= 1:
            return "indel_moderate_strand"
    if (m2_filter in {"PASS", "."} and mt[1] >= 3 and number(m2.get("ECNT")) == 1
            and ds_filter in {"PASS", ".", "RefCall"} and quality > 0):
        return "indel_reciprocal_single_event"
    return None


def trace(data, branch):
    """Compact INFO-safe measurements for the policy decision provenance."""
    parts = [f"rule:{POLICY}", f"branch:{branch}"]
    filters = dict(zip(data.get("callers", []), data.get("filters_original", [])))
    for caller in ("deepsomatic", "mutect2"):
        filt = str(filters.get(caller, ".")).replace(";", "+")
        parts.append(f"{caller}_filter:{filt}")
        qual = number(data.get("qualities_by_caller", {}).get(caller))
        parts.append(f"{caller}_qual:{qual if qual is not None else '.'}")
        for key, value in sorted(data.get("native_evidence", {}).get(caller, {}).items()):
            values = value if isinstance(value, (list, tuple)) else [value]
            encoded = "/".join(str(number(v)) if number(v) is not None else "." for v in values)
            parts.append(f"{caller}_{key}:{encoded}")
    return "|".join(parts)
