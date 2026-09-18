"""Separate native Germline/Reference nominations from Somatic admission.

DeepSomatic's somatic writer recodes GT: GERMLINE can have GT=0/0.
Use explicit native FILTER and model confidence, not ordinary diploid GT.
Mutect2/Strelka rejected Somatic calls never nominate a negative class alone.
Thresholds are provisional; candidate labels are not biological training approval.
"""
from .refined_native_policy import counts, number

POLICY = "native_three_class_v1"
CALLERS = {"mutect2", "strelka", "deepsomatic"}

def read_counts(evidence, role):
    ad = counts(evidence.get(role + "_AD"), 2)
    dp = number(evidence.get(role + "_DP"))
    if (ad is None or dp is None or not dp.is_integer() or dp <= 0
            or sum(ad) <= 0 or sum(ad) > dp):
        return None
    return ad, dp

def native_negative(data):
    """Nominate a class from explicit DeepSomatic model output, never failed PASS."""
    ref, alt = data.get("REF", ""), data.get("ALT", "")
    if (not ref or not alt or set(ref + alt) - set("ACGT") or "," in alt
            or data.get("is_multiallelic")):
        return None, "unsupported_allele"
    filters = dict(zip(data.get("callers", []), data.get("filters_original", [])))
    label = str(filters.get("deepsomatic", "")).upper()
    target = {"GERMLINE":"Germline", "REFCALL":"Reference"}.get(label)
    if not target:
        return None, "no_native_negative_nomination"
    ev = data.get("native_evidence", {}).get("deepsomatic", {})
    quality = number(ev.get("tumor_GQ"))
    measured = read_counts(ev, "tumor")
    if quality is None or quality < 30 or measured is None:
        return None, "insufficient_native_class_confidence"
    ad, dp = measured
    if dp < 20 or sum(ad) < 20:
        return None, "insufficient_native_class_depth"
    # PL describes model class probabilities before GT recoding. If supplied,
    # require the biallelic class ordering to agree; never reinterpret GT=0/0.
    if ev.get("tumor_PL") is not None:
        pl = counts(ev["tumor_PL"], 3)
        index = 1 if target == "Germline" else 0
        if pl is None or any(pl[j] <= pl[index] for j in range(3) if j != index):
            return None, "inconsistent_native_class_likelihoods"
    if target == "Germline" and ad[1] < 3:
        return None, "insufficient_native_germline_alt"
    if target == "Reference" and (ad[1] > 2 or ad[1] / sum(ad) > 0.05):
        return None, "conflicting_native_reference_alt"
    return target, "deepsomatic_native_" + target.lower()

def evaluate(data):
    label, reason = native_negative(data)
    if label is None:
        return None, reason
    corroborators = []
    for caller, ev in data.get("native_evidence", {}).items():
        if caller not in CALLERS or caller == "deepsomatic":
            continue
        normal = read_counts(ev, "normal")
        tumor = read_counts(ev, "tumor")
        if label == "Germline" and normal:
            ad, dp = normal
            if sum(ad) >= 60 and ad[1] == 0:
                return None, "conflicting_normal_reference_evidence:" + caller
            if sum(ad) >= 20 and ad[1] >= 5 and ad[1] / sum(ad) >= 0.2:
                corroborators.append(caller)
        if label == "Reference":
            for measured in (normal, tumor):
                if measured:
                    ad, dp = measured
                    if sum(ad) >= 10 and ad[1] >= 3 and ad[1] / sum(ad) > 0.05:
                        return None, "conflicting_reference_alt_evidence:" + caller
            if normal and tumor and all(sum(ad) >= 60 and ad[1] == 0 for ad, dp in (normal, tumor)):
                corroborators.append(caller)
    suffix = "corroborated:" + "+".join(sorted(corroborators)) if corroborators else "native_only"
    return label, reason + ":" + suffix

def classify(data, somatic_allowed):
    negative, reason = evaluate(data)
    if reason.startswith("conflicting_"):
        return "NoConsensus", reason
    if somatic_allowed and negative:
        return "NoConsensus", "somatic_negative_conflict:" + negative
    if somatic_allowed:
        return "Somatic", "refined_somatic_admission"
    if negative:
        return negative, reason
    return "NoConsensus", reason
