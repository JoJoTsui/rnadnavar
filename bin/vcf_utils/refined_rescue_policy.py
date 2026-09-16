"""Experimental rescue gate primitives; not enabled in the workflow.

DNA nomination is observed positive tumor alternate evidence, NOT a Somatic
vote. RNA votes must have unique caller identities within one alignment round.
Biological vetoes apply to additions, not to the retained DNA baseline.
"""
import math

CALLERS = frozenset({"deepsomatic", "mutect2", "strelka"})


def integer(value):
    try:
        n = float(value)
        return int(n) if math.isfinite(n) and n >= 0 and n.is_integer() else None
    except (TypeError, ValueError, OverflowError):
        return None


def nominates(caller, original_filter, tumor_alt):
    count = integer(tumor_alt)
    return (caller in CALLERS and count is not None and count > 0
            and (caller != "deepsomatic" or original_filter in {"PASS", "."}))


def rna_supports(caller, original_filter, tumor_alt):
    """A native PASS is insufficient: preserve the eligible ALT-read floor."""
    count = integer(tumor_alt)
    return (caller in CALLERS and original_filter in {"PASS", ".", "Somatic"}
            and count is not None and count >= 3)


def biological_veto(info):
    """Return a reason or None; malformed available AF fails closed."""
    for key, value in info.items():
        if key.casefold() != "gnomad_af" or value in (None, "", "."):
            continue
        try:
            values = value if isinstance(value, (list, tuple)) else str(value).split(",")
            afs = [float(v) for v in values]
            if not afs or any(not math.isfinite(v) or not 0 <= v <= 1 for v in afs):
                return "invalid_population_af"
            if max(afs) > 0.001:
                return "common_population_af"
        except (ValueError, TypeError, OverflowError):
            return "invalid_population_af"
    if (info.get("REDI_ACCESSION") not in (None, "", ".")
            and info.get("REDI_CANONICAL") == "YES"
            and integer(info.get("N_DNA_CALLERS_SOMATIC")) == 0):
        return "canonical_editing_without_dna_somatic"
    return None


def decide(ref, alt, candidate_filter, dna_nominators, rna_pass_callers, info):
    """Return (admitted, reason) for a proposed addition, never the baseline.

No record-count or first/realignment-round summation is allowed. The caller
    sets are eligible evidence (rna_supports) supplied by the adapter, not
    inferred from file names or raw PASS presence alone.
"""
    if len(ref) != 1 or len(alt) != 1 or set(ref + alt) - set("ACGT"):
        return False, "not_biallelic_snp"
    if candidate_filter not in {"Somatic", "PASS", "."}:
        return False, "candidate_not_somatic"
    if not (set(dna_nominators) & CALLERS):
        return False, "no_dna_nomination"
    if len(set(rna_pass_callers) & CALLERS) < 2:
        return False, "insufficient_distinct_rna_callers"
    veto = biological_veto(info)
    return (False, veto) if veto else (True, "dna_nominated_rna_supported")
