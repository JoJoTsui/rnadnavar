"""Lossless, deterministic per-sample caller evidence for sampleless VCFs.

CALLER_EVIDENCE is URL-escaped JSON, one VCF String, version paired-v2.
Entries are observations, never independent support votes. Conflicting or
partial measurements remain separate; equal measurements share provenance.
"""

import copy
import json
import math
import re
from urllib.parse import quote, unquote

VERSION = "paired-v2"
KEY_FIELDS = ("modality", "caller", "sample_role", "sample_id", "allele", "alignment_round")
LEGACY_FIELDS = ("GT", "DP", "AD", "VAF", "VAF_SOURCE", "ALT_COUNT", "GQ")


def _plain(value):
    """Keep typed values and invalid measurements representable in strict JSON."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if hasattr(value, "tolist"):
        return _plain(value.tolist())
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f"Unsupported evidence value type: {type(value).__name__}")


def _json(value):
    return json.dumps(_plain(value), sort_keys=True, separators=(",", ":"), allow_nan=False)


def read_info(variant):
    """Read cyvcf2 INFO without lossy stringification of arrays and flags."""
    return {key: _plain(value) for key, value in variant.INFO}


def merge_evidence(*groups):
    """Union observations; coalesce exact measurements and accumulate sources."""
    merged = {}
    for entries in groups:
        for original in entries or []:
            entry = copy.deepcopy(original)
            for key in KEY_FIELDS:
                entry[key] = entry.get(key) or "unknown"
            entry["measurements"] = _plain(entry.get("measurements") or {})
            entry["available"] = bool(entry.get("available", bool(entry["measurements"])))
            sources = entry.pop("sources", [])
            if isinstance(sources, str):
                sources = [sources]
            signature = _json(entry)
            if signature not in merged:
                merged[signature] = (entry, set())
            merged[signature][1].update(str(s) for s in sources)
    return [dict(entry, sources=sorted(sources)) for _, (entry, sources) in sorted(merged.items())]


def encode_evidence(entries):
    return quote(_json({"version": VERSION, "entries": merge_evidence(entries)}), safe="")


def bind_evidence(entries, modality="unknown", alignment_round="unknown"):
    """Fill unknown identity from explicit context without rewriting origins.

    Round context can be a mapping such as {'DNA': 'first', 'RNA': 'realignment'}
    when a writer contains both modalities.
    """
    result = copy.deepcopy(entries)
    for entry in result:
        if entry.get("modality", "unknown") == "unknown":
            inferred_modality, _ = _identity(entry["caller"], modality)
            entry["modality"] = inferred_modality
        if entry.get("alignment_round", "unknown") == "unknown":
            entry["alignment_round"] = (
                alignment_round.get(entry.get("modality"), "unknown")
                if isinstance(alignment_round, dict) else alignment_round
            ) or "unknown"
    return merge_evidence(result)


def decode_evidence(value):
    """Reject corrupt/unknown schemas instead of silently dropping evidence."""
    if value is None or value == "" or value == ".":
        return []
    payload = json.loads(unquote(value))
    if not isinstance(payload, dict) or payload.get("version") != VERSION:
        raise ValueError("Unsupported CALLER_EVIDENCE version")
    entries = payload.get("entries")
    if not isinstance(entries, list):
        raise ValueError("CALLER_EVIDENCE entries must be a list")
    for entry in entries:
        if not isinstance(entry, dict) or not all(key in entry for key in KEY_FIELDS):
            raise ValueError("CALLER_EVIDENCE entry is missing identity fields")
        if not isinstance(entry.get("measurements"), dict):
            raise ValueError("CALLER_EVIDENCE measurements must be an object")
    return merge_evidence(entries)


def _identity(label, modality):
    match = re.match(r"^(DNA|RNA)[_:](.+)$", label, flags=re.I)
    if match:
        return match[1].upper(), match[2].lower()
    return modality or "unknown", label.lower()


def _legacy_pairs(value):
    # A pipe inside a phased genotype is data, not a caller separator.
    if isinstance(value, (tuple, list)):
        value = ",".join("." if item is None else str(item) for item in value)
    if value is None:
        return []
    result = []
    for part in re.split(r"\|(?=[A-Za-z][A-Za-z0-9_.-]*:)", str(value)):
        if ":" in part:
            result.append(part.split(":", 1))
    return result


def _legacy_matches(field, canonical, legacy):
    if legacy is None:
        return True
    if field in {"DP", "GQ", "ALT_COUNT", "VAF"}:
        try:
            left, right = float(canonical), float(legacy)
        except (TypeError, ValueError, OverflowError):
            return False
        if not math.isfinite(left) or not math.isfinite(right):
            return False
        # Legacy VAF is printed to four decimal places by the VCF writer.
        return abs(left - right) <= 0.000050000001 if field == "VAF" else left == right
    return str(canonical) == str(legacy)


def evidence_from_record(info, caller, allele, genotype=None, normal_genotype=None,
                         modality="unknown", alignment_round="unknown",
                         sample_id="unknown", normal_sample_id="unknown", source=None):
    """Read canonical evidence, enrich from legacy fields, or capture raw calls.

    Unknown identity is bound only to explicit enclosing context. Existing known
    rounds/modalities are never relabeled on rescue/realignment ingestion.
    ``allele`` is the normalized chrom:pos:ref:alt variant key. Sampleless legacy
    measurements have unknown sample identity, not an invented sample name.
    """
    modality = modality or "unknown"
    if modality == "unknown":
        inferred = re.match(r"^(DNA|RNA)[_:]", caller, flags=re.I)
        if inferred:
            modality = inferred[1].upper()
    alignment_round = alignment_round or "unknown"
    source = source or caller
    entries = decode_evidence(info.get("CALLER_EVIDENCE"))
    for entry in entries:
        for key, context in (("modality", modality), ("alignment_round", alignment_round)):
            if entry[key] == "unknown":
                entry[key] = context

    legacy = {}
    ambiguous = {}
    for role, prefix in (("tumor", ""), ("normal", "NORMAL_")):
        for field in LEGACY_FIELDS:
            info_key = f"{prefix}{field}_BY_CALLER"
            raw = info.get(info_key)
            if raw is not None and "||" in str(raw):
                # Old rescue concatenation did not define boundaries reliably.
                # Keep the complete source text without fabricating a caller GT.
                ambiguous[info_key] = _plain(raw)
                continue
            for label, value in _legacy_pairs(raw):
                entry_modality, entry_caller = _identity(label, modality)
                key = (entry_modality, entry_caller, role)
                legacy.setdefault(key, {})[field] = value if value not in ("", ".") else None
    for (entry_modality, entry_caller, role), measurements in legacy.items():
        # Canonical evidence wins only if it already contains this legacy
        # information. Partial or conflicting source measurements are retained.
        compatible = any(
            e["modality"] == entry_modality and e["caller"] == entry_caller
            and e["sample_role"] == role and e["allele"] == allele
            and (alignment_round == "unknown" or e["alignment_round"] == alignment_round)
            and all(_legacy_matches(k, e["measurements"].get(k), v)
                    for k, v in measurements.items()) for e in entries
        )
        if not compatible:
            entries.append(dict(modality=entry_modality, caller=entry_caller,
                                sample_role=role, sample_id="unknown", allele=allele,
                                alignment_round=alignment_round, measurements=measurements,
                                available=any(v is not None for v in measurements.values()),
                                sources=[f"{source}:legacy_INFO"]))
    if ambiguous:
        entries.append(dict(modality=modality, caller="unknown", sample_role="unknown",
                            sample_id="unknown", allele=allele, alignment_round=alignment_round,
                            measurements={"LEGACY_RAW_FIELDS": ambiguous}, available=False,
                            sources=[f"{source}:ambiguous_legacy_INFO"]))
    if "consensus" not in caller.lower() and "rescue" not in caller.lower():
        entry_modality, entry_caller = _identity(caller, modality)
        for role, sample, measurements in (("tumor", sample_id, genotype),
                                           ("normal", normal_sample_id, normal_genotype)):
            entries.append(dict(modality=entry_modality, caller=entry_caller,
                                sample_role=role, sample_id=sample or "unknown", allele=allele,
                                alignment_round=alignment_round,
                                measurements=_plain(measurements or {}),
                                available=measurements is not None, sources=[source]))
    else:
        # Explicit absence of normal measurements accompanies each legacy tumor
        # observation. This is a missing marker, never support evidence.
        for entry in list(entries):
            if entry["sample_role"] != "tumor":
                continue
            if not any(e["sample_role"] == "normal" and all(e[k] == entry[k] for k in
                       ("modality", "caller", "allele", "alignment_round")) for e in entries):
                missing = dict(entry, sample_role="normal", sample_id="unknown",
                               measurements={}, available=False)
                entries.append(missing)
    return merge_evidence(entries)
