"""
lib/common.py — shared constants, disease logic, sample classification, partitioning.
Single source of truth used by both parse_projects_to_json.py and run_batch_from_json.py.
"""

import re
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Disease classification
# ---------------------------------------------------------------------------

# set1: ONLY samples whose disease normalizes to exactly "colorectal cancer"
# Colon cancer and rectal cancer are NOT colorectal — they go into sets 2-4.
COLORECTAL_EXACT = {"colorectal cancer", "colorectal"}

def normalize_disease(raw: str) -> str:
    """Lowercase and strip surrounding whitespace. Returns 'unknown' if empty."""
    result = raw.strip().lower()
    if not result:
        return "unknown"
    return result

def is_colorectal(disease: str) -> bool:
    """True only for 'Colorectal cancer' (exact, case-insensitive). NOT colon/rectal."""
    return normalize_disease(disease) in COLORECTAL_EXACT

# ---------------------------------------------------------------------------
# Sample status classification
# ---------------------------------------------------------------------------

REQUIRED_MODALITIES = ["DN", "DT", "RT"]
MOD_STATUS_CODE = {"DN": 0, "DT": 1, "RT": 2}   # nextflow CSV status column

def classify_sample(modalities: dict) -> tuple[str, str]:
    """
    Returns (status, reason).
    - incomplete : any required modality missing or has 0 valid pairs
    - extra      : all required present, but any has >1 pair
    - standard   : all required present, each with exactly 1 pair
    """
    missing = [m for m in REQUIRED_MODALITIES
               if m not in modalities or modalities[m]["n_pairs"] == 0]
    if missing:
        return "incomplete", f"missing or empty modalities: {missing}"

    extra_mods = [m for m in REQUIRED_MODALITIES if modalities[m]["n_pairs"] > 1]
    if extra_mods:
        return "extra", f"modalities with >1 pair: {extra_mods}"

    return "standard", "all required modalities have exactly 1 pair"

def is_eligible(sample: dict) -> bool:
    return sample["status"] in ("standard", "extra")

# ---------------------------------------------------------------------------
# FASTQ pair grouping
# ---------------------------------------------------------------------------

def pairs_from_paths(paths: list) -> list:
    """Group _1/_2 fastq paths into [{"r1": ..., "r2": ...}] sorted by filename."""
    r1 = sorted(p for p in paths if re.search(r'_1\.fastq', p))
    r2 = sorted(p for p in paths if re.search(r'_2\.fastq', p))
    return [{"r1": a, "r2": b} for a, b in zip(r1, r2)]

# ---------------------------------------------------------------------------
# Sample key (unique across projects)
# ---------------------------------------------------------------------------

def sample_key(project_id: str, patient_id: str) -> str:
    return f"{project_id}::{patient_id}"

def skey(s: dict) -> str:
    return sample_key(s["_project_id"], s["patient_id"])

# ---------------------------------------------------------------------------
# Partitioner
# ---------------------------------------------------------------------------

def partition_eligible(all_samples: list) -> dict:
    """
    Partition eligible samples into 4 sets:
      set1 : strict colorectal cancer only (disease == 'colorectal cancer')
      set2-4: all other eligible diseases, disease-exclusive across sets,
              near-equal size as secondary objective.

    Colon cancer and rectal cancer are NOT colorectal — they go into sets 2-4.

    Returns:
      {
        "assignment": {skey -> set_number},
        "exceptions": [...],          # cases where exclusivity was violated
        "bin_diseases": {2: [...], 3: [...], 4: [...]},
      }
    """
    eligible = [s for s in all_samples if is_eligible(s)]

    colorectal_samples = [s for s in eligible if is_colorectal(s["disease"])]
    other_samples      = [s for s in eligible if not is_colorectal(s["disease"])]

    assignment: dict = {}

    # set1 — strict colorectal only
    for s in colorectal_samples:
        assignment[skey(s)] = 1

    # group remaining by normalized disease name
    disease_groups: dict = defaultdict(list)
    for s in other_samples:
        disease_groups[normalize_disease(s["disease"])].append(s)

    # greedy bin-packing: largest disease group first, disease-exclusive per bin
    sorted_diseases = sorted(disease_groups.items(), key=lambda x: -len(x[1]))

    bins        = {2: [], 3: [], 4: []}
    bin_diseases: dict = {2: set(), 3: set(), 4: set()}
    exceptions  = []

    for disease, samples in sorted_diseases:
        # prefer a bin that doesn't already have this disease
        candidates = [b for b in [2, 3, 4] if disease not in bin_diseases[b]]
        if candidates:
            target = min(candidates, key=lambda b: len(bins[b]))
        else:
            # fallback: smallest bin, record exception
            target = min([2, 3, 4], key=lambda b: len(bins[b]))
            exceptions.append({
                "disease": disease,
                "bin": target,
                "reason": "disease already present in all bins; placed in smallest",
            })
        bins[target].extend(samples)
        bin_diseases[target].add(disease)
        for s in samples:
            assignment[skey(s)] = target

    return {
        "assignment": assignment,
        "exceptions": exceptions,
        "bin_diseases": {k: sorted(v) for k, v in bin_diseases.items()},
    }
