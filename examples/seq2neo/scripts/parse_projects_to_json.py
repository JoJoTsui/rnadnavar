#!/usr/bin/env python3
"""
parse_projects_to_json.py
─────────────────────────
Parse three project manifests into one unified JSON and produce per-set TSV lists.

Inputs (all under --seq2neo root):
  PRJNA298376.txt          tree-style manifest (├── box-drawing, \xa0 spaces)
  PRJNA298330.txt          pipe-style manifest (|-> markers)
  PRJNA298310.txt          pipe-style manifest (|-> markers), all melanoma
  PRJNA298330.disease.tsv  per-patient disease annotation from NCBI BioSample

Outputs:
  data/processed/merged.json          unified sample database
  data/processed/set{1-4}_samples.tsv partition membership lists

Sample status:
  standard   — exactly 1 DN + 1 DT + 1 RT pair
  extra      — all 3 modalities present, ≥1 has >1 pair  (first pair used at runtime)
  incomplete — missing DN, DT, or RT entirely

Partition sets:
  set1 — strict Colorectal cancer only
  set2-4 — all other eligible diseases, disease-exclusive across sets

Usage:
  python3 scripts/parse_projects_to_json.py
  python3 scripts/parse_projects_to_json.py --seq2neo /path/to/seq2neo
  python3 scripts/parse_projects_to_json.py --out /custom/path/merged.json
"""

import re
import sys
import json
import argparse
from pathlib import Path
from collections import defaultdict

# allow running from any directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.common import (
    normalize_disease, is_colorectal, classify_sample,
    pairs_from_paths, partition_eligible, skey,
    REQUIRED_MODALITIES,
)

SEQ2NEO_DEFAULT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Manifest parsers
# ---------------------------------------------------------------------------

def _parse_header(lines: list, project: dict):
    """Extract Study Title, Paper, Platform from first 10 lines."""
    for line in lines[:10]:
        if "Study Title--" in line:
            project["title"] = line.split("Study Title--", 1)[1].strip()
        elif "Paper--" in line:
            project["paper"] = line.split("Paper--", 1)[1].strip()
        elif "Sequencing platform--" in line:
            project["platform"] = line.split("Sequencing platform--", 1)[1].strip()


def _build_sample(project_id: str, patient_id: str, disease: str,
                  raw_paths: dict) -> dict:
    """Convert raw path lists → structured modality dict → classified sample."""
    modalities = {}
    for mod, paths in raw_paths.items():
        pairs = pairs_from_paths(paths)
        modalities[mod] = {"pairs": pairs, "n_pairs": len(pairs)}

    status, reason = classify_sample(modalities)
    return {
        "_project_id":       project_id,
        "patient_id":        patient_id,
        "disease":           disease,
        "disease_normalized": normalize_disease(disease),
        "modalities":        modalities,
        "status":            status,
        "status_reason":     reason,
    }


def parse_prjna298376(txt_path: Path) -> dict:
    """
    Tree-style format.
    Lines use ├── / └── / │ box-drawing chars and \xa0 non-breaking spaces.
    Modality headers appear as standalone tokens after stripping tree chars.
    Fastq paths are extracted by regex from anywhere in the line.
    """
    text = txt_path.read_text()
    project = {"project_id": "PRJNA298376", "title": "", "paper": "", "platform": "", "samples": []}
    _parse_header(text.splitlines(), project)

    for block in re.split(r'\nPatientID--', text)[1:]:
        lines       = block.splitlines()
        patient_id  = lines[0].strip()
        disease     = ""
        raw_paths: dict = defaultdict(list)
        current_mod = None

        for line in lines[1:]:
            if "README--" in line:
                disease = line.split("README--", 1)[1].strip()
                continue
            # strip all tree/box/whitespace chars; if result is a modality token → new section
            clean = re.sub(r'[├└│─\s\xa0]+', '', line)
            if re.match(r'^(DN|DT|RT|RN|T)$', clean):
                current_mod = "DT" if clean == "T" else clean   # "T" is a typo for "DT" in patient 4107
                continue
            # extract fastq path from anywhere in the line
            m = re.search(r'(/\S+\.fastq\.gz)', line)
            if m and current_mod:
                raw_paths[current_mod].append(m.group(1))

        project["samples"].append(
            _build_sample("PRJNA298376", patient_id, disease, raw_paths)
        )
    return project


def _parse_pipe_style(txt_path: Path, project_id: str,
                      disease_map: dict, default_disease: str) -> dict:
    """
    Shared parser for pipe-style manifests (PRJNA298330 and PRJNA298310).
    Format uses |->WES / |->RNA / |->Tumor / |->Normal section markers.
    Modality mapping:
      WES + Tumor  → DT
      WES + Normal → DN
      RNA + Tumor  → RT
      RNA + Normal → RN  (optional, not required for eligibility)
    """
    text = txt_path.read_text()
    project = {"project_id": project_id, "title": "", "paper": "", "platform": "", "samples": []}
    _parse_header(text.splitlines(), project)

    for block in re.split(r'\nPatientID--', text)[1:]:
        lines      = block.splitlines()
        patient_id = lines[0].strip()
        raw_paths: dict = defaultdict(list)
        section     = None   # "WES" | "RNA"
        current_mod = None

        for line in lines[1:]:
            stripped = line.strip()
            if   "|->WES"    in stripped: section = "WES"
            elif "|->RNA"    in stripped: section = "RNA"
            elif "|->Tumor"  in stripped:
                current_mod = {"WES": "DT", "RNA": "RT"}.get(section)
            elif "|->Normal" in stripped:
                current_mod = {"WES": "DN", "RNA": "RN"}.get(section)
            elif stripped.startswith("/") and ".fastq" in stripped:
                if current_mod:
                    raw_paths[current_mod].append(stripped)

        disease = disease_map.get(patient_id, default_disease)
        project["samples"].append(
            _build_sample(project_id, patient_id, disease, raw_paths)
        )
    return project


def parse_prjna298330(txt_path: Path, disease_tsv: Path) -> dict:
    """Pipe-style. Disease per patient from PRJNA298330.disease.tsv."""
    disease_map = {}
    if disease_tsv and disease_tsv.exists():
        for row in disease_tsv.read_text().splitlines():
            if row.startswith("#") or not row.strip():
                continue
            parts = row.split("\t")
            if len(parts) >= 4 and parts[0] != "patient_id":
                disease_map[parts[0]] = parts[3]
    return _parse_pipe_style(txt_path, "PRJNA298330", disease_map,
                             default_disease="Gastrointestinal cancer")


def parse_prjna298310(txt_path: Path) -> dict:
    """Pipe-style. All patients are Melanoma."""
    return _parse_pipe_style(txt_path, "PRJNA298310", disease_map={},
                             default_disease="Melanoma")


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _summarize(all_samples: list) -> dict:
    from collections import Counter
    status_counts = Counter(s["status"] for s in all_samples)
    set_counts    = Counter(s.get("partition_set") for s in all_samples
                            if s.get("partition_set") is not None)
    return {
        "total":     len(all_samples),
        "by_status": dict(status_counts),
        "by_set":    {f"set{k}": v for k, v in sorted(set_counts.items())},
    }


def _write_partition_tsv(all_samples: list, out_dir: Path):
    sets: dict = defaultdict(list)
    for s in all_samples:
        ps = s.get("partition_set")
        if ps is not None:
            sets[ps].append(s)
    for set_num, members in sorted(sets.items()):
        path = out_dir / f"set{set_num}_samples.tsv"
        rows = ["project_id\tpatient_id\tdisease\tstatus"]
        for s in members:
            rows.append(f"{s['_project_id']}\t{s['patient_id']}\t{s['disease']}\t{s['status']}")
        path.write_text("\n".join(rows) + "\n")
        print(f"  Written: {path}")


def _print_report(output: dict, all_samples: list):
    summ = output["summary"]
    part = output["partition"]

    print("\n=== VALIDATION REPORT ===")
    print(f"Total samples : {summ['total']}")

    print("\nBy status:")
    for k, v in summ["by_status"].items():
        print(f"  {k:12s}: {v}")

    print("\nBy partition set:")
    for k, v in summ["by_set"].items():
        print(f"  {k}: {v}")

    print("\nDisease assignment (set2-4):")
    for b, diseases in part["bin_diseases"].items():
        print(f"  set{b}: {', '.join(diseases) if diseases else '(empty)'}")

    if part["exceptions"]:
        print("\nPartition exceptions (exclusivity violated):")
        for ex in part["exceptions"]:
            print(f"  {ex}")
    else:
        print("\nNo partition exceptions — full disease exclusivity achieved.")

    print("\nIncomplete samples:")
    for s in all_samples:
        if s["status"] == "incomplete":
            print(f"  [{s['_project_id']}] patient {s['patient_id']:6s}"
                  f" | {s['disease']:45s} | {s['status_reason']}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Parse PRJNA manifests into unified JSON + partition TSVs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--seq2neo", default=str(SEQ2NEO_DEFAULT),
                    help="seq2neo root directory (default: parent of scripts/)")
    ap.add_argument("--out", default=None,
                    help="Output JSON path (default: <seq2neo>/data/processed/merged.json)")
    args = ap.parse_args()

    root     = Path(args.seq2neo)
    out_path = Path(args.out) if args.out else root / "data" / "processed" / "merged.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("Parsing manifests...")
    projects = [
        parse_prjna298376(root / "PRJNA298376.txt"),
        parse_prjna298330(root / "PRJNA298330.txt", root / "PRJNA298330.disease.tsv"),
        parse_prjna298310(root / "PRJNA298310.txt"),
    ]

    # flatten with project tag already set in _build_sample
    all_samples = [s for proj in projects for s in proj["samples"]]

    # partition
    part = partition_eligible(all_samples)
    for s in all_samples:
        s["partition_set"] = part["assignment"].get(skey(s))

    output = {
        "projects":  projects,
        "partition": {
            "exceptions":  part["exceptions"],
            "bin_diseases": part["bin_diseases"],
        },
        "summary": _summarize(all_samples),
    }

    out_path.write_text(json.dumps(output, indent=2))
    print(f"\nWritten: {out_path}")
    _write_partition_tsv(all_samples, out_path.parent)
    _print_report(output, all_samples)


if __name__ == "__main__":
    main()
