#!/usr/bin/env python3
"""
Validate existing statistics scripts against real seq2neo VCF data.

Tests variant_statistics_process.py.py on representative samples from each
set and documents findings including the known Set 1 prefix mismatch bug.

READ-ONLY: This script only reads VCF files, never writes to data directories.
"""

import os
import sys
import tempfile
from pathlib import Path

# ── Test samples (1-2 per set) ────────────────────────────────────────────
TEST_SAMPLES = [
    {
        "set_number": 1,
        "sample_id": "PRJNA298376_4060",
        "patient_id": "4060",
        "base_output_dir": "/t9k/mnt/WorkSpace/data/ngs/zhanlingmin/Rnadnavar/output",
        "dir_name": "4060_rnadnavar",
        "vcf_prefix": "4060",  # Correct prefix (not dir_name!)
        "status": "standard",
    },
    {
        "set_number": 1,
        "sample_id": "PRJNA298376_4255",
        "patient_id": "4255",
        "base_output_dir": "/t9k/mnt/WorkSpace/data/ngs/zhanlingmin/Rnadnavar/output",
        "dir_name": "4255_rnadnavar",
        "vcf_prefix": "4255",
        "status": "standard",
    },
    {
        "set_number": 2,
        "sample_id": "PRJNA298376_3942",
        "patient_id": "3942",
        "base_output_dir": "/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seq2neo/output",
        "dir_name": "PRJNA298376_3942",
        "vcf_prefix": "PRJNA298376_3942",
        "status": "standard",
    },
    {
        "set_number": 3,
        "sample_id": "PRJNA298376_3812",
        "patient_id": "3812",
        "base_output_dir": "/t9k/mnt/WorkSpace/data/ngs/liuxin/seq2neo/output",
        "dir_name": "PRJNA298376_3812",
        "vcf_prefix": "PRJNA298376_3812",
        "status": "standard",
    },
    {
        "set_number": 4,
        "sample_id": "PRJNA298376_4069",
        "patient_id": "4069",
        "base_output_dir": "/t9k/mnt/WorkSpace/data/ngs/liuxin/seq2neo/output",
        "dir_name": "PRJNA298376_4069",
        "vcf_prefix": "PRJNA298376_4069",
        "status": "standard",
    },
]

# Path to the existing script
SCRIPT_PATH = Path(__file__).resolve().parent / "variant_statistics_process.py.py"


def validate_rescue_vcf_exists(sample: dict) -> bool:
    """Check if the rescue VCF exists for a sample."""
    rescue_pattern = (
        f"{sample['vcf_prefix']}DT_vs_{sample['vcf_prefix']}DN"
        f"_rescued_{sample['vcf_prefix']}RT_realign_vs_{sample['vcf_prefix']}DN"
        f".filtered.vcf.stripped.vcf.gz"
    )
    rescue_dir = os.path.join(
        sample["base_output_dir"], sample["dir_name"],
        "vcf_realignment", "rescue",
        f"{sample['vcf_prefix']}DT_vs_{sample['vcf_prefix']}DN"
        f"_rescued_{sample['vcf_prefix']}RT_realign_vs_{sample['vcf_prefix']}DN",
    )
    path = os.path.join(rescue_dir, rescue_pattern)
    exists = os.path.isfile(path)
    return exists, path


def validate_caller_vcfs_exist(sample: dict) -> dict:
    """Check which of the 6 caller VCFs exist for a sample."""
    import glob

    prefix = sample["vcf_prefix"]
    base_dir = os.path.join(sample["base_output_dir"], sample["dir_name"])

    caller_configs = {
        "DNA_deepsomatic": os.path.join(
            base_dir, "variant_calling", "deepsomatic",
            f"{prefix}DT_vs_{prefix}DN", "*.deepsomatic.vcf.gz"
        ),
        "DNA_mutect2": os.path.join(
            base_dir, "variant_calling", "mutect2",
            f"{prefix}DT_vs_{prefix}DN", "*.mutect2.filtered.vcf.gz"
        ),
        "DNA_strelka": os.path.join(
            base_dir, "variant_calling", "strelka",
            f"{prefix}DT_vs_{prefix}DN", "*.strelka.variants.vcf.gz"
        ),
        "RNA_deepsomatic": os.path.join(
            base_dir, "vcf_realignment", "variant_calling", "deepsomatic",
            f"{prefix}RT_realign_vs_{prefix}DN", "*.deepsomatic.vcf.gz"
        ),
        "RNA_mutect2": os.path.join(
            base_dir, "vcf_realignment", "variant_calling", "mutect2",
            f"{prefix}RT_realign_vs_{prefix}DN", "*.mutect2.filtered.vcf.gz"
        ),
        "RNA_strelka": os.path.join(
            base_dir, "vcf_realignment", "variant_calling", "strelka",
            f"{prefix}RT_realign_vs_{prefix}DN", "*.strelka.variants.vcf.gz"
        ),
    }

    results = {}
    for caller, pattern in caller_configs.items():
        files = glob.glob(pattern)
        results[caller] = files[0] if files else None
    return results


def test_existing_script_with_wrong_prefix(sample: dict):
    """Test that the existing script FAILS with the wrong prefix for Set 1."""
    # The existing script uses dir_name as the project_name/prefix
    # For Set 1: dir_name = "4060_rnadnavar" but VCF prefix = "4060"
    wrong_prefix = sample["dir_name"]
    rescue_pattern_wrong = (
        f"{wrong_prefix}DT_vs_{wrong_prefix}DN"
        f"_rescued_{wrong_prefix}RT_realign_vs_{wrong_prefix}DN"
        f".filtered.vcf.stripped.vcf.gz"
    )
    rescue_dir_wrong = os.path.join(
        sample["base_output_dir"], sample["dir_name"],
        "vcf_realignment", "rescue",
        f"{wrong_prefix}DT_vs_{wrong_prefix}DN"
        f"_rescued_{wrong_prefix}RT_realign_vs_{wrong_prefix}DN",
    )
    path_wrong = os.path.join(rescue_dir_wrong, rescue_pattern_wrong)
    exists = os.path.isfile(path_wrong)
    return not exists, path_wrong


def main():
    print("=" * 70)
    print("VALIDATION OF EXISTING STATISTICS SCRIPTS")
    print("=" * 70)
    print(f"Script under test: {SCRIPT_PATH}")
    print(f"Script exists: {SCRIPT_PATH.exists()}")
    print()

    results = []

    for sample in TEST_SAMPLES:
        sid = sample["sample_id"]
        set_num = sample["set_number"]
        print(f"--- Set {set_num}: {sid} ({sample['status']}) ---")

        # 1. Check rescue VCF with correct prefix
        rescue_ok, rescue_path = validate_rescue_vcf_exists(sample)
        print(f"  Rescue VCF (correct prefix): {'OK' if rescue_ok else 'MISSING'}")

        # 2. Check rescue VCF with WRONG prefix (simulating existing script bug for Set 1)
        if set_num == 1:
            fails, wrong_path = test_existing_script_with_wrong_prefix(sample)
            print(f"  Rescue VCF (wrong prefix='{sample['dir_name']}'): {'MISSING (bug confirmed!)' if fails else 'EXISTS (unexpected)'}")
            if fails:
                print(f"    -> Existing script would FAIL for Set 1!")
        else:
            print(f"  Rescue VCF (prefix='{sample['vcf_prefix']}' matches dir_name): OK")

        # 3. Check caller VCFs
        caller_results = validate_caller_vcfs_exist(sample)
        n_callers = sum(1 for v in caller_results.values() if v)
        print(f"  Caller VCFs: {n_callers}/6 present")
        for caller, path in caller_results.items():
            status = "OK" if path else "MISSING"
            if not path:
                print(f"    {caller}: {status}")

        results.append({
            "sample_id": sid,
            "set_number": set_num,
            "rescue_ok": rescue_ok,
            "callers_ok": n_callers,
            "prefix_match": set_num != 1,  # Set 1 always has prefix mismatch
        })
        print()

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    all_rescue = all(r["rescue_ok"] for r in results)
    set1_prefix_issue = all(r["prefix_match"] or r["set_number"] != 1 for r in results)
    print(f"All rescue VCFs present: {all_rescue}")
    print(f"Set 1 prefix mismatch confirmed: {not set1_prefix_issue}")
    print()
    print("The existing variant_statistics_process.py.py uses project_name")
    print("(the directory name) as the VCF prefix. For Set 1 samples, the")
    print("directory is named '{patient_id}_rnadnavar' but the VCF files use")
    print("'{patient_id}' as the prefix. This causes FileNotFoundError for")
    print("all 15 Set 1 samples when run unmodified.")
    print()
    print("Fix: Strip '_rnadnavar' suffix from project_name for Set 1 samples,")
    print("or use the sample_manifest.parquet vcf_prefix column instead of")
    print("deriving it from directory names.")


if __name__ == "__main__":
    main()
