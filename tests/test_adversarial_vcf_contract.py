import gzip
import importlib.util
from pathlib import Path

import pytest
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

def load():
    from vcf_utils import aggregation
    return aggregation

def write_vcf(path, samples, rows):
    with gzip.open(path, "wt") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=GT>\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description=AD>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(samples)+"\n")
        for row in rows: out.write(row+"\n")

def test_sample_order_and_missing_zero_ad_are_not_conflated(tmp_path):
    mod=load(); path=tmp_path/"mutect2.vcf.gz"
    write_vcf(path, ["NORMAL","TUMOR"], ["chr1\t10\t.\tA\tC\t30\tPASS\t.\tGT:AD\t0/0:20,0\t0/1:20,3"])
    rows=mod.read_variants_from_vcf(path, "mutect2")
    assert rows[("1:10:A:C")]["genotype"]["AD"] == "20,3"
    missing=tmp_path/"missing.vcf.gz"
    write_vcf(missing, ["NORMAL","TUMOR"], ["chr1\t10\t.\tA\tC\t30\tPASS\t.\tGT:AD\t0/0:.\t0/1:."])
    row=mod.read_variants_from_vcf(missing, "mutect2")[("1:10:A:C")]
    assert row["genotype"].get("AD") in (None, ".")

def test_duplicate_allele_records_fail_closed(tmp_path):
    mod=load(); path=tmp_path/"duplicate.vcf.gz"
    write_vcf(path, ["TUMOR"], ["chr1\t10\t.\tA\tC\t30\tPASS\t.\tGT:AD\t0/1:20,3", "chr1\t10\t.\tA\tC\t31\tPASS\t.\tGT:AD\t0/1:20,4"])
    with pytest.raises(ValueError, match="Duplicate"):
        mod.read_variants_from_vcf(path, "mutect2")
