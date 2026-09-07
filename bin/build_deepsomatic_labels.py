#!/usr/bin/env python3
"""Build deterministic labels from DeepSomatic PASS and verified RNA nominations."""
import argparse
import gzip
import json
import subprocess
import tempfile
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def key(fields):
    return fields[0], int(fields[1]), fields[3].upper(), fields[4].upper()


def add_rationale(info, rationale):
    return rationale if info in ("", ".") else info + ";" + rationale


def read_records(path):
    headers, records = [], {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                headers.append(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"VCF record has fewer than 8 columns in {path}: {line.strip()!r}")
            records[key(fields)] = fields
    return headers, records


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("--deepsomatic-vcf", required=True)
    parser.add_argument("--rna-nominations")
    parser.add_argument("--verification-json")
    parser.add_argument("--out", required=True)
    parser.add_argument("--index", action="store_true", help="Create a tabix index for .vcf.gz output")
    args = parser.parse_args()
    try:
        headers, records = read_records(args.deepsomatic_vcf)
        verified = {}
        if args.verification_json:
            data = json.loads(Path(args.verification_json).read_text())
            if not isinstance(data, dict) or not isinstance(data.get("results", []), list):
                raise ValueError("verification JSON must contain a results list")
            for evidence in data["results"]:
                if evidence.get("status") != "confirmed":
                    continue
                evidence_key = (str(evidence["chrom"]), int(evidence["pos"]), evidence["ref"].upper(), evidence["alt"].upper())
                verified[evidence_key] = evidence
        for variant_key, fields in list(records.items()):
            filters = {item.upper() for item in fields[6].split(";")}
            if filters & {"PASS", ".", "SOMATIC"}:
                fields[6] = "Somatic"
                fields[7] = add_rationale(fields[7], "CLASSIFICATION_RATIONALE=rule:deepsomatic_pass_starting_set|class:Somatic")
        if args.rna_nominations:
            _, nominations = read_records(args.rna_nominations)
            for variant_key, fields in nominations.items():
                evidence = verified.get(variant_key)
                if evidence is None:
                    continue
                fields[6] = "Somatic"
                fields[7] = add_rationale(fields[7], "CLASSIFICATION_RATIONALE=rule:verified_dna_rna_nomination|class:Somatic|tumor_alt:%s|normal_alt:%s" % (evidence.get("tumor_alt"), evidence.get("normal_alt")))
                records[variant_key] = fields
        output = Path(args.out)
        output.parent.mkdir(parents=True, exist_ok=True)
        body = "".join(headers) + "".join("\t".join(records[item]) + "\n" for item in sorted(records))
        if args.index:
            if not str(output).endswith(".vcf.gz"):
                raise ValueError("--index requires an output ending in .vcf.gz")
            with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as temporary:
                temporary.write(body)
                temporary_path = temporary.name
            try:
                with output.open("wb") as target:
                    subprocess.run(["bgzip", "-c", temporary_path], check=True, stdout=target)
                subprocess.run(["tabix", "-f", "-p", "vcf", str(output)], check=True)
            finally:
                Path(temporary_path).unlink(missing_ok=True)
        elif str(output).endswith(".gz"):
            with gzip.open(output, "wt") as target:
                target.write(body)
        else:
            output.write_text(body)
    except (OSError, ValueError, KeyError, json.JSONDecodeError, subprocess.CalledProcessError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
