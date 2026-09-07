#!/usr/bin/env python3
"""Build a content-bound provenance manifest for a label artifact."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def identity(path, label):
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"{label} does not exist or is not a file: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": digest.hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--metadata", help="JSON object containing stage, selectors, tools and effective_args")
    parser.add_argument("--reference")
    parser.add_argument("--regions")
    parser.add_argument("--command", action="append", default=[])
    args = parser.parse_args()
    try:
        metadata = json.loads(Path(args.metadata).read_text()) if args.metadata else {}
        if not isinstance(metadata, dict):
            raise ValueError("metadata must be a JSON object")
        result = {
            "schema": "seqc2-artifact-provenance.v1",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "artifact": identity(args.artifact, "artifact"),
            "inputs": {},
            "stage": metadata.get("stage", "unknown"),
            "sample": metadata.get("sample", "unknown"),
            "modality": metadata.get("modality", "unknown"),
            "library": metadata.get("library", "unknown"),
            "selector": metadata.get("selector", "unknown"),
            "effective_args": metadata.get("effective_args", {}),
            "tools": metadata.get("tools", {}),
            "model": metadata.get("model", "unknown"),
            "databases": metadata.get("databases", {}),
            "commands": args.command,
        }
        if not isinstance(result["effective_args"], dict) or not isinstance(result["tools"], dict) or not isinstance(result["databases"], dict):
            raise ValueError("effective_args, tools and databases must be JSON objects")
        if args.reference:
            result["inputs"]["reference"] = identity(args.reference, "reference")
        if args.regions:
            result["inputs"]["regions"] = identity(args.regions, "regions")
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.error(str(error))
    Path(args.out).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
