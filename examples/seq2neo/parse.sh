#!/usr/bin/env bash
# Parse all three manifests → merged.json + set TSVs.
# seq2neo_root (data output dir) is read from config/runner.yaml.
#
# Usage:
#   bash parse.sh
#   bash /any/path/rnadnavar/examples/seq2neo/parse.sh

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "$HERE/scripts/parse_projects_to_json.py" \
    --config-yaml "$HERE/config/runner.yaml"
