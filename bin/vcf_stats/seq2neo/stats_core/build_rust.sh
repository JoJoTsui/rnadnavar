#!/usr/bin/env bash
# Build the Rust stats_core module for the rnadnavar seq2neo statistics pipeline.
#
# Usage:
#   ./build_rust.sh          # Release build (optimized)
#   ./build_rust.sh debug    # Debug build (with assertions, faster compile)
#   ./build_rust.sh release  # Explicit release build
#
# Requirements:
#   - Rust toolchain (cargo, rustc) — https://rustup.rs
#   - maturin >= 1.0 — pip install maturin
#   - Python 3.10+ with the target venv/conda environment active
#
# The built .so is placed alongside this script (in stats_core/), and a symlink
# or copy is placed in the parent directory (seq2neo/) where Python imports it.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MODE="${1:-release}"

case "$MODE" in
    release|--release|-r)
        echo "=== Building stats_core (release) ==="
        maturin develop --release
        ;;
    debug|--debug|-d)
        echo "=== Building stats_core (debug) ==="
        maturin develop
        ;;
    *)
        echo "Usage: $0 [release|debug]"
        echo "  release  — optimized build (default)"
        echo "  debug    — debug build with assertions, faster compile"
        exit 1
        ;;
esac

# Verify the build
TARGET_DIR="target/release"
if [ "$MODE" = "debug" ]; then
    TARGET_DIR="target/debug"
fi

SO_FILE=$(find "$TARGET_DIR" -name "stats_core*.so" 2>/dev/null | head -1)
if [ -n "$SO_FILE" ]; then
    echo "=== Build successful ==="
    echo "Output: $SO_FILE"
else
    echo "=== Build completed but .so not found in $TARGET_DIR ==="
    echo "maturin may have installed it directly into the Python environment."
fi
