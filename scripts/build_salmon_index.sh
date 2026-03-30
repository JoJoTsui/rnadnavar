#!/usr/bin/env bash
# ==============================================================================
# build_salmon_index.sh
#
# Downloads GENCODE reference files and builds a decoy-aware Salmon index.
# Compatible with Salmon v1.11.x (SSHash-based index).
#
# Usage:
#   bash build_salmon_index.sh [OPTIONS]
#
# Options:
#   --skip-download     Skip downloading reference files (use existing local files)
#   --dry-run           Print resolved paths and commands without executing anything
#   --threads N         Number of CPU threads (default: 12)
#   --aria-conn N       aria2c connections per download (default: 16)
#   --gencode-version N GENCODE release version (default: 49)
#   --outdir DIR        Directory to write index and intermediate files (default: .)
#   --cleanup           Remove intermediate files after successful index build
#   -h, --help          Show this help message
#
# Examples:
#   # Preview all paths and commands without running anything
#   bash build_salmon_index.sh --dry-run --outdir /data/salmon
#
#   # Full run: download + build
#   bash build_salmon_index.sh --outdir /data/salmon
#
#   # Rebuild index only (files already downloaded)
#   bash build_salmon_index.sh --skip-download --outdir /data/salmon
#
#   # Use a different GENCODE version
#   bash build_salmon_index.sh --gencode-version 48 --outdir /data/salmon
# ==============================================================================
set -euo pipefail

# ==============================================================================
# Defaults
# ==============================================================================
GENCODE_VERSION="49"
THREADS=12
ARIA_CONN=16
OUTDIR="."
SKIP_DOWNLOAD=false
CLEANUP=false
DRY_RUN=false

# ==============================================================================
# Argument Parsing
# ==============================================================================
usage() {
    sed -n '/^# Usage:/,/^# ==/p' "$0" | grep "^#" | sed 's/^# \{0,2\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-download)  SKIP_DOWNLOAD=true ;;
        --dry-run)        DRY_RUN=true ;;
        --cleanup)        CLEANUP=true ;;
        --threads)        THREADS="$2";        shift ;;
        --aria-conn)      ARIA_CONN="$2";      shift ;;
        --gencode-version) GENCODE_VERSION="$2"; shift ;;
        --outdir)         OUTDIR="$2";         shift ;;
        -h|--help)        usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
    shift
done

# ==============================================================================
# Derived paths (all relative to OUTDIR)
# ==============================================================================
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}"
TRANSCRIPTS_FILE="${OUTDIR}/gencode.v${GENCODE_VERSION}.transcripts.fa.gz"
GENOME_FILE="${OUTDIR}/GRCh38.primary_assembly.genome.fa.gz"
GENTROME_FILE="${OUTDIR}/gentrome.v${GENCODE_VERSION}.fa.gz"
DECOYS_FILE="${OUTDIR}/decoys.v${GENCODE_VERSION}.txt"

# Detect salmon version (major.minor only, e.g. "1.11") for index naming.
# Indices built with v1.11.x use a new SSHash format incompatible with v1.10.x.
_detect_salmon_version() {
    if command -v salmon &>/dev/null; then
        salmon --version 2>&1 | grep -oP '\d+\.\d+' | head -1
    else
        echo "unknown"
    fi
}
SALMON_MAJOR_MINOR="$(_detect_salmon_version)"
INDEX_DIR="${OUTDIR}/salmon_index_gencode_v${GENCODE_VERSION}_salmon_v${SALMON_MAJOR_MINOR}"

# ==============================================================================
# Pre-flight checks
# ==============================================================================
preflight_checks() {
    local missing=0

    if ! command -v salmon &>/dev/null; then
        echo "Error: salmon is not installed or not in PATH."
        missing=1
    fi

    if [[ "${SKIP_DOWNLOAD}" == false ]] && ! command -v aria2c &>/dev/null; then
        echo "Error: aria2c is not installed (required for download). Install with: sudo apt install aria2"
        missing=1
    fi

    [[ $missing -eq 0 ]] || exit 1

    mkdir -p "${OUTDIR}"
}

# ==============================================================================
# Step 1: Download reference files
# ==============================================================================
download_references() {
    echo "[1/2] Downloading GENCODE v${GENCODE_VERSION} references..."

    aria2c -c -x "${ARIA_CONN}" -s "${ARIA_CONN}" \
        -d "${OUTDIR}" \
        "${BASE_URL}/$(basename "${TRANSCRIPTS_FILE}")"

    aria2c -c -x "${ARIA_CONN}" -s "${ARIA_CONN}" \
        -d "${OUTDIR}" \
        "${BASE_URL}/$(basename "${GENOME_FILE}")"

    echo "    -> Download complete."
}

# ==============================================================================
# Step 2: Build decoy-aware Salmon index
# ==============================================================================
build_index() {
    echo "[2/2] Building decoy-aware Salmon index (v${GENCODE_VERSION})..."

    # Validate inputs exist
    for f in "${TRANSCRIPTS_FILE}" "${GENOME_FILE}"; do
        if [[ ! -f "$f" ]]; then
            echo "Error: Required file not found: $f"
            echo "       Run without --skip-download, or place the file manually."
            exit 1
        fi
    done

    echo "    [2a] Extracting decoy sequence IDs from genome..."
    zcat "${GENOME_FILE}" \
        | grep "^>" \
        | cut -d " " -f 1 \
        | sed 's/>//g' \
        > "${DECOYS_FILE}"

    echo "    [2b] Concatenating gentrome (transcriptome first, then genome)..."
    cat "${TRANSCRIPTS_FILE}" "${GENOME_FILE}" > "${GENTROME_FILE}"

    echo "    [2c] Running salmon index (threads=${THREADS})..."
    salmon index \
        -t "${GENTROME_FILE}" \
        -d "${DECOYS_FILE}" \
        -p "${THREADS}" \
        -i "${INDEX_DIR}" \
        --gencode

    echo "    -> Index written to: ${INDEX_DIR}"
}

# ==============================================================================
# Optional cleanup
# ==============================================================================
cleanup_intermediates() {
    echo "[cleanup] Removing intermediate files..."
    rm -f "${TRANSCRIPTS_FILE}" "${GENOME_FILE}" "${GENTROME_FILE}" "${DECOYS_FILE}"
    echo "    -> Done."
}

# ==============================================================================
# Dry-run: print all resolved paths and commands, then exit
# ==============================================================================
dry_run_report() {
    echo "========================================"
    echo " DRY RUN — nothing will be executed"
    echo "========================================"
    echo ""
    echo "Resolved configuration:"
    echo "  GENCODE version : ${GENCODE_VERSION}"
    echo "  Salmon version  : ${SALMON_MAJOR_MINOR}  ($(command -v salmon 2>/dev/null || echo 'not found'))"
    echo "  Threads         : ${THREADS}"
    echo "  aria2c conns    : ${ARIA_CONN}"
    echo "  Skip download   : ${SKIP_DOWNLOAD}"
    echo "  Cleanup         : ${CLEANUP}"
    echo ""
    echo "Resolved paths:"
    echo "  OUTDIR          : $(realpath -m "${OUTDIR}")"
    echo "  TRANSCRIPTS     : $(realpath -m "${TRANSCRIPTS_FILE}")"
    echo "  GENOME          : $(realpath -m "${GENOME_FILE}")"
    echo "  GENTROME        : $(realpath -m "${GENTROME_FILE}")"
    echo "  DECOYS          : $(realpath -m "${DECOYS_FILE}")"
    echo "  INDEX_DIR       : $(realpath -m "${INDEX_DIR}")"
    echo ""

    if [[ "${SKIP_DOWNLOAD}" == false ]]; then
        echo "Commands that would run:"
        echo ""
        echo "  [download transcripts]"
        echo "  aria2c -c -x ${ARIA_CONN} -s ${ARIA_CONN} \\"
        echo "      -d $(realpath -m "${OUTDIR}") \\"
        echo "      ${BASE_URL}/$(basename "${TRANSCRIPTS_FILE}")"
        echo ""
        echo "  [download genome]"
        echo "  aria2c -c -x ${ARIA_CONN} -s ${ARIA_CONN} \\"
        echo "      -d $(realpath -m "${OUTDIR}") \\"
        echo "      ${BASE_URL}/$(basename "${GENOME_FILE}")"
        echo ""
    else
        echo "  [download] SKIPPED (--skip-download)"
        echo ""
    fi

    echo "  [extract decoys]"
    echo "  zcat $(realpath -m "${GENOME_FILE}") | grep '^>' | cut -d ' ' -f 1 | sed 's/>//g' \\"
    echo "      > $(realpath -m "${DECOYS_FILE}")"
    echo ""
    echo "  [build gentrome]"
    echo "  cat $(realpath -m "${TRANSCRIPTS_FILE}") $(realpath -m "${GENOME_FILE}") \\"
    echo "      > $(realpath -m "${GENTROME_FILE}")"
    echo ""
    echo "  [salmon index]"
    echo "  salmon index \\"
    echo "      -t $(realpath -m "${GENTROME_FILE}") \\"
    echo "      -d $(realpath -m "${DECOYS_FILE}") \\"
    echo "      -p ${THREADS} \\"
    echo "      -i $(realpath -m "${INDEX_DIR}") \\"
    echo "      --gencode"
    echo ""

    if [[ "${CLEANUP}" == true ]]; then
        echo "  [cleanup] Would remove:"
        echo "      $(realpath -m "${TRANSCRIPTS_FILE}")"
        echo "      $(realpath -m "${GENOME_FILE}")"
        echo "      $(realpath -m "${GENTROME_FILE}")"
        echo "      $(realpath -m "${DECOYS_FILE}")"
        echo ""
    fi

    echo "========================================"
    echo " End of dry run. No files were created."
    echo "========================================"
}

# ==============================================================================
# Main
# ==============================================================================
main() {
    if [[ "${DRY_RUN}" == true ]]; then
        dry_run_report
        exit 0
    fi

    echo "========================================"
    echo " Salmon Index Builder"
    echo " GENCODE v${GENCODE_VERSION} | Salmon v${SALMON_MAJOR_MINOR}"
    echo " outdir: ${OUTDIR}"
    echo " index:  ${INDEX_DIR}"
    echo " skip-download: ${SKIP_DOWNLOAD}"
    echo "========================================"

    preflight_checks

    if [[ "${SKIP_DOWNLOAD}" == false ]]; then
        download_references
    else
        echo "[1/2] Skipping download (--skip-download set)."
    fi

    build_index

    if [[ "${CLEANUP}" == true ]]; then
        cleanup_intermediates
    fi

    echo ""
    echo "Done. Index ready at: ${INDEX_DIR}"
}

main
