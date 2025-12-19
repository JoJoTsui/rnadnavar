#!/usr/bin/env bash



# illumina bams
DEST="./google_illumina_bams/"
mkdir -p "${DEST}"
gsutil -m rsync -r \
  "gs://brain-genomics-public/publications/park2024_deepsomatic/bams/illumina/bwa_mem2_grch38_bams/" \
  "${DEST}"



# benchmark vcfs and beds
BENCHMARK="./benchmarking/"
mkdir -p "${BENCHMARK}"
gsutil -m rsync -r \
  "gs://brain-genomics-public/publications/park2024_deepsomatic/benchmarking/" \
  "${BENCHMARK}"
