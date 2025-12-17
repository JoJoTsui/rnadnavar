#!/usr/bin/env bash


micromamba run -n rnadnavar \
  python ~/hdd/work/Vax/pipeline/rnadnavar/scripts/prepare_cosmic_database.py \
  --input Cosmic_GenomeScreensMutant_Normal_v103_GRCh38.vcf.gz \
  --output ucsc_chr_cosmic_GenomeScreensMutant_Normal_v103_GRCh38.vcf.gz
