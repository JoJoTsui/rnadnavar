#!/usr/bin/bash


micromamba run -n rnadnavar python \
  /t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/scripts/convert_rediportal_to_vcf.py \
  -i TABLE1_hg38_v3.txt.gz \
  -o REDIportal_hg38_v3.vcf.gz
