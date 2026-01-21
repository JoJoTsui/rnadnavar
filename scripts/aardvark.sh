#!/usr/bin/env bash

FA="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
HC_TRUTH="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/benchmark/seqc2/somatic/release/latest/high-confidence_sSNV+INDEL_in_HC_regions_v1.2.1.vcf.gz"
HC_BED="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/benchmark/seqc2/somatic/release/latest/High-Confidence_Regions_v1.2.bed"

DS_VCF="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/output/seqc2.wes.ll/variant_calling/deepsomatic/WES_LL_T_1_vs_WES_LL_N_1/WES_LL_T_1_vs_WES_LL_N_1.deepsomatic.vcf.gz"
M2_VCF="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/output/seqc2.wes.ll/variant_calling/mutect2/WES_LL_T_1_vs_WES_LL_N_1/WES_LL_T_1_vs_WES_LL_N_1.mutect2.filtered.vcf.gz"
S2_VCF="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/output/seqc2.wes.ll/variant_calling/strelka/WES_LL_T_1_vs_WES_LL_N_1/WES_LL_T_1_vs_WES_LL_N_1.strelka.variants.vcf.gz"

OD="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/comparison/wes/aardvark"
mkdir -p "${OD}"

# aardvark compare \
#     --reference "${FA}" \
#     --truth-vcf "${HC_TRUTH}" \
#     --query-vcf "${DS_VCF}" \
#     --regions "${HC_BED}" \
#     --output-dir "${OD}" \
#     --query-sample "WES_LL_T_1"

# aardvark compare \
#     --reference "${FA}" \
#     --truth-vcf "${HC_TRUTH}" \
#     --query-vcf "${M2_VCF}" \
#     --regions "${HC_BED}" \
#     --output-dir "${OD}" \
#     --query-sample "WES_LL_T_1"

# aardvark compare \
#     --reference "${FA}" \
#     --truth-vcf "${HC_TRUTH}" \
#     --query-vcf "${S2_VCF}" \
#     --regions "${HC_BED}" \
#     --output-dir "${OD}" \
#     --query-sample "WES_LL_T_1"

# micromamba run -n happy som.py \
#     "${HC_TRUTH}" \
#     "${DS_VCF}" \
#     -f "${HC_BED}" \
#     -o "${OD}/DS" \
#     -r "${FA}"

# micromamba run -n happy som.py \
#     "${HC_TRUTH}" \
#     "${M2_VCF}" \
#     -f "${HC_BED}" \
#     -o "${OD}/M2" \
#     -r "${FA}"

# micromamba run -n happy som.py \
#     "${HC_TRUTH}" \
#     "${S2_VCF}" \
#     -f "${HC_BED}" \
#     -o "${OD}/S2" \
#     -r "${FA}"

micromamba run -n happy som.py \
    "${HC_TRUTH}" \
    "/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/comparison/wes/0000.vcf.gz" \
    -f "${HC_BED}" \
    -o "${OD}/C3" \
    -r "${FA}"