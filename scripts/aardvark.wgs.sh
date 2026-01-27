#!/usr/bin/env bash

FA="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/references/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
HC_TRUTH="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/benchmark/seqc2/somatic/release/latest/high-confidence_sSNV+INDEL_in_HC_regions_v1.2.1.vcf.gz"
HC_BED="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/benchmark/seqc2/somatic/release/latest/High-Confidence_Regions_v1.2.bed"
# WES_BED="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/bio_db/intervals/illumina_beds/hg38_Twist_Bioscience_for_Illumina_Exome_2_5_CSPGx.bed"
# HC_WES_BED="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/benchmark/seqc2/somatic/release/latest/High-Confidence_hg38_Twist_Bioscience_for_Illumina_Exome_2_5_CSPGx.bed"

DS_VCF="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/output/seqc2.wgs.ll/variant_calling/deepsomatic/WGS_LL_T_1_vs_WGS_LL_N_1/WGS_LL_T_1_vs_WGS_LL_N_1.deepsomatic.vcf.gz"
M2_VCF="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/output/seqc2.wgs.ll/variant_calling/mutect2/WGS_LL_T_1_vs_WGS_LL_N_1/WGS_LL_T_1_vs_WGS_LL_N_1.mutect2.filtered.vcf.gz"
S2_VCF="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/output/seqc2.wgs.ll/variant_calling/strelka/WGS_LL_T_1_vs_WGS_LL_N_1/WGS_LL_T_1_vs_WGS_LL_N_1.strelka.variants.vcf.gz"


OD="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/comparison/wgs/aardvark.wgs"
mkdir -p "${OD}"


# comparison using high confidence region
bcftools isec -n +3 -p wgs -f 'PASS' -O z \
  "${DS_VCF}" \
  "${M2_VCF}" \
  "${S2_VCF}"

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

micromamba run -n happy som.py \
    "${HC_TRUTH}" \
    "${DS_VCF}" \
    -R "${HC_BED}" \
    -o "${OD}/DS" \
    -r "${FA}" -N

micromamba run -n happy som.py \
    "${HC_TRUTH}" \
    "${M2_VCF}" \
    -R "${HC_BED}" \
    -o "${OD}/M2" \
    -r "${FA}" -N

micromamba run -n happy som.py \
    "${HC_TRUTH}" \
    "${S2_VCF}" \
    -R "${HC_BED}" \
    -o "${OD}/S2" \
    -r "${FA}" -N

# micromamba run -n happy som.py \
#     "${HC_TRUTH}" \
#     "/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/comparison/wes/0000.vcf.gz" \
#     -f "${HC_BED}" \
#     -o "${OD}/C3" \
#     -r "${FA}"