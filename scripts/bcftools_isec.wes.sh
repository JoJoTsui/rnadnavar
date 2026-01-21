#!/usr/bin/env bash

# IO
WD="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/benchmark/seqc2/somatic/release/latest"
AMPLISEQ_D="${WD}/ampliseq"
AMPLISEQ_BED="${AMPLISEQ_D}/AmpliSeq.ComprehensiveCancer.dna_manifest.20220908.bed"
AMPLISEQ_TRUTH="${AMPLISEQ_D}/ampliseq.high-confidence_sSNV+INDEL_in_HC_regions_v1.2.1.vcf.gz"


DS_VCF="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/output/seqc2.wes.ll/variant_calling/deepsomatic/WES_LL_T_1_vs_WES_LL_N_1/WES_LL_T_1_vs_WES_LL_N_1.deepsomatic.vcf.gz"
M2_VCF="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/output/seqc2.wes.ll/variant_calling/mutect2/WES_LL_T_1_vs_WES_LL_N_1/WES_LL_T_1_vs_WES_LL_N_1.mutect2.filtered.vcf.gz"
S2_VCF="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/output/seqc2.wes.ll/variant_calling/strelka/WES_LL_T_1_vs_WES_LL_N_1/WES_LL_T_1_vs_WES_LL_N_1.strelka.variants.vcf.gz"

COMP_DIR="/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seqc2_benchmark/comparison"
DS_VCF_A="${COMP_DIR}/HCC1395_vs_HCC1395BL.deepsomatic.ampliseq.vcf.gz"
M2_VCF_A="${COMP_DIR}/HCC1395_vs_HCC1395BL.mutect2.ampliseq.vcf.gz"
S2_VCF_A="${COMP_DIR}/HCC1395_vs_HCC1395BL.strelka.ampliseq.vcf.gz"

# comparison using high confidence region
bcftools isec -n +3 -p wes -f 'PASS' -O z \
  "${DS_VCF}" \
  "${M2_VCF}" \
  "${S2_VCF}"



# # subset the truth VCF based on given region bed
# bcftools view -R "${AMPLISEQ_BED}" \
#   "${DS_VCF}" \
#   -O z -o "${DS_VCF_A}"
# tabix -p vcf "${DS_VCF_A}"

# bcftools view -R "${AMPLISEQ_BED}" \
#   "${M2_VCF}" \
#   -O z -o "${M2_VCF_A}"
# tabix -p vcf "${M2_VCF_A}"

# bcftools view -R "${AMPLISEQ_BED}" \
#   "${S2_VCF}" \
#   -O z -o "${S2_VCF_A}"
# tabix -p vcf "${S2_VCF_A}"

# # comparison
# bcftools isec -n +3 -p ampliseq -f 'PASS' -O z \
#   "${DS_VCF_A}" \
#   "${M2_VCF_A}" \
#   "${S2_VCF_A}"
