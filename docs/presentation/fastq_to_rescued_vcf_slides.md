# rnadnavar: FASTQ to Rescued VCF
Conference overview for mixed bioinformatics and AI-agent audience

## Why paired DNA and RNA
- DNA and RNA carry complementary variant evidence
- Consensus reduces single-caller noise
- Rescue keeps cross-modality supported variants

## Input labeling drives workflow routing
- status 0: DNA normal (COO8801DN)
- status 1: DNA tumor (COO8801DT)
- status 2: RNA tumor (COO8801RT)
- Same patient id links all rows

Run example input:
- /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/input/test.rdv.shared.csv

## End-to-end flow
- FASTQ -> alignment
- preprocessing -> variant calling
- normalization -> consensus
- rescue -> annotation/filtering
- optional RNA realignment -> second rescue

Reference diagram:
- ../images/rnadnavar_schemav3.png

## What consensus does
- Aggregates caller VCFs within one modality
- Applies thresholds by type
  - SNV threshold: 2 in this run
  - indel threshold: 2 in this run
- Assigns unified FILTER class per record

## What rescue does
- Merges DNA and RNA consensus evidence
- Uses individual caller support for disagreement handling
- Produces rescued VCF for cross-modality-supported records
- Optional second rescue uses realigned RNA consensus

## FILTER interpretation (slide-level)
In this presentation, focus on:
- Somatic
- Germline
- Reference

Detailed classes (Artifact, NoConsensus, RNAedit) are documented in:
- ../FASTQ_TO_RESCUED_VCF_FILTER_GUIDE.md

## Real output evidence (COO8801.shared)
Consensus:
- /t9k/mnt/hdd/work/Vax/sequencing/aim_exp/rdv_test/output/COO8801.shared/consensus/COO8801DT_vs_COO8801DN/COO8801DT_vs_COO8801DN.consensus.vcf.gz

First rescue:
- /t9k/mnt/hdd/work/Vax/sequencing/aim_exp/rdv_test/output/COO8801.shared/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_vs_COO8801DN.rescued.vcf.gz

Second rescue (realigned RNA):
- /t9k/mnt/hdd/work/Vax/sequencing/aim_exp/rdv_test/output/COO8801.shared/vcf_realignment/rescue/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN/COO8801DT_vs_COO8801DN_rescued_COO8801RT_realign_vs_COO8801DN.rescued.vcf.gz

## Run settings used for this verification
From params_2026-03-18_20-24-15.json:
- rna = true
- dna = true
- tools include consensus,rescue,realignment,vep
- realignment_mode = vcf
- rescue_snv_thr = 2
- rescue_indel_thr = 2

## Takeaway
- Pipeline logic and docs are now aligned
- No algorithm or code changes in this documentation update
- Use the canonical guide for full FILTER taxonomy and troubleshooting
