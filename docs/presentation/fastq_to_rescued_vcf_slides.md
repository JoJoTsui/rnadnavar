# rnadnavar: FASTQ to Final FILTER Labels
Presentation focus: how each variant obtains labeled classification from raw input to final VCF FILTER.

## 1. Workflow introduction and routing
- status 0 = DNA normal, status 1 = DNA tumor, status 2 = RNA tumor
- FASTQ -> alignment -> calling -> consensus -> rescue -> final FILTER
- DNA and RNA are processed in parallel then merged in rescue

## 2. Per-caller label normalization
- DeepSomatic, Mutect2, and Strelka are normalized to shared biological labels
- Strelka uses FILTER + NT + normal depth for mapping
- These per-caller labels become consensus input

## 3. Within-modality consensus logics
- Aggregate by normalized variant key
- Apply SNV/indel thresholds (2/2 in this run)
- Majority class becomes consensus label
- Tie -> Artifact; below threshold -> NoConsensus

## 4. Cross-modality rescue logics
- Combine DNA-consensus and RNA-consensus labels
- Resolve agreement/disagreement using support-aware rules
- Optional second rescue uses realigned RNA consensus

## 5. Final VCF FILTER categories in this deck
- Somatic
- Germline
- Reference

## 6. Consistency and traceability summary
- Final FILTER is unified classification label for interpretation
- Label path is traceable from routing to caller mapping to consensus/rescue
- Slide notes are bilingual (English + Chinese) on every slide
