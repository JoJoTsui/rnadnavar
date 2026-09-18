# Three-class rerun review and development policy (2026-09-18)

Status: candidate exports and opt-in development rules; biological training
approval remains open. Neither the published SEQC2/HG008 performance claims nor
the completed 66-sample VCFs are changed by these code additions.

## Why DNA consensus and rescue have different negative-label counts

The frozen 2026-09-16 rerun used seqc2_refined_v2 and seqc2_refined_gate_v1.

| FILTER | Original rescue | Refined DNA consensus | Refined rescue |
|---|---:|---:|---:|
| Somatic | 117859 | 32725 | 32933 |
| Germline | 986236 | 15839 | 986181 |
| Reference | 7197707 | 38480 | 7197503 |
| NoConsensus | 114109911 | 105963208 | 114226688 |

The similar absolute sizes of the three positive/negative label groups in DNA
consensus are incidental. There is no class-balancing rule.
Large negative counts alone are not proof of bad labels: somatic candidate
unions naturally contain many more non-Somatic loci, and the rescue union also
contains RNA/realignment-only candidates. The concern is unsupported confidence,
not matching the sizes of the three classes.
The legacy support rule also uses a tumor ALT-read floor, which is appropriate
for ALT support but not for establishing a true reference genotype with zero
ALT reads. The new Reference evidence rule bypasses that Somatic voting floor.
In classification.py, refined SNPs either pass the Somatic rule or become
NoConsensus before biological majority classification. Rejected indels still
reach the legacy classifier. Thus ordinary negative SNP classifications are
suppressed. A synthetic two-caller unanimous Germline or Reference vote
reproduces the suppression.

apply_refined_rescue.py takes the union of DNA consensus and the pre-existing,
annotated realignment rescue. Its retained_rescue_negative branch preserves
125048812 records across all non-Somatic labels, including NoConsensus and
Artifact. It does not independently establish negative-label confidence.
It retains 32725 Somatic DNA baseline records and adds 208 gated Somatic records.
Reported baseline conflicts: 11451 with rescue negatives and 3607 with biological
annotations. These groups may overlap; neither count alone is a measured FP count.

All 132 stored structural-audit counts matched the prior metadata. These audits
check representation/provenance, not biological validity. Prior metadata copied
state hashes; matching copied hashes alone was not a fresh integrity check.
The corrected export hashes its VCF/audit inputs and compares actual exported
class counts against audits.

## Export correction

The old Parquet contains only 32933 Somatic rows, and all five exported rescue/
caller-support/gnomAD fields were null. It is not the complete training table.
The new exporter retains all 8216617 Somatic/Germline/Reference candidate rows
expected from the frozen rescue VCFs. It keeps actual source evidence (including
normal/tumor caller fields), full ALT, decision rationale, and rescue round.
Missing source fields stay null. It does not turn missing evidence into zero.

Heavy outputs are under:
`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/examples/seq2neo/output_refined_native_v2_20260916/variant_parquet/`

- refined_three_class_candidates.parquet
- refined_three_class_manifest.parquet

Lightweight manifest and summary are generated under
`examples/seq2neo/data/processed/refined_three_class_20260918/`.
The sample manifest retains project/patient, BAM, caller and original-rescue
paths. It points rescue_vcf_path to the new candidate VCF and leaves
training_label_vcf empty. Split assignment and previous QC approvals are not
silently inherited. All rows have training_eligible=false until separately
validated. The original Somatic-only export is retained as historical evidence.

Reproduction:

```bash
.venv/bin/python examples/seq2neo/scripts/build_refined_rerun_artifacts.py
```

Use fresh --outdir and --parquet paths for another execution. Outputs are not
overwritten. Data are streamed in bounded batches; reducing the class set to
save memory is not acceptable.

## Truth: what we can and cannot infer

The [SEQC2 reference study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8532138/)
describes matched-normal HCC1395BL germline discovery and validation using
multiple replicates, aligners and germline callers. Therefore it is incorrect
to assert that SEQC2 has no germline reference work. A usable, versioned germline
VCF plus its applicable confidence region must still be identified and checked.
The [currently inspected v1.2.1 release](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/v1.2.1/)
contains somatic SNV/indel benchmarks and their region BED; that BED must not
be assumed to be a germline-confidence mask.

[NIST's HG008 page](https://www.nist.gov/programs-projects/cancer-genome-bottle)
describes v0.3 as a clonal/truncal somatic small-variant benchmark.
The [official PacBio analysis directory](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/PacBio_Revio_20240125/pacbio-wgs-wdl_germline_20240206/)
contains HG008-N-P.GRCh38.deepvariant.phased.vcf.gz and the corresponding
g.vcf.gz, alongside tumor files. These are available orthogonal normal-genotype
evidence, not automatically adjudicated three-class truth. Tissue identity
(N-P versus N-D), reference dictionary, callability and confidence need checking.

For any benchmark, absence from a somatic VCF is not evidence that the person's
genotype is Reference. It may be germline, an artifact, uncallable, an unreported
subclone, or outside the benchmark's detection scope. Somatic truth can reveal
negative labels wrongly covering known Somatic alleles after equivalent-allele
matching, but cannot measure Germline-versus-Reference precision by itself.
Similarly gnomAD frequency is a prior, not this individual's genotype.

A validated Reference training label requires evidence against a particular
candidate ALT in both DNA samples, at a stated detection limit. The experimental
native-only Reference nomination below does not yet certify that requirement.
Neither certifies all other alleles or all unreported genome positions.

## Opt-in development rules

New CLI: run_consensus_vcf.py --experimental-refined-native
--experimental-three-class. Retain RefCall/Germline inputs.
Rescue: apply_refined_rescue.py --experimental-three-class, consuming DNA
consensus generated by that mode. Existing defaults remain unchanged.

```text
Normalized biallelic DNA evidence
  + refined Somatic admission
  + explicit native Germline/Reference nomination and evidence
            |
            + conflicting supported classes --> NoConsensus / review
            + Somatic only ------------------> Somatic candidate
            + Germline only -----------------> Germline candidate
            + Reference only ----------------> Reference candidate
            + missing/inadequate evidence ---> NoConsensus

Rescue + new DNA labels
  + annotation / verification / class conflict --> NoConsensus
  + supported DNA negative --------------------> retain negative candidate
  + inherited negative without native support -> NoConsensus
  + new gated Somatic addition -----------------> also require DNA verification
```

Provisional evidence thresholds, applied to both SNPs and indels:

| Class | Development evidence requirements |
|---|---|
| Germline | Explicit DeepSomatic GERMLINE, native GQ >=30, DP and observed allele-read total >=20, ALT >=3. When present, three-value PL must uniquely favor its Germline class. Mutect2/Strelka normal evidence with >=20 observed reads, ALT >=5 and ALT fraction >=0.2 corroborates; >=60 observed normal reads with zero ALT conflicts. Tumor VAF is not required to be 0.5: purity, LOH and CNV can change it. |
| Reference | Explicit DeepSomatic RefCall, native GQ >=30, DP and observed allele-read total >=20, ALT <=2 and ALT fraction <=0.05. When present, PL must uniquely favor Reference. Mutect2/Strelka normal or tumor evidence with >=10 observed reads, ALT >=3 and ALT fraction >0.05 conflicts. A caller with >=60 observed reads and zero ALT in both samples corroborates. |
| Somatic | Existing refined admission provides nomination; supported negative-class conflict abstains. Annotated rescue conflicts are no longer silently overridden in the new mode. |
| Uncertain | Missing native GQ/AD/DP, malformed/unsupported alleles, low depth, contradictory native likelihoods, class conflict or unsupported inherited negatives abstain. Missing normal GQ alone does not discard a native DeepSomatic class. |

Here GQ is DeepSomatic's native model confidence, not normal-sample genotype
confidence, QUAL, TLOD or GERMQ. A Phred score of 30 conventionally corresponds
to a model error probability of 0.001, not calibrated truth accuracy here.
DP is sample read depth;
AD supplies allele-specific counts. Shared-read caller counts must not be summed
or interpreted as independent experiments.
With zero ALT among 60 independent observations the one-sided binomial 95%
upper bound is 1 - 0.05^(1/60), about 4.87%. This does not exclude a 1% subclone;
doing that under the same ideal model requires about 299 observations.
Mapping, strand, repeats and correlated reads impose further limitations.

These rules are an explicit development baseline, not optimized or biologically
validated confidence thresholds. A native-only candidate is distinguished from
a read-corroborated candidate in CLASSIFICATION_RATIONALE; neither is training
approval. Shared reads do not make corroborating callers independent truth.
Only biallelic A/C/G/T alleles and three-value model PL, when provided, are
supported. No global default change is justified yet.

Policy identity: native_three_class_v1; the rescue counterpart is
native_three_class_gate_v1. Somatic failure is never a negative nomination:
Mutect2 germline/contamination/low-evidence filters and Strelka rejected Somatic
candidates are evidence for review, not independent Germline/Reference truth.
The legacy classifier can map such failures into biological categories; those
mappings are not used to nominate negatives in this new opt-in policy.

DeepSomatic's FILTER, sample-column layout and GT must be distinguished.
In the cohort's 1.9 implementation, its
[somatic VCF writer](https://github.com/google/deepvariant/blob/r1.9/third_party/nucleus/io/vcf_writer.cc)
recodes the native Germline prediction to GT=0/0 and FILTER=GERMLINE. Therefore
GT=0/0 is not sufficient evidence for Reference in these outputs. An actual
cohort record is chr1:63735 CCTA>C, GERMLINE, GT=0/0, GQ=7, DP=2, AD=0,2,
PL=5,0,29: its explicit class is Germline, but its confidence/depth do not pass
the provisional development gate. The policy never requires an ordinary
non-reference GT from DeepSomatic to recognize its Germline nomination.

## Evidence availability found in the cohort

PRJNA298376_3812 DNA input inspection:

- DeepSomatic emits one tumor-named FORMAT column, including GERMLINE and
  RefCall records. This does not mean Somatic-only output or tumor-only calling.
- Strelka somatic FORMAT has base/tier counts, but no paired GT/GQ.
- Mutect2 declares GQ in the header, but all 5579 inspected records had missing
  normal GQ (normal GT was present).

An earlier provisional two-caller normal-GQ requirement was inappropriate for
these somatic-caller outputs and is superseded by the native nomination policy
above. DeepSomatic native model confidence can nominate its explicit classes;
it must not be described as normal-sample genotype confidence. Mutect2/Strelka
read evidence can corroborate or contradict those nominations without normal
GQ. Existing paired BAMs support additional validation without repeating
mapping/variant calling; independent normal gVCFs add orthogonal evidence.

The bounded cohort inventory subsequently sampled the first 200 records of
each of 198 DNA VCFs (39600 records total). All 66 DeepSomatic files were
single-tumor-column VCFs, not Somatic-only VCFs. In 13200 Mutect2 records, normal GT/DP were present in all, AD in
12454, and GQ in zero. Strelka had normal DP but neither GT nor GQ nor generic
AD; it uses native base/tier counts instead. These are prefix samples, not a
random or exhaustive estimate. Results and exact file paths are recorded in
examples/seq2neo/data/processed/three_class_review_20260918/evidence_inventory.json.
Some source indexes had older timestamps than VCFs; this inventory read
sequentially without region queries and did not modify indexes or source files.

The follow-up native_class_inventory.json uses the same bounded prefix design.
Among 13200 DeepSomatic records, 1360 are GERMLINE, 11828 RefCall and 12 PASS;
all 1360 GERMLINE records have recoded GT=0/0. The provisional native-only gate
nominates 255 Germline and 2445 Reference candidates, before any cross-caller
conflict check. The remainder fail nomination/confidence/depth/ALT checks.
These figures validate evidence availability and parsing, not class accuracy,
and must not be extrapolated as cohort retention estimates.

Reproduce in a fresh file (the inventory refuses overwrite):

```bash
.venv/bin/python examples/seq2neo/scripts/audit_three_class_evidence.py --out /tmp/native_class_inventory.json
```

## Evaluation and remaining work

1. Export all existing candidate classes, retain original labels and provenance,
   and stratify counts, missing evidence and annotation conflicts.
2. Validate genotype extraction, missingness, label conflicts, and CLI behavior
   with synthetic paired records and unchanged-policy regressions.
3. Locate/reference-check matched-normal germline benchmarks; retain HG008 for
   evaluation rather than tuning SEQC2 thresholds against its results.
4. Validate negative candidates using available normal genotypes and targeted
   normal/tumor BAM checks. Include map/base quality, strand, local repeats,
   allele representation, depth and tumor purity/CNV strata.
5. Evaluate per-class precision/recall/F1 where true class labels are available;
   report abstention/coverage and the full class confusion matrix. Report
   Somatic-truth overlap of negatives as contamination, not complete negative
   class accuracy. Benchmark SNPs and indels separately in both existing UKB
   and MedExome regions. Group patient-level train/evaluation partitions.
6. Assess the effect of conflict abstention on the 3607/11451 previously reported
   conflict groups and on Somatic TP before changing defaults or cohort labels.

The corrected export is a recovery of the frozen rerun's missing classes.
It is not a new three-class VCF rerun or a claim that all exported negatives
meet the experimental policy. A new training-approved cohort requires the
remaining biological validation.

## Completed export (2026-09-18)

The streaming export finished for all 66 samples: 8216617 rows. Somatic:
31571 SNPs + 1362 indels; Germline: 916950 SNPs + 69231 indels; Reference:
7140338 SNPs + 57165 indels. All rows remain candidate-only. The recorded
Parquet SHA256 is
1ea4fb32b8a7673ed8d36ad5a2621e6eb7cba2159c2119d5e50d415a21300b7b.
This export does not apply the experimental native_three_class_v1 policy.

A second streaming read verified all 8216617 rows, 66 sample IDs and zero
training_eligible=true rows. Non-null evidence counts: RESCUED and
RESCUE_PROMOTED 8216617 each; N_DNA_CALLERS_SUPPORT and N_RNA_CALLERS_SUPPORT
8216522 each; GNOMAD_AF 3449066; CLASSIFICATION_RATIONALE 8216617.
AD_BY_CALLER and NORMAL_AD_BY_CALLER occur in only 95 records each. Inspection
confirmed these fields are absent from an inherited source rescue record even
though declared in its header; this is upstream missing evidence, not a license
to manufacture counts. Original caller VCFs/BAMs remain necessary for validating
most negative candidates. Do not treat this complete-class export as a complete
per-caller-evidence reconstruction.
