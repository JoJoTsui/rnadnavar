# SEQC2 rescue FP investigation — 2026-09-14

Status: **exploratory, not promoted to production**. Starting code: `70e0db58`.
Extends the [corrected historical replay](SEQC2_HISTORICAL_NATIVE_GATE_REPLAY.md)
using the [verified source bundle](SEQC2_VERIFIED_OUTPUTS_AND_ARCHIVE.md).
Authoritative historical comparisons remain unchanged.

## Result

A DNA variant-nomination requirement plus biological vetoes removes two WES
rescue FPs and all three WGS rescue FPs without losing any scored rescue TP.
This VCF-only experiment uses the historical native baseline and completed
second/realignment-rescue candidates, not a new workflow run.

| Dataset / target | DeepSomatic TP/FP/FN | Candidate TP/FP/FN | DeepSomatic F1 | Candidate F1 |
| --- | ---: | ---: | ---: | ---: |
| WES-LL / UKB | 1048/38/1252 | 1062/38/1238 | 0.619019 | 0.624706 |
| WES-LL / MedExome | 563/19/266 | 570/19/259 | 0.798016 | 0.803949 |
| WGS-IL / UKB | 2168/19/132 | 2169/19/131 | 0.966347 | 0.966578 |
| WGS-IL / MedExome | 702/6/127 | 702/6/127 | 0.913468 | 0.913468 |

These are observed improvements over this comparator, not proof of SOTA.
WGS UKB gains only one TP; MedExome ties. Both datasets represent the same
cell-line pair, not independent biological validation.

## Why the historical gate regressed

The historical gate requires `N_DNA_CALLERS_SUPPORT >= 1` and
`N_RNA_CALLERS_SOMATIC >= 2`. The actual VCF header and writer define
`N_DNA_CALLERS_SUPPORT` as **observed** DNA callers, not eligible Somatic
votes. Eligibility and Somatic agreement have separate fields. The preceding
replay report incorrectly called this eligible support; that wording is now
corrected without changing the reproduced allele sets.

All three scored WGS rescue FPs have only a DNA DeepSomatic `RefCall`, QUAL=0,
no DNA Mutect2/Strelka record, and zero DNA Somatic votes. RNA agreement plus
a DNA reference observation was enough to pass the historical export gate.

| WGS false positive | DNA tumor AD / DP | Other evidence in saved VCFs |
| --- | --- | --- |
| chr1:10460321 T>C | 39,3 / 42 | REDIportal EDHSAAAK8017, canonical editing annotation, AluSx1 |
| chr17:18662096 A>G | 79,4 / 83 | RNA Mutect2 absent; RNA Strelka LowEVS becomes PASS after realignment; normal Strelka has 0 tier-1 versus 5 tier-2 alternate reads |
| chr2:169812059 C>A | 85,3 / 88 | Initial RNA Mutect2 haplotype/orientation/strand-bias rejection; realignment PASS still has alternate reads on one strand; phased with adjacent site |

The chr17 site also lies in MedExome. Editing/mapping/haplotype mechanisms are
plausible explanations, not read-level diagnoses: no BAM pileup or independent
validation was performed here. Realigning the same RNA reads is not an
independent biological observation.

Two removable WES FPs expose complementary weaknesses:

- chr7:140456191 C>T: only DNA DeepSomatic RefCall, QUAL=0, AD=40,2;
  no other DNA nomination. RNA strand/orientation warnings precede realignment.
- chr3:128620215 A>T: `GNOMAD_AF=0.0995979` despite the source rescue export's
  Somatic label. The historical manual gate does not apply a population veto.
  This identifies a gate-level omission, not the exact upstream stage that
  created the contradictory label.

Remaining WES rescue FPs are chr17:76353649 C>T and chr19:33206715 T>C.
They have low-level DNA Strelka nominations without a decisive common-AF/editing
veto. A VAF threshold chosen solely to remove these two would be post-hoc tuning.

## Tested rule and scope

```text
Entire verified historical native baseline ----------------------> retain

Historical gated-rescue additions (Somatic/PASS SNPs only)
  |
  +-- DNA nomination with positive tumor alternate count?
  |     DeepSomatic: raw FILTER=PASS, or
  |     Mutect2/Strelka: variant record, even if caller-filtered
  |     Missing tumor evidence / RefCall alone -------------------> reject
  |
  +-- GNOMAD_AF > 0.001? -----------------------------------------> reject
  |
  +-- REDI_ACCESSION present AND REDI_CANONICAL=YES
  |     AND N_DNA_CALLERS_SOMATIC=0? ------------------------------> reject
  |
  +--------------------------------------------------------------> add
```

The RNA >=2 and observed DNA >=1 conditions are inherited from the historical
gate input. Nomination is deliberately weaker than an eligible Somatic vote.
Requiring a DNA Somatic vote would retain only 1 of 11 scored WES rescue TPs;
most have low alternate counts or caller-filtered DNA evidence.

Tumor AD uses the matching ALT index; Strelka uses tier-1 base counts. Normal
counts never substitute for tumor evidence. This forensic tool recognizes the
verified SEQC2 sample names, not arbitrary cohort sample roles. Missing AF does
not mean AF=0 and does not trigger this veto; malformed AF fails the assay.
Annotation absence does not establish biological safety. Vetoes apply only to
additions, never prune the retained native baseline. RaVeX_FILTER is not used.

Indels remain the historical DeepSomatic-derived baseline for parity. This is
not validation or acceptance of a new evidence-based indel consensus rule.

## Ablations and per-type metrics

Incremental scored additions relative to native, UKB TP/FP:

| Gate | WES | WGS |
| --- | ---: | ---: |
| Historical | 11/4 | 0/3 |
| DNA nomination only | 11/3 | 0/0 |
| DNA nomination + biological veto | 11/2 | 0/0 |

Both nomination variants were benchmarked by som.py in both regions: eight
completed runs. Final candidate metrics:

| Dataset / target / type | TP/FP/FN | Precision | Recall | F1 |
| --- | ---: | ---: | ---: | ---: |
| WES / UKB / SNP | 1021/34/1184 | 0.967773 | 0.463039 | 0.626380 |
| WES / UKB / indel | 41/4/54 | 0.911111 | 0.431579 | 0.585714 |
| WES / UKB / records | 1062/38/1238 | 0.965455 | 0.461739 | 0.624706 |
| WES / MedExome / SNP | 547/17/248 | 0.969858 | 0.688050 | 0.805004 |
| WES / MedExome / indel | 23/2/11 | 0.920000 | 0.676471 | 0.779661 |
| WES / MedExome / records | 570/19/259 | 0.967742 | 0.687575 | 0.803949 |
| WGS / UKB / SNP | 2084/8/121 | 0.996176 | 0.945125 | 0.969979 |
| WGS / UKB / indel | 85/11/10 | 0.885417 | 0.894737 | 0.890052 |
| WGS / UKB / records | 2169/19/131 | 0.991316 | 0.943043 | 0.966578 |
| WGS / MedExome / SNP | 674/2/121 | 0.997041 | 0.847799 | 0.916383 |
| WGS / MedExome / indel | 28/4/6 | 0.875000 | 0.823529 | 0.848485 |
| WGS / MedExome / records | 702/6/127 | 0.991525 | 0.846803 | 0.913468 |

Selection operates on **all** historical additions before truth scoring:
WES 303 additions become 77; WGS 610 become 131. Only 15 WES and 3 WGS historical
additions are scored in HC intersect UKB. Most excluded sites are outside that
domain; their correctness and potential TP losses remain unresolved.
Baseline retention is exact: 1087 WES and 2188 WGS alleles.

## Evidence and reproduction

Output root: `examples/seqc2/comparison/rescue_fp_investigation_20260914/`.

- `evidence.json`: original DNA, RNA, RNA-realignment caller fields and source
  SHA-256 identities for scored additions; missing records remain absent.
- `sites.tsv`: allele truth attribution, checked against historical som.py deltas.
- `gate_tests/evaluation.json`: source hashes, selection/removal lists, exact
  commands, and SNP/indel/record metrics for all eight tests.
- `gate_tests/{wes_ll,wgs_il}/nomination_biological/query.vcf.gz`: candidate
  PASS allele-only benchmark queries, **not annotated training labels**.
- Corresponding `{ukb,medexome}/benchmark.metrics.json` and `benchmark.log`:
  original benchmark evidence. Nomination-only controls remain alongside.

Repeat from the repository root, using fresh output directories:

```bash
.venv/bin/python examples/seqc2/scripts/investigate_verified_rescue.py --outdir examples/seqc2/comparison/rescue_evidence_NEW
.venv/bin/python examples/seqc2/scripts/test_rescue_nomination_gate.py --outdir examples/seqc2/comparison/rescue_gate_NEW
.venv/bin/python -m pytest tests/test_rescue_nomination_experiment.py tests/test_historical_native_gate_replay.py -q
```

Both scripts read the fixed verified 20260914 bundle: forensic assays, not
portable cohort launchers. Candidate selection does not read truth, but the
hypothesis was informed by these errors and is not held-out validation.
som.py uses the replay's identical truth, assembly38 FASTA, HC `-R`, UKB/MedExome
`-T`, and `-N`, without `-P`. Historical native candidate scope is retained,
including WGS's 2190-position universe and HC/UKB preselection; see the replay.
Twelve focused assay/parity tests passed.

## Next decision

Freeze this candidate for independent validation rather than tune thresholds
to the remaining FPs. Before deployment, resolve out-of-domain losses,
candidate-universe portability, generic sample-role parsing, provenance and
FILTER/rationale preservation, and the separate indel policy. Confirm parity
through consensus/rescue-only execution before changing workflow defaults.
HG008 was not inspected or tuned. Mapping, calling, caches, original inputs,
archives, and authoritative historical results were not modified.
