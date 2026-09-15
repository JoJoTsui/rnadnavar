# Native consensus policy optimization

Date: 2026-09-15  \  
Branch: `seqc2-consolidated`

## Decision

The default native SNV policy is now a DeepSomatic-preserving backbone with a narrowly qualified Mutect2 addition. Native mode is post-calling VCF classification; it does not alter mapping, preprocessing, or variant-caller processes.

### SNV decision logic

```text
                 candidate SNV
                       |
          DeepSomatic FILTER PASS/. ?
                 | yes          | no
                 v              v
             Somatic      DeepSomatic QUAL > 0
                                  |
                           Mutect2 TLOD >= 12
                           Mutect2 GERMQ >= 60
                                  |
                    Mutect2 veto combination present?
                         | no                 | yes
                         v                    v
                      Somatic            NoConsensus
```

A failed native-evidence SNV is explicitly `NoConsensus`; it must not fall through to ordinary majority voting. A candidate without DeepSomatic evidence is not admitted by native mode.

The two Mutect2 veto combinations are:

- `contamination;germline;haplotype;panel_of_norms`
- `contamination;orientation;weak_evidence`

The population-frequency and DNA-verification vetoes remain enforced by the existing classifier.

### Indel decision logic

Indels remain on the ordinary configured caller-threshold rule (`2-of-3` by default). They are not promoted by native SNV evidence. The HG008 result supports retaining this separation: threshold consensus produced better indel F1 than a DeepSomatic-only indel backbone.

## Defect corrected

Previously, native evidence was only a promotion path. When native evidence failed, the classifier continued into ordinary majority voting, so a non-native SNV could still become Somatic through 2-of-3 agreement. In addition, normalized DeepSomatic biological labels were treated as equivalent to native `FILTER=PASS`, allowing non-PASS records into the backbone.

The implementation now records `native_evidence_enabled`, requires native `FILTER=PASS`/`.` for baseline retention, and closes the majority-vote fallback for failed native SNVs. Legacy behavior remains available when native mode is disabled.

## Validation

Focused tests and the complete `vcf_utils` suite passed:

```text
6 focused consensus/historical tests passed
207 vcf_utils tests passed, 54 warnings
```

HG008 standalone validation used the completed DeepSomatic, Mutect2, and Strelka VCFs without rerunning mapping or calling. The generated artifacts are retained in `/tmp/hg008-native-consensus-20260915-fixed/`; the benchmark used the authoritative HG008 somatic truth VCF, HG008 truth BED, the shared UKB target BED, the pipeline GRCh38 FASTA, and `som.py -N`.

| Method | SNP F1 | Indel F1 | Record F1 |
| --- | ---: | ---: | ---: |
| Previous native policy | 0.9025 | 0.5369 | 0.8061 |
| Corrected native policy | 0.9566 | 0.5369 | 0.8401 |
| DeepSomatic | 0.9566 | 0.4051 | 0.8103 |

Corrected HG008 counts were 408/15/22 (SNP TP/FP/FN), 91/5/152 (indel TP/FP/FN), and 499/20/170 (record TP/FP/FN).

The previously corrected SEQC2 replay remains the cross-dataset reference: WES-LL UKB F1 0.6206 vs DeepSomatic 0.6190, WES-LL MedExome 0.8000 vs 0.7980, WGS-IL UKB 0.9666 vs 0.9663, and WGS-IL MedExome 0.9135 vs 0.9135 (tie).

## Cache and workflow boundary

The change is consumed after caller VCF production. It does not change process inputs, mapping arguments, caller arguments, work directories, conda environments, or caller cache keys. The ordinary three-FASTQ workflow remains structurally runnable; only the default native consensus labeling semantics are tightened.
