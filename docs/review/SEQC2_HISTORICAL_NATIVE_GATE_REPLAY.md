# Historical native consensus and gated rescue: corrected comparison

Active result paths and superseded-output cleanup are recorded in
[the authoritative result index](SEQC2_AUTHORITATIVE_RESULTS.md).

Audit date: 2026-09-14. Branch: `seqc2-consolidated`, starting commit `04bde89d`.

## Conclusion

The earlier manual WES improvement is reproducible. The recent three-policy
experiment did not test that method: it tested consensus-only CLI variants,
with different indel and SNP selection semantics, and no gated-rescue query.
Its negative conclusion must not be generalized to the manual combined method.

With the historical WGS candidate scope recovered, native consensus gains one
TP with no extra FP on UKB and ties DeepSomatic on MedExome. The historical
gate applied to the completed WGS hybrid realignment-rescue export adds three
FP and zero TP on UKB, including one FP on MedExome. The combined method
therefore does not outperform DeepSomatic on WGS in these comparisons.

No production policy is selected or changed by this audit. In particular,
reproducing a historical DeepSomatic-derived indel set is not acceptance of
that rule as the future consensus algorithm.

## Recovered provenance and discrepancies

1. The historical WES native metrics point to
   `/tmp/consensus_experiments/wesll/policy_relaxed/query.vcf.gz`, which still
   exists. The winning gated query is
   `examples/seqc2/hybrid/comparison/authoritative_medexome_bed_20260910/native_gated_rescue.vcf.gz`.
   Benchmarking these recovered alleles reproduces the MedExome winning counts.
2. Reconstructing the WES gate from that native baseline plus the historical
   second-rescue candidate VCF gives **exact allele-set equality** with the
   historical gated VCF. The gate is a union retaining every baseline allele;
   additions must be SNVs with `N_DNA_CALLERS_SUPPORT >= 1` and
   `N_RNA_CALLERS_SOMATIC >= 2`. The first field is eligible support, not a
   requirement that the DNA caller label itself be Somatic. These semantics
   must not be silently replaced by a stricter description.
3. The historical WGS evidence BED `/tmp/wgsil_sites.bed` contains 2,190
   positions. It is exactly the union of positions in its cached DeepSomatic
   PASS and original-consensus VCFs. The raw DeepSomatic and Mutect2 evidence
   snapshots were queried only at these positions. This candidate restriction
   was missing from the prose rule description.
4. Applying the reconstructed WGS thresholds to the **current hybrid caller
   VCFs**, restricted to that recovered candidate universe, exactly reproduces
   all 2,188 alleles of the historical WGS native query. There are zero added
   and zero missing alleles. This verifies the historical WGS result against
   current caller evidence, not merely its old metric JSON.
5. Applying those thresholds to all available caller sites admits 20 additional
   SNP alleles. This is a separate ablation, not an exact historical replay.
   Their DeepSomatic records were absent from the old position-restricted
   evidence snapshot, not necessarily absent from the full caller output.
6. The deployed `native_evidence_snv` override differs from the manual selection:
   indels remain threshold-based; qualified DeepSomatic SNVs bypass the two
   Mutect2 combination vetoes; and ordinary classification remains available
   when the native override fails. The manual WGS reconstruction applies its
   exact veto combinations to the selected baseline as well as additions.
   These differences must be resolved explicitly before claiming workflow
   parity with the manual experiment.

## Common benchmark contract

All new comparisons use the same merged truth VCF, SEQC2 HC BED, assembly38
FASTA, installed som.py, `-N`, and region-specific `-T`. No `-P` is used.
Selected alleles are written as PASS in separate, allele-only benchmark VCFs;
these are not training-label outputs. RNA candidates are selected from
Somatic/PASS exports. RaVeX_FILTER is not used as a selection criterion.

The historical native inputs were already preselected to HC and UKB. The
replay gives DeepSomatic the same preselection and then applies either UKB or
MedExome as `-T`; this preserves the historical domain restriction instead of
claiming unrestricted MedExome candidate discovery. DeepSomatic's counts
match the previous common matrix in all four cells. The recovered WES native
UKB query gives 1,051 TP under this normalization contract, rather than the
1,050 TP printed in one earlier table; the current replay value is used below.

The WGS gated query uses the mandatory completed second-rescue artifact:

```text
<shared repo>/examples/seqc2/hybrid/output/seqc2.wgs.il.hybrid/
  vcf_realignment/rescue/
  WGS_IL_T_1_vs_WGS_IL_N_1_rescued_WGS_IL_RT_1_realign_vs_WGS_IL_N_1/
  WGS_IL_T_1_vs_WGS_IL_N_1_rescued_WGS_IL_RT_1_realign_vs_WGS_IL_N_1.rescue.filtered.stripped.vep.vcf.gz
```

This is a VCF-only union/gate replay using the completed second-rescue
candidates. It is not a new full workflow execution or proof of downstream
annotation parity with a newly generated native baseline.

## Verified metrics

TP/FP/FN are shown separately for SNP, indel, and all records. P/R/F1 refers to
records. All per-type precision, recall, and F1 values are retained in the
machine-readable provenance and native som.py metrics.

| Dataset / target / method | SNP TP/FP/FN | Indel TP/FP/FN | Records TP/FP/FN | Records P / R / F1 |
| --- | ---: | ---: | ---: | ---: |
| WES-LL / UKB / DeepSomatic | 1007/34/1198 | 41/4/54 | 1048/38/1252 | 0.965009 / 0.455652 / 0.619019 |
| WES-LL / UKB / manual native | 1010/32/1195 | 41/4/54 | 1051/36/1249 | 0.966881 / 0.456957 / 0.620608 |
| WES-LL / UKB / manual native + gate | 1021/36/1184 | 41/4/54 | 1062/40/1238 | 0.963702 / 0.461739 / 0.624339 |
| WES-LL / MedExome / DeepSomatic | 540/17/255 | 23/2/11 | 563/19/266 | 0.967354 / 0.679131 / 0.798016 |
| WES-LL / MedExome / manual native | 541/15/254 | 23/2/11 | 564/17/265 | 0.970740 / 0.680338 / 0.800000 |
| WES-LL / MedExome / manual native + gate | 547/17/248 | 23/2/11 | 570/19/259 | 0.967742 / 0.687575 / 0.803949 |
| WGS-IL / UKB / DeepSomatic | 2083/8/122 | 85/11/10 | 2168/19/132 | 0.991312 / 0.942609 / 0.966347 |
| WGS-IL / UKB / manual native | 2084/8/121 | 85/11/10 | 2169/19/131 | 0.991316 / 0.943043 / 0.966578 |
| WGS-IL / UKB / manual native + gate | 2084/11/121 | 85/11/10 | 2169/22/131 | 0.989959 / 0.943043 / 0.965932 |
| WGS-IL / MedExome / DeepSomatic | 674/2/121 | 28/4/6 | 702/6/127 | 0.991525 / 0.846803 / 0.913468 |
| WGS-IL / MedExome / manual native | 674/2/121 | 28/4/6 | 702/6/127 | 0.991525 / 0.846803 / 0.913468 |
| WGS-IL / MedExome / manual native + gate | 674/3/121 | 28/4/6 | 702/7/127 | 0.990127 / 0.846803 / 0.912874 |

The WES UKB combined method improves F1 but has slightly lower precision than
DeepSomatic. MedExome improves both. A higher observed F1 on these assays is
not proof of general SOTA performance or independent biological validation.

The separate all-site WGS ablation gives 2,177/31/123 records before the gate
and 2,177/34/123 after it on UKB; MedExome gives 703/8/126 and 703/9/126.
Those results must not be substituted for the recovered candidate-scope replay.

## Reproduction and retained evidence

Run from the repository root into a fresh output directory:

```bash
.venv/bin/python examples/seqc2/scripts/replay_historical_native_gate.py --outdir examples/seqc2/comparison/historical_manual_replay_NEW
.venv/bin/python examples/seqc2/scripts/replay_historical_candidate_scope.py --replay-dir examples/seqc2/comparison/historical_manual_replay_NEW
```

These forensic replay scripts prefer the durable snapshots in
`examples/seqc2/verified/20260914/provenance_inputs/`, with the recovered `/tmp`
paths as a fallback. They are not portable production launchers.
Source SHA-256 identities, exact
commands, equality checks, metrics, logs, and derived queries are preserved in:

- `examples/seqc2/comparison/historical_manual_replay_20260914/provenance.json`
  contains WES parity and WGS DeepSomatic controls; all-site metrics are retired.
- `examples/seqc2/comparison/historical_manual_replay_20260914/wgs_historical_scope/provenance.json`
  contains the exact WGS historical-scope native replay and its gated rescue.
- `tests/test_historical_native_gate_replay.py` checks baseline retention,
  indel non-promotion, veto scope, positive-quality requirements, missing
  evidence, and exclusion of Artifact rescue records. Four additional
  benchmark-contract tests also passed (8 total).

Next policy work must start from these corrected comparisons: examine the
WGS rescue-only false positives and the historical candidate-universe rule,
then specify one explicit portable algorithm. Preserve the distinction
between experimental indel retention and an accepted evidence-based indel
consensus. No mapping, caller execution, input modification, cache change, or
HG008 tuning was performed.
