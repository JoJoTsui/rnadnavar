# Verified VCFs, source workflow outputs, and archive

Organized 2026-09-14. Optimization is paused. No workflow, consensus, rescue,
or benchmark was rerun as part of this organization.

## One entry point

From the current repository root, use:

```text
examples/seqc2/verified/20260914/
```

Each dataset contains these verified links:

```text
wes_ll/ or wgs_il/
  native_consensus.vcf.gz           verified manual native allele set
  native_gated_rescue.vcf.gz        verified native + second-rescue gate
  dna_callers/                      named DeepSomatic, Mutect2, Strelka VCFs
  rna_callers/                      first RNA caller VCFs
  rna_realign_callers/              realigned RNA caller VCFs
  alignments/
    dna_normal.bam + .bai           actual caller-ready reordered input
    dna_tumor.bam + .bai            actual caller-ready reordered input
    rna_tumor.cram + .crai          first RNA caller input
    rna_realign.cram + .crai        relevant second-rescue caller input
  original_inputs/                 original external DNA BAM/BAI files
```

The selected VCFs are PASS-labelled, allele-only benchmark queries. They are
the sets whose reported performance was verified. They are **not** newly
generated, fully annotated workflow consensus/rescue training-label VCFs.
Original workflow consensus/rescue exports remain separately linked under
`original_workflow_consensus/` and `original_workflow_rescue/`; do not substitute
them for the selected manual query merely because they carry richer annotations.

## Exact upstream source mapping

| Verified product | Actual workflow source | Link under the bundle |
| --- | --- | --- |
| WES-LL native DNA evidence | Current repo: `examples/seqc2/hybrid/output/seqc2.wes.ll.hybrid.realign.policy-default` | `wes_ll/workflow_dna_source` |
| WES-LL winning gated-rescue candidates | Current repo: historical `seqc2.wes.ll.hybrid.realign.latest`, now archived | `wes_ll/workflow_rescue_source` |
| WGS-IL DNA and gated-rescue evidence | Shared repo: `examples/seqc2/hybrid/output/seqc2.wgs.il.hybrid` | `wgs_il/workflow_source` |

Shared repo means
`/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar`.

The WES rescue source correction matters: historical `realign.latest` has an
exact 2,668-allele Somatic match to the cached winning experiment's rescue
candidate export. `policy-default` has only 2,633 of those alleles. Their
realigned RNA CRAMs are different. The selected WES RNA-realignment link is
therefore the `fe/475e12...` task CRAM, not the later `29/b3b340...` CRAM.
Their DNA caller tasks are shared.

The task traces bind the actual alignment files to all three DNA/RNA callers.
`OUTPUT_MANIFEST.json` contains resolved paths, staged inputs, separately
resolved index files, task hashes, executed caller commands, source VCF hashes,
archive inventories, and the alignment quickcheck results. A CRAM index may
come from a separate indexing task; the bundle uses the exact index staged in
the caller task rather than assuming it sits beside the original CRAM.

## Archive, not deletion

Thirteen old publication directories were moved into archives: six earlier
WES hybrid iterations in each repo, plus the current repo's non-authoritative
WGS hybrid publication copy. The archived WES names are:

- `seqc2.wes.ll.hybrid`
- `seqc2.wes.ll.hybrid.modality`
- `seqc2.wes.ll.hybrid.pooling_fix`
- `seqc2.wes.ll.hybrid.realign.full`
- `seqc2.wes.ll.hybrid.realign.latest`
- `seqc2.wes.ll.hybrid.realign.preview`

Archive locations:

```text
<current repo>/.artifacts/archived_seqc2_workflow_outputs_20260914/
<shared repo>/.artifacts/archived_seqc2_workflow_outputs_20260914/
```

Old output paths are compatibility symlinks to their archives. This preserves
existing absolute-path provenance and downstream readers. Those old names
are not an active-results registry; use the verified bundle and its manifest.
Every archived entry retained its device/inode and size. Nothing was deleted,
compressed, or copied over another output; archiving does not free disk space.

Retained in place: the source WES `policy-default` publications, authoritative
shared WGS hybrid output, independent DNA-only WES-LL/WES-IL/WGS-IL datasets,
all work/Conda/Nextflow caches, and external original inputs. The running HG008
job and its outputs remain untouched. A dataset is not obsolete merely because
it has lower coverage or a lower benchmark score.

## Durable evidence and acceptance checks

Historical native VCFs, the original WGS consensus and DeepSomatic PASS VCF,
and the candidate-position BED were copied from `/tmp` into
`provenance_inputs/`, with checksums. Replay scripts now prefer those durable
copies. Source copies retain their original contents; no annotations were
invented or rewritten.

Verified checks include all archive inventories, all curated links and indexes,
source VCF checksums, the 12 authoritative metric hashes, and `samtools
quickcheck` on nine distinct task-input BAM/CRAM files across the retained
lineages. These are integrity/provenance checks, not claims of read-level
biological validation.

The current metric index remains
[SEQC2_AUTHORITATIVE_RESULTS.md](SEQC2_AUTHORITATIVE_RESULTS.md), also linked
from the bundle's `comparison/`. The bundle's native/gated VCF links point to
those same verified queries. Any future rule optimization remains a separate
task after this organization.
