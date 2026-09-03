# Normalize external BAMs at the shared caller boundary

Status: accepted

## Context

External BAMs can carry a sequence dictionary that differs from the selected
reference even when their shared human contigs use the same assembly.  The
SEQC2 WES BAMs contain 199 appended viral contigs and omit reference-only
alternate/HLA contigs.  GATK can be told to ignore this, but Strelka2 requires
exact dictionary equality.  Caller-specific workarounds would give the three
callers different evidence and therefore undermine consensus-label provenance.

FASTQ-derived BAMs are created against the selected reference and must not be
routed through an external-input repair path.

## Decision

At `--step variant_calling`, every external BAM is audited once before any
caller.  `bam_dictionary_policy=normalize` is the default: exact BAMs pass
through, while safely resolvable name/order differences are normalized with
Picard ReorderSam.  The normalized tumor/normal pair is the sole input pair for
Mutect2, Strelka2, and DeepSomatic.  `bam_dictionary_policy=strict` validates
but never rewrites.  The legacy `reorder_bam_contigs` boolean maps to these two
policies while deprecated.

Safe normalization requires every shared contig name to have the same length.
Extra BAM contigs must be absent from the selected reference; missing
reference contigs may be added to the output dictionary.  Duplicate/malformed
records, same-name length conflicts, and assembly incompatibility are fatal.
Requested intervals remain reference-scoped, so an off-reference contig
cannot be a valid requested interval.

External CRAMs are validated strictly and are not rewritten.  Mapping-derived
BAMs bypass this boundary entirely.

Each audit records dictionary differences, off-reference read counts, input
and reference checksums, the selected policy, and the decision.  The normalized
alignment is checked again for exact name/length/order equality.  Normalized
CRAM and audit artifacts are retained; the intermediate BAM is optional.

## Consequences

Normalizing an external BAM changes the caller input checksum and intentionally
invalidates all three caller caches.  It preserves records from off-reference
contigs as unmapped and reports the affected counts explicitly.  It allows strict validation to remain
enabled in every caller and prevents caller-specific reference views.
