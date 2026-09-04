# Converge hybrid inputs at the caller-ready alignment boundary

Status: accepted

The pipeline must reuse caller-ready DN/DT alignments while aligning raw RT reads in one run, but its global `step` currently either remaps every BAM or rejects FASTQ input. Hybrid sheets therefore declare an all-or-none per-row input stage: raw reads follow the existing mapping and preprocessing route, caller-ready BAMs follow the existing audit and safe-normalization route, and both converge as the same caller-ready alignment contract before pairing and variant calling. Legacy sheets without input stages retain their current behavior; invalid caller-ready inputs fail rather than silently falling back to remapping, and ingress provenance remains separate from downstream metadata so seq2neo naming, joins, and output contracts stay stable.

## Consequences

RNA library repeats retain distinct library identities through read groups and are pooled only as one logical RT sample. Hybrid manifests are validated completely before tasks start, external CRAMs remain strict-only, and the existing global interval and caller configuration continue to apply. The first SEQC2 hybrid run is a cell-line-matched engineering benchmark rather than specimen-matched training truth, and its DNA calls must remain content-equivalent to the existing DNA-only baseline.
