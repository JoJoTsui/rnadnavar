# Validate the frozen rescue experiment before claiming default-policy parity

Status: accepted validation boundary, 2026-09-14; detailed validation design remains under discussion.

The successful SEQC2 standalone nomination-gate experiment is a distinct policy
from the similarly named native-consensus and rescue flags in the workflow.
We will validate the frozen experiment and audit production parity separately,
without changing global defaults before independent validation as required by
[ADR-0007](0007-default-native-consensus-and-gated-rescue.md). Reproducing a
benchmark score alone cannot establish correct training labels: allele selection,
FILTER, sample identity, annotations, evidence provenance, and completed
realignment-rescue lineage must also agree.

Validation may read existing WES-LL/WGS-IL BAMs/CRAMs and execute isolated
consensus/rescue/annotation checks into fresh output directories. Mapping,
variant calling, cohort reruns, HG008 interference, and modification of source
outputs or caches are out of scope. BAM, hybrid, and FASTQ-triplet routing will
be checked through configuration and focused tests. SNP/indel/overall scores
must reproduce in both target regions, with no TP loss or FP increase relative
to the frozen candidate during integration; DeepSomatic remains a separate
comparator, ties are acceptable, and historical indel retention stays provisional.

This deliberately separates implementation correctness from independent
biological validation, avoiding the irreversible propagation of an unvalidated
policy into cohort training labels. See the [validation interview record](../review/SEQC2_UPDATED_POLICY_VALIDATION_INTERVIEW.md).
