"""Validate the independent real-BAM indel audit against known synthetic reads."""
import pysam

import export_variant_evidence as e
import verify_variant_evidence as v


def test_independent_indel_audit(tmp_path):
    path = tmp_path / 'events.bam'
    with pysam.AlignmentFile(str(path), 'wb', header={'HD': {'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 1000}]}) as bam:
        for name, seq, cigar in [('ref', 'AAAAAA', '6M'), ('ins', 'AACTAAAA', '2M2I4M'),
                                 ('del', 'AAAA', '2M2D2M'), ('splice', 'AAAA', '2M2N2M'),
                                 ('anchor_only', 'AA', '2M')]:
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = seq
            read.query_qualities = [30] * len(seq)
            read.reference_id = 0
            read.reference_start = 0
            read.cigarstring = cigar
            read.mapping_quality = 60
            bam.write(read)
    pysam.index(str(path))
    with pysam.AlignmentFile(str(path), 'rb') as bam:
        insertion = v.independent_indel_counts({'CHROM': 'chr1', 'POS': 2, 'REF': 'A', 'ALT': 'ACT'}, bam, False)
        deletion = v.independent_indel_counts({'CHROM': 'chr1', 'POS': 2, 'REF': 'AAA', 'ALT': 'A'}, bam, False)
    assert insertion == {'depth': 2, 'ref_count': 1, 'alt_count': 1}
    assert deletion == {'depth': 3, 'ref_count': 1, 'alt_count': 1}


def test_invalid_numeric_measurement_is_rejected():
    import polars as pl
    import pytest
    row = {}
    for tag in e.MODALITIES:
        row.update({f'bam_{tag}_{k}': val for k, val in e.unavailable('bam_absent').items()})
    row['bam_DN_alt_count'] = 0
    with pytest.raises(ValueError, match='fabricated'):
        e.verify_measurements(pl.DataFrame([row]))
