"""Read-only characterization of the collector, not training-policy approval."""
from pathlib import Path
import sys

import pysam

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'examples/seqc2/scripts'))
from pilot_three_class_bam_evidence import evidence


def collect(tmp_path, descriptions):
    sequence = 'ACGT' * 250
    fasta = tmp_path/'ref.fa'
    fasta.write_text('>chr1\n'+sequence+'\n')
    pysam.faidx(str(fasta))
    path = tmp_path/'reads.bam'
    with pysam.AlignmentFile(str(path), 'wb', header={
            'HD':{'VN':'1.6','SO':'coordinate'}, 'SQ':[{'SN':'chr1','LN':1000}]}) as bam:
        for name,flag,mq,bq,start,mate in sorted(descriptions, key=lambda x:x[4]):
            read = pysam.AlignedSegment(bam.header)
            read.query_name=name
            read.flag=flag
            read.reference_id=0
            read.reference_start=start
            read.mapping_quality=mq
            read.cigarstring='100M'
            read.query_sequence=sequence[start:start+100]
            read.query_qualities=[bq]*100
            if flag & 1:
                read.next_reference_id=0
                read.next_reference_start=mate
                read.template_length=120 if start < mate else -120
            bam.write(read)
    pysam.index(str(path))
    with pysam.AlignmentFile(str(path)) as bam, pysam.FastaFile(str(fasta)) as ref:
        return evidence(bam, ref, ('chr1',151,'G','T'))


def test_duplicate_secondary_qcfail_low_quality_are_excluded(tmp_path):
    reads=[('good',0,60,40,100,0)]
    reads += [(str(flag),flag,60,40,100,0) for flag in (256,512,1024)]
    reads += [('low_mapq',0,10,40,100,0), ('low_bq',0,60,10,100,0)]
    assert collect(tmp_path, reads) == dict(depth=1,ref=1,alt=0,other=0)


def test_overlapping_pair_is_one_observation(tmp_path):
    reads=[('pair',99,60,40,100,120), ('pair',147,60,40,120,100)]
    assert collect(tmp_path, reads) == dict(depth=1,ref=1,alt=0,other=0)


def test_characterize_supplementary_alignment_independence_gap(tmp_path):
    # The current default filter does not exclude supplementary alignments.
    # Repeated sequence from one molecule can therefore contribute twice.
    # This passing characterization is a documented safety finding, NOT proof
    # that two such observations satisfy the binomial independence assumption.
    reads=[('one_molecule',0,60,40,100,0), ('one_molecule',2048,60,40,100,0)]
    assert collect(tmp_path, reads) == dict(depth=2,ref=2,alt=0,other=0)
