"""Guard evidence roles and read/template counting in SEQC2-only assays."""
from pathlib import Path
import sys

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'examples/seqc2/scripts'))
from explore_indel_corroboration import gates
from audit_retained_rescue_templates import summarize_reads
from refine_indel_corroboration import refined_gates


def evidence():
    ds = {'filter': 'RefCall', 'qual': '1', 'samples': {'WES_LL_T_1': {'AD': '7,3'}}}
    m2 = {'filter': 'PASS', 'info': 'GERMQ=20;TLOD=5', 'samples': {
        'WES_LL_T_1': {'AD': '7,3', 'SB': '3,4,1,2'},
        'WES_LL_N_1': {'AD': '10,0', 'DP': '10'}}}
    st = {'filter': 'PASS', 'samples': {
        'TUMOR': {'TIR': '3,4'}, 'NORMAL': {'TIR': '0,0', 'DP': '10'}}}
    return ds, m2, st


def test_reciprocal_support_and_strand_ablation():
    ds, m2, st = evidence()
    assert gates(ds, m2, st) == {'mutect_pass', 'strelka_pass'}
    ds['filter'] = 'PASS'
    assert gates(ds, m2, st) == {'mutect_pass', 'strelka_pass', 'ds_strand'}
    m2['samples']['WES_LL_T_1']['SB'] = '3,4,0,3'
    assert 'ds_strand' not in gates(ds, m2, st)


def test_normal_evidence_and_artifact_rejection():
    ds, m2, st = evidence()
    m2['samples']['WES_LL_N_1']['AD'] = '9,1'
    assert not gates(ds, m2, st)
    ds, m2, st = evidence()
    m2['filter'] = 'normal_artifact'
    assert not gates(ds, m2, st)
    ds, m2, st = evidence()
    ds['filter'] = 'GERMLINE'
    assert 'mutect_pass' not in gates(ds, m2, {})


def test_refinement_event_scope_and_moderate_confidence():
    ds, m2, _ = evidence()
    m2['info'] += ';ECNT=1'
    assert refined_gates(ds, m2) == {'single_event'}
    m2['info'] = m2['info'].replace('ECNT=1', 'ECNT=2')
    assert not refined_gates(ds, m2)
    ds['filter'], ds['qual'] = 'PASS', '10'
    m2['samples']['WES_LL_T_1']['AD'] = '8,2'
    assert refined_gates(ds, m2) == {'moderate_strand'}
    ds['qual'] = '9.9'
    assert not refined_gates(ds, m2)
    ds['qual'] = '10'
    m2['samples']['WES_LL_T_1']['SB'] = '3,4,0,2'
    assert not refined_gates(ds, m2)


def read(name, sequence, flag):
    r = pysam.AlignedSegment()
    r.query_name, r.query_sequence, r.flag = name, sequence, flag
    r.reference_id, r.reference_start, r.mapping_quality = 0, 9, 60
    r.cigarstring = '3M'
    r.query_qualities = pysam.qualitystring_to_array('III')
    return r


def test_mates_count_as_one_template_and_conflicts_are_not_alt_support():
    result = summarize_reads([read('same', 'ACA', 65), read('same', 'ACA', 129)], 11, 'A', 'C')
    assert result['counts']['alt_reads'] == 2
    assert result['counts']['alt_templates'] == 1
    result = summarize_reads([read('same', 'ACA', 65), read('same', 'AAA', 129),
                              read('duplicate', 'ACA', 1024)], 11, 'A', 'C')
    assert result['counts']['alt_templates'] == 0
    assert result['counts']['discordant_templates'] == 1
    assert result['counts']['excluded_flags'] == 1
