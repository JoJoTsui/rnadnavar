"""Truth-confidence partition and veto attribution for the forensic assay."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'examples/seqc2/scripts'))
from audit_rescue_gate_domain import domain, reasons
from test_rescue_nomination_gate import biological_veto


def test_target_membership_does_not_confer_truth_confidence():
    key = ('chr1', 10, 'A', 'C')
    assert domain(key, set(), {key}, {key}) == 'outside_hc_unassessed'
    assert domain(key, {key}, set(), {key}) == 'hc_outside_ukb_medexome'
    assert domain(key, {key}, {key}, set()) == 'hc_ukb_outside_medexome'


def test_reason_attribution_matches_fixed_experiment_and_retains_overlap():
    examples = [{}, {'GNOMAD_AF': '.'}, {'GNOMAD_AF': '0.001'},
                {'GNOMAD_AF': '0.02'},
                {'REDI_ACCESSION': 'site', 'REDI_CANONICAL': 'YES', 'N_DNA_CALLERS_SOMATIC': '0'},
                {'REDI_ACCESSION': 'site', 'REDI_CANONICAL': 'YES', 'N_DNA_CALLERS_SOMATIC': '1'}]
    for info in examples:
        rows = [{'info': info}]
        assert bool(reasons(rows)) == biological_veto(rows)
    assert reasons([{'info': examples[3]}, {'info': examples[4]}]) == {
        'common_af', 'editing_without_dna_somatic'}
