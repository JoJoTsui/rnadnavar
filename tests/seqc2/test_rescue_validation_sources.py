from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'examples/seqc2/scripts'))
from validate_refined_rescue_output import rescue_pattern


def test_round_specific_rescue_source_patterns():
    assert rescue_pattern('wgs', 'first') == '*/*.filtered.vcf.gz'
    assert rescue_pattern('wgs', 'realignment') == '*/*.rescue.filtered.stripped.vep.vcf.gz'
    assert rescue_pattern('wes', 'first') == rescue_pattern('wes', 'realignment') == '*/*.filtered.vcf.gz'
