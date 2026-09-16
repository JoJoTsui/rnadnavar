"""Standalone validation must keep selection separate from truth/target scoring."""
import importlib.util
import json
from pathlib import Path
import sys

import pysam
import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / 'examples/seqc2/scripts'
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location('frozen_validation', SCRIPTS / 'validate_frozen_hybrid_policy.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_query_is_streamed_sampleless_pass_copy(tmp_path):
    source = tmp_path / 'source.vcf'
    source.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n'
                      '##FILTER=<ID=Somatic,Description="label">\n'
                      '##FILTER=<ID=Artifact,Description="label">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                      'chr1\t1\t.\tA\tC\t.\tSomatic\t.\n'
                      'chr1\t2\t.\tA\tG\t.\tArtifact\t.\n'
                      'chr1\t3\t.\tA\tT\t.\tPASS\t.\n'
                      'chr1\t4\t.\tA\tC\t.\t.\t.\n')
    before = source.read_bytes()
    for biological, expected in ((True, [1]), (False, [3, 4])):
        dest = tmp_path / f'query.{biological}.vcf.gz'
        module.pass_query(source, dest, biological)
        with pysam.VariantFile(str(dest)) as reader:
            assert not reader.header.samples
            rows = list(reader)
            assert [r.pos for r in rows] == expected
            assert all(set(r.filter) == {'PASS'} for r in rows)
        assert Path(str(dest) + '.tbi').is_file()
    assert source.read_bytes() == before


def test_manifest_requires_both_rounds_and_existing_distinct_callers(tmp_path):
    paths = []
    for i in range(3):
        p = tmp_path / str(i)
        p.touch()
        paths.append(str(p))
    panel = dict(zip(module.CALLERS, paths))
    data = dict(dna=panel, rna_first=panel, rna_realignment=panel,
                rescues=dict(first=paths[0], realignment=paths[1]),
                targets=dict(ukb=paths[0], medexome=paths[1]),
                truth=paths[0], hc=paths[1], fasta=paths[2])
    path = tmp_path / 'manifest.json'
    path.write_text(json.dumps(data))
    assert module.load_manifest(path)[0] == data
    data['rescues'].pop('realignment')
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='both rescue'):
        module.load_manifest(path)


def test_no_truth_region_filter_before_scoring():
    source = (SCRIPTS / 'validate_frozen_hybrid_policy.py').read_text()
    normalization = source.split("run(['bcftools', 'norm'", 1)[1].split('normalize.{caller}.log', 1)[0]
    assert "'-R'" not in normalization and "'-T'" not in normalization
    assert '--experimental-refined-native' in source
    assert "'nextflow'" not in source
