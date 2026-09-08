"""Production annotation routing with fake executables; no annotation accuracy claim.

RUN_NEXTFLOW_VEP_TESTS=1 .venv/bin/python -m pytest tests/seqc2/test_vep_routing.py -q
"""
import os
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.skipif(
    os.environ.get('RUN_NEXTFLOW_VEP_TESTS') != '1', reason='opt-in Nextflow integration'
)


@pytest.mark.parametrize('mode', ['hybrid', 'seq2neo', 'disabled', 'missing_genome', 'missing_species', 'missing_version'])
def test_vep_routing(tmp_path, mode):
    cache = tmp_path / 'cache'
    cache.mkdir()
    fasta = tmp_path / 'reference.fa'
    fasta.write_text('>chr1\nA\n')
    vcf = tmp_path / 'input.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tC\t.\tSomatic\t.\n')
    binaries = tmp_path / 'bin'
    binaries.mkdir()
    scripts = {
        'vep': '''#!/usr/bin/env python3
import gzip, pathlib, sys
args = sys.argv[1:]
if '--help' in args:
    print('ensembl-vep : 115')
else:
    with gzip.open(args[args.index('-o') + 1], 'wb') as out:
        out.write(pathlib.Path(args[args.index('-i') + 1]).read_bytes())
''',
        'tabix': '''#!/usr/bin/env python3
import pathlib, sys
if '-h' in sys.argv:
    print('Version: 1.0')
else:
    pathlib.Path(sys.argv[-1] + '.tbi').touch()
''',
    }
    for name, script in scripts.items():
        path = binaries / name
        path.write_text(script)
        path.chmod(0o755)
    tools = 'rescue,vep' if mode == 'seq2neo' else 'rescue,realignment'
    genome = 'null' if mode in ('missing_genome', 'disabled') else "'GRCh38'"
    species = 'null' if mode == 'missing_species' else "'homo_sapiens'"
    version = 'null' if mode == 'missing_version' else '115'
    config = tmp_path / 'test.config'
    config.write_text(f"""
params {{
    input = '{vcf}'
    fasta = '{fasta}'
    vep_cache = '{cache}'
    outdir = '{tmp_path / 'out'}'
    step = 'mapping'
    tools = '{tools}'
    test_realignment = {str(mode != 'disabled').lower()}
    vep_include_fasta = true
    vep_genome = {genome}
    vep_species = {species}
    vep_cache_version = {version}
    dbnsfp = null
    spliceai_snv = null
}}
process.executor = 'local'
process.ext.args = '--vcf'
process.ext.when = true
conda.enabled = false
docker.enabled = false
singularity.enabled = false
""")
    result = subprocess.run(
        ['nextflow', '-C', str(config), 'run', str(ROOT / 'tests/seqc2/vep_routing.nf'),
         '-ansi-log', 'false', '-with-trace', str(tmp_path / 'trace.tsv')],
        cwd=tmp_path, env={**os.environ, 'NXF_OFFLINE': 'true', 'PATH': f"{binaries}:{os.environ['PATH']}"},
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120,
    )
    if mode.startswith('missing_'):
        parameter = {'missing_genome': 'vep_genome', 'missing_species': 'vep_species',
                     'missing_version': 'vep_cache_version'}[mode]
        assert result.returncode != 0, result.stdout
        assert f'--{parameter}' in result.stdout, result.stdout
    else:
        assert result.returncode == 0, result.stdout
        trace = (tmp_path / 'trace.tsv').read_text()
        assert trace.count('ENSEMBLVEP_VEP') == {'hybrid': 2, 'seq2neo': 4, 'disabled': 0}[mode], result.stdout
        assert 'ANNOTATED_SAMPLES=[RT1, RT2]' in result.stdout, result.stdout
