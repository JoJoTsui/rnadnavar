"""Opt-in real Nextflow regression: RUN_NEXTFLOW_RESOURCE_TESTS=1.

Requires nextflow (or micromamba's nextflow environment), samtools, hisat2,
and hisat2_extract_splice_sites.py on PATH. No human data or full indices used.
"""
import os
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.skipif(
    os.environ.get('RUN_NEXTFLOW_RESOURCE_TESTS') != '1',
    reason='opt-in Nextflow resource integration tests',
)


@pytest.mark.parametrize('mode', ['supplied', 'generated_fai', 'generated_index', 'seq2neo', 'disabled', 'missing_fai'])
def test_reference_resource_routing(tmp_path, mode):
    fasta = tmp_path / 'tiny.fa'
    fasta.write_text('>chr1\n' + 'A' * 100 + '\n')
    fai = tmp_path / 'tiny.fa.fai'
    fai.write_text('chr1\t100\t6\t100\t101\n')
    splice = tmp_path / 'input.splice_sites.txt'
    splice.write_text('chr1\t9\t39\t+\nunknown\t9\t39\t+\n')
    gtf = tmp_path / 'tiny.gtf'
    gtf.write_text(
        'chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
        'chr1\ttest\texon\t40\t60\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
    )
    # The consumer only stages the resource; index contents are not aligned here.
    index = tmp_path / 'index'
    index.mkdir()
    (index / 'tiny.1.ht2').write_text('fixture\n')
    fai_param = 'null' if mode == 'generated_fai' else f"'{fai}'"
    if mode == 'missing_fai':
        fai_param = f"'{tmp_path / 'missing.fai'}'"
    splice_param = 'null' if mode == 'seq2neo' else f"'{splice}'"
    index_param = 'null' if mode == 'generated_index' else f"'{index}'"
    tools = 'consensus' if mode == 'disabled' else 'consensus,realignment'
    config = tmp_path / 'test.config'
    config.write_text(f"""
params {{
    fasta = '{fasta}'
    fasta_fai = {fai_param}
    splicesites = {splice_param}
    gtf = '{gtf}'
    gff = null
    dna = true
    rna = true
    read_length = 76
    step = 'mapping'
    tools = '{tools}'
    bwa = '{index}'
    aligner = 'bwa-mem'
    dict = 'supplied.dict'
    star_index = '{index}'
    hisat2_index = {index_param}
    hisat2_build_memory = '1 GB'
    save_reference = false
    build_only_index = false
    publish_dir_mode = 'copy'
    outdir = '{tmp_path / 'out'}'
}}
includeConfig '{ROOT / 'conf/modules/prepare_resources/prepare_genome.config'}'
conda.enabled = false
docker.enabled = false
singularity.enabled = false
process.executor = 'local'
""")
    command = ['nextflow'] if shutil.which('nextflow') else ['micromamba', 'run', '-n', 'nextflow', 'nextflow']
    result = subprocess.run(
        command + ['-C', str(config), 'run', str(ROOT / 'tests/seqc2/resource_routing.nf'),
                   '-ansi-log', 'false', '-with-trace', str(tmp_path / 'trace.tsv')],
        cwd=tmp_path, env={**os.environ, 'NXF_OFFLINE': 'true'},
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180,
    )
    if mode == 'missing_fai':
        assert result.returncode != 0, result.stdout
        assert 'missing.fai' in result.stdout, result.stdout
        return
    assert result.returncode == 0, result.stdout
    trace = (tmp_path / 'trace.tsv').read_text()
    if mode == 'disabled':
        assert 'USE_RESOURCES' not in trace
        assert 'HISAT' not in trace
    else:
        assert 'RESOURCE_SAMPLES=[RT1, RT2]' in result.stdout, result.stdout
        assert trace.count('USE_RESOURCES') == 2
        assert ('SAMTOOLS_FAIDX' in trace) == (mode == 'generated_fai')
        assert ('FILTER_HISAT_SPLICESITES' in trace) == (mode != 'seq2neo')
        assert ('HISAT2_EXTRACTSPLICESITES' in trace) == (mode == 'seq2neo')
        assert ('HISAT2_BUILD' in trace) == (mode == 'generated_index')
    assert fasta.read_text() == '>chr1\n' + 'A' * 100 + '\n'
    assert splice.read_text() == 'chr1\t9\t39\t+\nunknown\t9\t39\t+\n'
