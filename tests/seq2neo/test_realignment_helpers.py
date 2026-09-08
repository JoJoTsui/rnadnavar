"""Run real helper modules and the production reference-selection expression."""
import json
import pytest
from test_alignment_pooling import ROOT, run_nextflow


@pytest.mark.parametrize('module,process,payload,output', [
    ('validate_read_ids', 'VALIDATE_READ_IDS', 'readA\n', 'read_ids'),
    ('sort_merge_bed', 'SORT_MERGE_BED', 'chr1\t20\t30\nchr1\t10\t25\n', 'bed'),
    ('filter_hisat_splicesites', 'FILTER_HISAT_SPLICESITES', 'chr1\t10\t20\t+\n', 'splicesites'),
])
def test_helper_executes_and_reports_versions(tmp_path, module, process, payload, output):
    source = tmp_path / 'input.txt'
    source.write_text(payload)
    fai = tmp_path / 'ref.fai'
    fai.write_text('chr1\t100\t0\t100\t101\n')
    extra = f", file('{fai}')" if process == 'FILTER_HISAT_SPLICESITES' else ''
    result = run_nextflow(tmp_path, f"""
include {{ {process} }} from '{ROOT}/modules/local/{module}/main'
workflow {{
    {process}(Channel.of([[id:'sample'], file('{source}'){extra}]))
    {process}.out.versions.view {{ 'VERSIONS=' + it.text }}
    {process}.out.{output}.view {{ meta, path -> 'DATA=' + groovy.json.JsonOutput.toJson(path.text) }}
}}
""", [])
    assert result.returncode == 0, result.stdout + result.stderr
    assert f'"{process}":' in result.stdout
    expected = 'chr1\t10\t30\n' if process == 'SORT_MERGE_BED' else payload
    assert 'DATA=' + json.dumps(expected, separators=(',', ':')) in result.stdout


@pytest.mark.parametrize('supplied', [True, False])
def test_splice_filter_receives_supplied_or_generated_fai(tmp_path, supplied):
    source = (ROOT / 'subworkflows/local/prepare_genome/main.nf').read_text()
    selection = next((line.strip() for line in source.splitlines()
                      if line.strip().startswith('splice_fai =')), '')
    invocation = next(line.strip() for line in source.splitlines()
                      if 'FILTER_HISAT_SPLICESITES(supplied_splicesites.combine(' in line)
    expression = invocation[len('FILTER_HISAT_SPLICESITES('):-1]
    fai = tmp_path / 'ref.fai'
    fai.write_text('chr1\t100\t0\t100\t101\n')
    configured = f"'{fai}'" if supplied else 'null'
    generated = 'Channel.empty()' if supplied else f"Channel.of([[id:'ref'], file('{fai}')])"
    result = run_nextflow(tmp_path, f"""
params.fasta_fai = {configured}
workflow {{
    SAMTOOLS_FAIDX = [out: [fai: {generated}]]
    supplied_splicesites = Channel.of([[id:'splice'], 'sites.txt'])
    {selection}
    ({expression}).toList().view {{ rows -> 'COUNT=' + rows.size() }}
}}
""", [])
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'COUNT=1' in result.stdout, result.stdout


def test_conversion_sizes_follow_staged_symlinks(tmp_path):
    """Exercise the production shell assignment with Nextflow-style staging."""
    import subprocess
    source = (ROOT / 'modules/local/samtools_convert_enhanced/main.nf').read_text()
    assignment = next(line.strip() for line in source.splitlines()
                      if line.strip().startswith('input_size='))
    target = tmp_path / 'original.cram'
    target.write_bytes(b'x' * 4096)
    staged = tmp_path / 'staged.cram'
    staged.symlink_to(target)
    command = assignment.replace('${input}', str(staged)).replace('\\$', '$')
    result = subprocess.run(['bash', '-eu', '-c', command + '\nprintf "%s" "$input_size"'],
                            capture_output=True, text=True, check=True)
    assert int(result.stdout) == target.stat().st_size


@pytest.mark.parametrize('extra', [
    {},
    {'strandedness': 'reverse', 'library': 'RT_lib', 'input_stage': 'raw_reads'},
    {'strandedness': 'forward', 'libraries': ['RT_A', 'RT_B']},
])
def test_realignment_conversion_preserves_hisat_metadata(tmp_path, extra):
    source = (ROOT / 'subworkflows/local/enhanced_cram2bam_conversion/main.nf').read_text()
    mapping = source.split('.map { meta, cram, crai ->', 1)[1].split(
        '// Skip reference file validation', 1)[0]
    meta = dict(id='RT', patient='P', sample='RT', status=2, single_end=False,
                data_type='cram', readsid_path='reads.txt', **extra)
    result = run_nextflow(tmp_path, '''
workflow {
    def rows = new groovy.json.JsonSlurper().parseText(params.rows)
    Channel.fromList(rows).map { meta, cram, crai ->
''' + mapping + '''
    .view { meta, cram, crai -> 'META=' + groovy.json.JsonOutput.toJson(meta) }
}
''', [[meta, 'RT.cram', 'RT.cram.crai']])
    assert result.returncode == 0, result.stdout + result.stderr
    observed = json.loads(next(line[5:] for line in result.stdout.splitlines()
                               if line.startswith('META=')))
    for key, value in meta.items():
        assert observed.get(key) == value, (key, observed)
