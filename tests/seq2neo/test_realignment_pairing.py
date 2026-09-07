"""Exercise production realignment joins with actual Nextflow channels."""
import json

import pytest

from test_alignment_pooling import ROOT, run_nextflow


def pair(tmp_path, vcfs, crams):
    source = (ROOT / 'subworkflows/local/prepare_realignment_vcf/main.nf').read_text()
    operators = source.split('// === STEP 3: CHANNEL JOIN ===', 1)[1].split(
        '// === STEP 4: VCF TO BED CONVERSION ===', 1)[0]
    return run_nextflow(tmp_path, '''
workflow {
    def rows = new JsonSlurper().parseText(params.rows)
    vcf_to_realign = [realign: Channel.fromList(rows[0])]
    reads_to_realign_branch = [realign: Channel.fromList(rows[1])]
''' + operators + '''
    joined_data.toList().view { 'PAIRS=' + JsonOutput.toJson(it) }
}
''', [vcfs, crams])


def inputs(patient='WES_LL', sample='WES_LL_RT_1'):
    # Active caller outputs deliberately carry pair ID, patient and status only.
    vcf = [dict(patient=patient, id=f'{sample}_vs_DN', status=2), 'rt.vcf.gz', 'rt.vcf.gz.tbi']
    cram = [dict(patient=patient, sample=sample, id=sample, status=2), 'rt.cram', 'rt.cram.crai']
    return vcf, cram


@pytest.mark.parametrize('sample', ['WES_LL_RT_1', 'seq2neo_RT'])
def test_pair_id_metadata_matches_logical_rna_sample(tmp_path, sample):
    vcf, cram = inputs(sample=sample)
    result = pair(tmp_path, [vcf], [cram])
    assert result.returncode == 0, result.stdout + result.stderr
    rows = json.loads(next(line[6:] for line in result.stdout.splitlines() if line.startswith('PAIRS=')))
    assert len(rows) == 1
    assert rows[0][0]['sample'] == sample
    assert rows[0][0]['vcf_sample'] == sample


@pytest.mark.parametrize('missing', ['vcf', 'cram'])
def test_missing_partner_fails_clearly_not_as_zero_candidates(tmp_path, missing):
    vcf, cram = inputs()
    result = pair(tmp_path, [] if missing == 'vcf' else [vcf], [] if missing == 'cram' else [cram])
    assert result.returncode != 0
    output = result.stdout + result.stderr
    assert 'Join mismatch' in output
    assert 'Invalid method invocation' not in output


def test_duplicate_pair_fails(tmp_path):
    vcf, cram = inputs()
    result = pair(tmp_path, [vcf, vcf], [cram])
    assert result.returncode != 0
    assert 'duplicate' in (result.stdout + result.stderr).lower()


def test_patients_cannot_cross_pair(tmp_path):
    vcf, cram = inputs()
    cram[0]['patient'] = 'OTHER'
    result = pair(tmp_path, [vcf], [cram])
    assert result.returncode != 0
    assert 'Join mismatch' in result.stdout + result.stderr


@pytest.mark.parametrize('empty', [False, True])
def test_real_vcf2bed_preserves_all_classes_and_stops_empty_candidates(tmp_path, empty):
    vcf = tmp_path / 'candidate.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr7,length=1000>\n'
                   '##FILTER=<ID=Somatic,Description="Somatic">\n'
                   '##FILTER=<ID=Artifact,Description="Artifact">\n'
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                   + ('' if empty else 'chr7\t10\t.\tA\tT\t.\tSomatic\t.\nchr7\t20\t.\tC\tG\t.\tArtifact\t.\n'))
    tbi = tmp_path / 'candidate.vcf.tbi'
    tbi.touch()  # VCF2BED performs a sequential read, not an indexed query.
    source = (ROOT / 'subworkflows/local/prepare_realignment_vcf/main.nf').read_text()
    gate = next(line.strip() for line in source.splitlines() if 'bed = VCF2BED.out.bed.filter' in line)
    # Exercise the real qualified process publisher, not a module-only default.
    (tmp_path / 'nextflow.config').write_text(
        f"params.outdir = '{tmp_path}/output'\n"
        "params.publish_dir_mode = 'copy'\n"
        "params.save_output_as_bam = false\n"
        "params.save_align_intermeds = false\n"
        "params.bam_csi_index = false\n"
        f"includeConfig '{ROOT}/conf/modules/prepare_realignment/vcf_realignment.config'\n")
    result = run_nextflow(tmp_path, f'''
params.outdir = '{tmp_path}/output'
params.publish_dir_mode = 'copy'
include {{ VCF2BED }} from '{ROOT}/modules/local/vcf2bed/main'
workflow PREPARE_REALIGNMENT_VCF {{
    VCF2BED(Channel.of([[id:'RT', patient:'P', sample:'RT', status:2], file('{vcf}'), file('{tbi}')]))
    {gate}
    bed.toList().view {{ 'CANDIDATES=' + it.size() }}
}}
workflow TEST {{ PREPARE_REALIGNMENT_VCF() }}
workflow {{ TEST() }}
''', [])
    assert result.returncode == 0, result.stdout + result.stderr
    assert f'CANDIDATES={0 if empty else 1}' in result.stdout
    published = tmp_path / 'output/vcf_realignment/vcf2bed/RT/RT.bed'
    assert published.read_text() == ('' if empty else 'chr7\t9\t10\nchr7\t19\t20\n')
