#!/usr/bin/env python3
"""One-time, non-deleting organization of the audited SEQC2 artifacts.

Archives only enumerated superseded publications, never work/cache/input data.
Old paths become compatibility symlinks. Manifests bind all retained links to
their real targets and record source task commands and alignment identities.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


REPO = Path(__file__).resolve().parents[3]
SHARED = Path('/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar')
ARCHIVE = '.artifacts/archived_seqc2_workflow_outputs_20260914'
VERIFIED = REPO / 'examples/seqc2/verified/20260914'
COMPARISON = REPO / 'examples/seqc2/comparison/historical_manual_replay_20260914'
OLD_NAMES = ('seqc2.wes.ll.hybrid', 'seqc2.wes.ll.hybrid.modality',
             'seqc2.wes.ll.hybrid.pooling_fix', 'seqc2.wes.ll.hybrid.realign.full',
             'seqc2.wes.ll.hybrid.realign.latest', 'seqc2.wes.ll.hybrid.realign.preview')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def inventory(root):
    result = {}
    for folder, dirs, files in os.walk(root):
        for name in dirs + files:
            path = Path(folder) / name
            # Preflight found no symlinks in these publication directories.
            # Fail rather than change the meaning of a relative link on move.
            if path.is_symlink():
                raise ValueError(f'Unexpected publication symlink: {path}')
            stat = path.stat()
            result[str(path.relative_to(root))] = [stat.st_dev, stat.st_ino, stat.st_size]
    return result


def task_manifest(output, work):
    trace = sorted((output / 'pipeline_info').glob('execution_trace_*.txt'))[-1]
    result = []
    for row in csv.DictReader(trace.open(), delimiter='\t'):
        if not any(':' + name + ' (' in row['name'] for name in
                   ('DEEPSOMATIC', 'MUTECT2_PAIRED', 'STRELKA_SOMATIC')):
            continue
        prefix, suffix = row['hash'].split('/')
        matches = list((work / prefix).glob(suffix + '*'))
        if len(matches) != 1:
            raise ValueError(f'Cannot resolve task: {row}')
        task = matches[0]
        command = task / '.command.sh'
        alignments = []
        for path in task.iterdir():
            if path.name.endswith(('.bam', '.cram')):
                target = path.resolve(strict=True)
                extension = '.bai' if path.name.endswith('.bam') else '.crai'
                index = Path(str(path) + extension).resolve(strict=True)
                alignments.append({'staged': str(path), 'target': str(target),
                                   'index': str(index),
                                   'size': target.stat().st_size})
        if not alignments:
            raise ValueError(f'Missing task alignment inputs: {task}')
        result.append({'name': row['name'], 'task_hash': row['hash'],
                       'workdir': str(task), 'status': row['status'],
                       'command': command.read_text(), 'command_sha256': sha(command),
                       'alignments': alignments})
    return {'trace': str(trace), 'trace_sha256': sha(trace), 'tasks': result}


def main():
    previous = VERIFIED / 'OUTPUT_MANIFEST.json'
    if previous.exists() and json.loads(previous.read_text()).get('status') == 'organized_and_verified':
        raise ValueError(f'Already organized: {VERIFIED}')
    # An interrupted organization can resume against its exact archive paths.
    wes_dna = REPO / 'examples/seqc2/hybrid/output/seqc2.wes.ll.hybrid.realign.policy-default'
    wes_rescue = REPO / 'examples/seqc2/hybrid/output/seqc2.wes.ll.hybrid.realign.latest'
    wgs = SHARED / 'examples/seqc2/hybrid/output/seqc2.wgs.il.hybrid'
    work = REPO / 'examples/seqc2/hybrid/work'
    shared_work = SHARED.parent / 'nf_work'
    lineage = {'wes_ll_dna': task_manifest(wes_dna, work),
               'wes_ll_rescue': task_manifest(wes_rescue, work),
               'wgs_il': task_manifest(wgs, shared_work)}
    plans = []
    for repo in (REPO, SHARED):
        names = OLD_NAMES + (('seqc2.wgs.il.hybrid',) if repo == REPO else ())
        for name in names:
            source = repo / 'examples/seqc2/hybrid/output' / name
            target = repo / ARCHIVE / name
            if source.is_symlink() and source.resolve() == target and target.is_dir():
                plans.append((source, target, inventory(target)))
            elif source.is_dir() and not source.is_symlink() and not target.exists():
                plans.append((source, target, inventory(source)))
            else:
                raise ValueError(f'Unexpected archive state: {source} -> {target}')
    # No caller execution; quickcheck checks existing BAM/CRAM headers and EOF.
    alignments = sorted({a['target'] for group in lineage.values()
                         for task in group['tasks'] for a in task['alignments']})
    subprocess.run(['samtools', 'quickcheck', '-v', *alignments], check=True)
    VERIFIED.mkdir(parents=True, exist_ok=True)
    manifest = {'purpose': 'Verified manual VCFs and exact upstream workflow lineage',
                'optimization': 'paused; organization only', 'lineage': lineage,
                'archive': [], 'links': [], 'snapshots': [],
                'alignment_quickcheck': {'passed': True, 'paths': alignments}}
    manifest_path = VERIFIED / 'OUTPUT_MANIFEST.json'

    def save():
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')

    save()
    for source, target, before in plans:
        target.parent.mkdir(parents=True, exist_ok=True)
        if not source.is_symlink():
            source.rename(target)
            source.symlink_to(target, target_is_directory=True)
        if inventory(target) != before:
            raise ValueError(f'Archive inventory changed: {target}')
        manifest['archive'].append({'original': str(source), 'archive': str(target),
                                    'legacy_path': 'compatibility symlink; archived, not active',
                                    'entries': len(before), 'inode_size_inventory': before})
        save()

    def link(relative, target):
        target = Path(target).resolve(strict=True)
        path = VERIFIED / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.is_symlink():
            assert path.resolve(strict=True) == target, str(path)
        else:
            path.symlink_to(target, target_is_directory=target.is_dir())
        manifest['links'].append({'link': str(path), 'target': str(target)})

    def alignment_links(dataset, group, realign_group=None):
        sources = {}
        for task in lineage[group]['tasks']:
            if 'DEEPSOMATIC (' not in task['name']:
                continue
            for a in task['alignments']:
                p = Path(a['target'])
                if '_N_1-' in p.name: sources['dna_normal.bam'] = p
                elif '_T_1-' in p.name: sources['dna_tumor.bam'] = p
                elif '_realign.' in p.name: sources['rna_realign.cram'] = p
                elif '.recal.cram' in p.name: sources['rna_tumor.cram'] = p
        if realign_group:
            for task in lineage[realign_group]['tasks']:
                if 'DEEPSOMATIC (' in task['name']:
                    for a in task['alignments']:
                        if '_realign.' in Path(a['target']).name:
                            sources['rna_realign.cram'] = Path(a['target'])
        for name, path in sources.items():
            link(f'{dataset}/alignments/{name}', path)
            extension = '.bai' if name.endswith('.bam') else '.crai'
            index = next(a['index'] for group_data in lineage.values()
                         for task in group_data['tasks'] for a in task['alignments'] if a['target'] == str(path))
            link(f'{dataset}/alignments/{name}{extension}', Path(index))

    link('comparison', COMPARISON)
    link('wes_ll/workflow_dna_source', wes_dna)
    link('wes_ll/workflow_rescue_source', wes_rescue)
    link('wgs_il/workflow_source', wgs)
    for dataset, folder in [('wes_ll', 'wes_ll'), ('wgs_il', 'wgs_historical_scope')]:
        for label, query in [('native_consensus', 'historical_native'),
                             ('native_gated_rescue', 'historical_native_gated')]:
            for ext in ('.vcf.gz', '.vcf.gz.tbi'):
                link(f'{dataset}/{label}{ext}', COMPARISON / folder / (query + ext))
        source = wes_dna if dataset == 'wes_ll' else wgs
        rescue_source = wes_rescue if dataset == 'wes_ll' else wgs
        link(f'{dataset}/caller_vcfs', source / 'variant_calling')
        link(f'{dataset}/normalized_caller_vcfs', source / 'normalized')
        link(f'{dataset}/realignment_caller_vcfs', rescue_source / 'vcf_realignment/variant_calling')
        link(f'{dataset}/original_workflow_consensus', source / 'consensus')
        link(f'{dataset}/original_workflow_rescue', rescue_source / 'vcf_realignment/rescue')
    alignment_links('wes_ll', 'wes_ll_dna', 'wes_ll_rescue')
    alignment_links('wgs_il', 'wgs_il')
    snapshots = {'wes_native_original': Path('/tmp/consensus_experiments/wesll/policy_relaxed/query.vcf.gz'),
                 'wgs_native_original': Path('/tmp/consensus_experiments/wgsil/policy_relaxed/query.vcf.gz'),
                 'wgs_original_consensus': Path('/tmp/consensus_experiments/wgsil/audit/consensus.vcf.gz'),
                 'wgs_deepsomatic_pass': Path('/tmp/consensus_experiments/wgsil/audit/ds.vcf.gz')}
    snaproot = VERIFIED / 'provenance_inputs'
    snaproot.mkdir(exist_ok=True)
    for name, source in snapshots.items():
        dest = snaproot / (name + '.vcf.gz')
        shutil.copy2(source, dest)
        if sha(source) != sha(dest):
            raise ValueError(f'Copy mismatch: {source}')
        if Path(str(source)+'.tbi').exists():
            shutil.copy2(Path(str(source)+'.tbi'), Path(str(dest)+'.tbi'))
        manifest['snapshots'].append({'source': str(source), 'copy': str(dest), 'sha256': sha(dest)})
    bed = Path('/tmp/wgsil_sites.bed')
    shutil.copy2(bed, snaproot / 'wgs_candidate_positions.bed')
    manifest['snapshots'].append({'source': str(bed), 'copy': str(snaproot/'wgs_candidate_positions.bed'), 'sha256': sha(bed)})
    for item in manifest['links']:
        assert Path(item['link']).exists(), item
    manifest['status'] = 'organized_and_verified'
    save()
    print(f'Archived {len(plans)} publications without deletion; verified {len(alignments)} alignments.')
    print(f'Verified bundle: {VERIFIED}')


if __name__ == '__main__':
    main()
