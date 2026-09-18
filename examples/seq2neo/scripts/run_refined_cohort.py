#!/usr/bin/env python3
"""Opt-in candidate-label rerun from six caller VCFs and annotated realignment rescue.

Preparation is read-only unless --plan is supplied. --execute is mandatory for
generation; no Nextflow, mapping, variant calling, new VEP or label cleaning.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import pysam

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / 'bin'))
from vcf_utils.aggregation import resolve_tumor_sample_index, _normal_sample_from_header

CALLERS = ('deepsomatic', 'mutect2', 'strelka')
GIB = 1024**3
POLICY = 'seqc2_refined_v2+seqc2_refined_gate_v1'
THREE_CLASS_POLICY = 'separated_three_class_v2'
EOF = bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save(path, data):
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(data, indent=2) + '\n')
    temp.replace(path)


def below(path, parent):
    return path == parent or parent in path.parents


def nearest_existing(path):
    while not path.exists():
        path = path.parent
    return path


def resource_plan(cfg, total_input_bytes):
    workers = int(cfg['workers'])
    memory = int(cfg['memory_gib_per_worker'])
    if workers < 1 or memory < 1:
        raise ValueError('Positive workers and memory required')
    caps = []
    for p in (Path('/sys/fs/cgroup/memory.max'), Path('/sys/fs/cgroup/memory/memory.limit_in_bytes')):
        if p.exists() and p.read_text().strip().isdigit():
            caps.append(int(p.read_text().strip()))
    required_memory = (workers * memory + int(cfg['reserve_memory_gib'])) * GIB
    if caps and required_memory > min(caps):
        raise ValueError('Worker budgets plus reserve exceed cgroup memory limit')
    quota_cpus = None
    if Path('/sys/fs/cgroup/cpu.max').exists():
        quota, period = Path('/sys/fs/cgroup/cpu.max').read_text().split()
        if quota != 'max': quota_cpus = int(quota) / int(period)
    elif Path('/sys/fs/cgroup/cpu/cpu.cfs_quota_us').exists():
        quota = int(Path('/sys/fs/cgroup/cpu/cpu.cfs_quota_us').read_text())
        if quota > 0: quota_cpus = quota / int(Path('/sys/fs/cgroup/cpu/cpu.cfs_period_us').read_text())
    if workers > (quota_cpus or os.cpu_count() or 1):
        raise ValueError('Workers exceed available CPU quota')
    required_disk = max(int(cfg['minimum_free_gib']) * GIB, int(total_input_bytes * cfg['disk_expansion_factor']))
    free = {key: shutil.disk_usage(nearest_existing(Path(cfg[key]))).free for key in ('output_root', 'work_root')}
    return dict(workers=workers, memory_gib_per_worker=memory, cpu_threads_per_job=1,
                cgroup_memory_bytes=min(caps) if caps else None, cpu_quota=quota_cpus,
                estimated_required_free_bytes=required_disk, free_bytes=free,
                disk_ready=all(v >= required_disk for v in free.values()))


def inspect_vcf(path, reference, caller=None):
    indices = [Path(str(path) + ext) for ext in ('.tbi', '.csi') if Path(str(path) + ext).is_file()]
    if not indices:
        raise ValueError(f'Missing index: {path}')
    with path.open('rb') as handle:
        handle.seek(-28, 2)
        if handle.read() != EOF: raise ValueError(f'Missing BGZF EOF: {path}')
    with pysam.VariantFile(str(path)) as reader:
        samples = list(reader.header.samples)
        if caller:
            if not 1 <= len(samples) <= 2: raise ValueError(f'Invalid sample count: {path}')
            tumor = resolve_tumor_sample_index(samples, caller, _normal_sample_from_header(str(reader.header)))
            if tumor is None: raise ValueError(f'Unresolved tumor sample: {path}')
            needed = {'AD', 'DP'} if caller != 'strelka' else {'AU', 'CU', 'GU', 'TU'}
            if not needed <= set(reader.header.formats): raise ValueError(f'Missing native evidence fields: {path}')
            if caller == 'mutect2' and not {'TLOD', 'GERMQ'} <= set(reader.header.info):
                raise ValueError(f'Missing Mutect2 native score schema: {path}')
        else:
            if samples: raise ValueError(f'Rescue must be sampleless: {path}')
            if 'GATE_POLICY' in reader.header.info: raise ValueError(f'Refuse already gated rescue: {path}')
            if not {'GNOMAD_AF', 'REDI_CANONICAL'} <= set(reader.header.info):
                raise ValueError(f'Missing rescue annotation schema: {path}')
        for name, entry in reader.header.contigs.items():
            if name not in reference or (entry.length is not None and entry.length != reference[name]):
                raise ValueError(f'Reference dictionary mismatch {name}: {path}')
        first = next(reader, None)
        if first is not None:
            # Exercise original index read-only; no rebuilding original indexes.
            fetched = list(reader.fetch(first.contig, first.start, first.start + 1))
            if not any(str(r) == str(first) for r in fetched):
                raise ValueError(f'Index does not recover first record: {path}')
        return dict(samples=samples, contigs=len(reader.header.contigs), has_vep_csq='CSQ' in reader.header.info,
                    indices=[str(p.resolve()) for p in indices], header_and_first_record_only=True,
                    index_older_than_vcf=any(p.stat().st_mtime_ns < path.stat().st_mtime_ns for p in indices))


def required_validation_code():
    return {str(p.relative_to(REPO)) for p in (
        REPO/'bin/run_consensus_vcf.py', REPO/'bin/apply_refined_rescue.py',
        REPO/'bin/apply_three_class_labels.py', REPO/'bin/assess_negative_label_evidence.py',
        *sorted((REPO/'bin/vcf_utils').glob('*.py')), *sorted((REPO/'bin/common').glob('*.py')))}


def validation_gate(cfg):
    """Allow candidate generation only, bound to completed checks and exact code."""
    if cfg.get('candidate_only') is not True:
        raise ValueError('Three-class policy requires candidate_only=true; no training approval')
    path=Path(cfg['validation_summary'])
    path=(path if path.is_absolute() else REPO/path).resolve(strict=True)
    report=json.loads(path.read_text())
    if (report.get('policy')!=THREE_CLASS_POLICY or report.get('candidate_execution_checks_pass') is not True
            or report.get('biological_training_approved') is not False):
        raise ValueError('Require completed candidate validation, not inferred biological approval')
    datasets=report.get('datasets',{})
    if set(datasets)!={'seqc2_wes_ll','seqc2_wgs_il','hg008_wgs'}:
        raise ValueError('Incomplete validation dataset coverage')
    for dataset in datasets.values():
        stages=dataset.get('stages',{})
        if (set(stages)!={'consensus','first','realignment'}
                or any(stage.get('somatic_parity') is not True for stage in stages.values())
                or dataset.get('negative_collision_coverage',{}).get('status')
                    !='rescue_negatives_subset_of_challenged_consensus_negatives'):
            raise ValueError('Incomplete Somatic/negative validation coverage')
    code=report.get('validated_code',{})
    if not required_validation_code() <= set(code):
        raise ValueError('Incomplete validated policy code provenance')
    for rel,expected in code.items():
        source=(REPO/rel).resolve(strict=True)
        if Path(rel).is_absolute() or not below(source,REPO) or digest(source)!=expected:
            raise ValueError('Policy code differs from completed validation: '+rel)
    if not report.get('archived_sha256'):
        raise ValueError('Missing archived validation evidence')
    for rel,expected in report['archived_sha256'].items():
        source=(path.parent/rel).resolve(strict=True)
        if Path(rel).is_absolute() or not below(source,path.parent) or digest(source)!=expected:
            raise ValueError('Archived validation evidence changed: '+rel)
    return dict(summary=str(path),sha256=digest(path),candidate_only=True,training_approved=False)


def prepare(cfg):
    if cfg.get('policy') not in (POLICY,THREE_CLASS_POLICY): raise ValueError('Unsupported policy')
    gate=validation_gate(cfg) if cfg['policy']==THREE_CLASS_POLICY else None
    for key in ('manifest', 'output_root', 'work_root', 'fasta'):
        if not Path(cfg[key]).is_absolute(): raise ValueError(f'{key} must be absolute')
    manifest, fasta = Path(cfg['manifest']), Path(cfg['fasta'])
    if not fasta.is_file(): raise ValueError('Missing reference FASTA')
    reference = {v[0]: int(v[1]) for line in Path(str(fasta) + '.fai').read_text().splitlines() if (v := line.split('\t'))}
    with manifest.open() as handle: rows = list(csv.DictReader(handle, delimiter='\t'))
    ids = [r['sample_id'] for r in rows]
    if len(ids) != cfg['expected_samples'] or len(set(ids)) != len(ids):
        raise ValueError('Cohort count mismatch or duplicate IDs')
    if not set(cfg['pilot_samples']) <= set(ids): raise ValueError('Pilot sample missing')
    out, work = Path(cfg['output_root']).resolve(), Path(cfg['work_root']).resolve()
    if below(out, work) or below(work, out): raise ValueError('Output and work namespaces must be separate')
    samples = []
    for row in rows:
        sid = row['sample_id']
        if not re.fullmatch(r'[A-Za-z0-9_-]+', sid): raise ValueError('Unsafe sample ID')
        source_root = Path(row['base_output_dir']).resolve()
        for dest in (out, work):
            if below(dest, source_root) or below(source_root, dest): raise ValueError('Output/work overlaps original source root')
        paths = {f'{modality}_{caller}': str(Path(row[f'caller_{modality}_{caller}']).resolve(strict=True))
                 for modality in ('dna', 'rna') for caller in CALLERS}
        paths['rescue'] = str(Path(row['rescue_vcf_path']).resolve(strict=True))
        if any(below(Path(p), d) for p in paths.values() for d in (out, work)):
            raise ValueError('Source file is inside output/work namespace')
        if len(set(paths.values())) != 7: raise ValueError(f'{sid}: reused source file')
        if any('vcf_realignment' not in paths[key] for key in ('rescue', 'rna_deepsomatic', 'rna_mutect2', 'rna_strelka')):
            raise ValueError(f'{sid}: require matched realignment source paths')
        checks, all_files = {}, set(paths.values())
        for key, value in paths.items():
            checks[key] = inspect_vcf(Path(value), reference, key.split('_')[1] if key != 'rescue' else None)
            if gate and key=='dna_deepsomatic':
                with pysam.VariantFile(value) as reader:
                    if 'GQ' not in reader.header.formats:
                        raise ValueError(f'{sid}: missing native DeepSomatic GQ schema')
            all_files.update(checks[key]['indices'])
        stats = {p: {'bytes': Path(p).stat().st_size, 'mtime_ns': Path(p).stat().st_mtime_ns} for p in sorted(all_files)}
        samples.append(dict(sample_id=sid, inputs=paths, files=stats, checks=checks))
    return dict(status='preflight_pass_header_checks_only', policy=cfg['policy'], samples=samples,validation_gate=gate,
                manifest_sha256=digest(manifest), reference=str(fasta),
                resources=resource_plan(cfg, sum(v['bytes'] for s in samples for v in s['files'].values())),
                limitations=['Full-record parsing occurs at execution', 'No new VEP; inherited annotations only',
                             'Candidate labels are not approved training labels'])


def code_hashes():
    paths = [Path(__file__), REPO / 'bin/run_consensus_vcf.py', REPO / 'bin/apply_refined_rescue.py',
             REPO / 'bin/apply_three_class_labels.py',
             REPO / 'examples/seqc2/scripts/audit_refined_label_contract.py',
             REPO / 'examples/seqc2/scripts/validate_refined_native_integration.py',
             *sorted((REPO / 'bin/vcf_utils').glob('*.py')), *sorted((REPO / 'bin/common').glob('*.py'))]
    return {str(p.relative_to(REPO)): digest(p) for p in paths}


def completed(sample_root, identity, source_hashes):
    state = sample_root / 'state.json'
    if not state.exists(): return False
    d = json.loads(state.read_text())
    if identity.get('policy')==THREE_CLASS_POLICY:
        paths=d.get('final_artifacts',{})
        dest=Path(d.get('output',''))
        if (set(paths)!={'consensus','rescue'} or dest.resolve().parent!=sample_root.resolve()
                or any(Path(paths[k])!=dest/f'three_class.{k}.vcf.gz'
                       or paths[k] not in d.get('outputs',{}) for k in paths)):
            return False
    return (d.get('status') == 'candidate_complete_not_training_approved' and d.get('identity') == identity
            and d.get('source_hashes') == source_hashes and bool(d.get('outputs'))
            and all(Path(p).is_file() and digest(p) == h for p, h in d['outputs'].items()))


def run_sample(sample, cfg, identity):
    if cfg['policy']==THREE_CLASS_POLICY and validation_gate(cfg)['sha256']!=identity.get('validation_gate'):
        raise ValueError('Validation gate changed since preparation')
    sid = sample['sample_id']
    root = Path(cfg['output_root']) / sid
    root.mkdir(parents=True, exist_ok=True)
    hashes = {p: digest(p) for p in sample['files']}
    if completed(root, identity, hashes): return sid, 'cached_candidate_complete'
    attempt = max([int(p.name[7:]) for p in root.glob('attempt[0-9]*') if p.name[7:].isdigit()] or [0]) + 1
    dest = root / f'attempt{attempt:03d}'
    dest.mkdir()
    namespace = hashlib.sha256(str(Path(cfg['output_root']).resolve()).encode()).hexdigest()[:16]
    work = Path(cfg['work_root']) / namespace / sid / dest.name
    work.mkdir(parents=True, exist_ok=False)
    state = dict(status='running', sample_id=sid, identity=identity, source_hashes=hashes,
                 work=str(work), output=str(dest), commands=[], started=time.time())
    save(root / 'state.json', state)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    def run(cmd, name):
        limited = ['prlimit', f'--as={int(cfg["memory_gib_per_worker"])*GIB}', '--', *map(str, cmd)]
        state['commands'].append(limited)
        save(root / 'state.json', state)
        with (dest / f'{name}.log').open('w') as log:
            subprocess.run(limited, env=env, cwd=work, stdin=subprocess.DEVNULL, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
    try:
        inputs = work / 'dna_inputs'; inputs.mkdir()
        for caller in CALLERS:
            target = inputs / f'sample.{caller}.vcf.gz'
            run(['bcftools','norm','-f',cfg['fasta'],'-c','e','-d','exact',sample['inputs'][f'dna_{caller}'],
                 '-Oz','-o',target], f'normalize_{caller}')
            run(['bcftools','index','-t',target], f'index_{caller}')
        run([sys.executable,REPO/'bin/run_consensus_vcf.py','--input_dir',inputs,'--expected_callers',','.join(CALLERS),
             '--experimental-refined-native','--out_prefix',dest/'refined'], 'consensus')
        baseline = dest / 'refined.vcf.gz'
        run(['bcftools','index','-t',baseline], 'index_consensus')
        cmd = [sys.executable,REPO/'bin/apply_refined_rescue.py','--dna-consensus',baseline,
               '--annotated-rescue',sample['inputs']['rescue'],'--alignment-round','realignment','--outdir',work/'rescue']
        for modality in ('dna','rna'):
            for caller in CALLERS: cmd += [f'--{modality}-vcf', f'{caller}={sample["inputs"][f"{modality}_{caller}"]}']
        run(cmd, 'rescue')
        report = json.loads((work/'rescue/report.json').read_text())
        if not report['sources_unchanged']: raise ValueError('Adapter input integrity failure')
        for suffix in ('', '.tbi'):
            shutil.copyfile(work/'rescue'/('refined.rescue.vcf.gz'+suffix), dest/('refined.rescue.vcf.gz'+suffix))
        shutil.copyfile(work/'rescue/report.json', dest/'rescue.report.json')
        final_paths=dict(consensus=baseline,rescue=dest/'refined.rescue.vcf.gz')
        if cfg['policy']==THREE_CLASS_POLICY:
            # The established Somatic algorithms above remain unchanged. The
            # three-class experiment is used only for independent nominations.
            run([sys.executable,REPO/'bin/run_consensus_vcf.py','--input_dir',inputs,
                 '--expected_callers',','.join(CALLERS),'--experimental-refined-native',
                 '--experimental-three-class','--out_prefix',dest/'native'], 'native_candidates')
            native=dest/'native.vcf.gz'
            run(['bcftools','index','-t',native], 'index_native')
            for name,source in list(final_paths.items()):
                folder=work/('three_class_'+name)
                run([sys.executable,REPO/'bin/apply_three_class_labels.py','--somatic-baseline',source,
                     '--native-candidates',native,'--expected-baseline-sha256',digest(source),
                     '--expected-native-sha256',digest(native),'--stage',
                     'consensus' if name=='consensus' else 'realignment','--outdir',folder], 'three_class_'+name)
                outcome=json.loads((folder/'report.json').read_text())
                if (outcome['status']!='candidate_complete_not_training_approved'
                        or outcome['somatic_membership_mismatches'] or not outcome['sources_unchanged']):
                    raise ValueError('Separated candidate contract failed')
                target=dest/f'three_class.{name}.vcf.gz'
                for suffix in ('','.tbi'):
                    shutil.copyfile(str(folder/'candidates.vcf.gz')+suffix,str(target)+suffix)
                if digest(target)!=outcome['output_sha256']:
                    raise ValueError('Candidate publication checksum mismatch')
                shutil.copyfile(folder/'report.json',dest/f'three_class.{name}.report.json')
                final_paths[name]=target
        state['final_artifacts']={k:str(v) for k,v in final_paths.items()}
        for name, path in final_paths.items():
            run([sys.executable,REPO/'examples/seqc2/scripts/audit_refined_label_contract.py','--vcf',path,
                 '--report',dest/f'{name}.audit.json'], f'audit_{name}')
            audit = json.loads((dest/f'{name}.audit.json').read_text())
            if audit['issues'] or not audit['sources_unchanged']: raise ValueError('Output structural audit failed')
        if any(digest(p) != h for p,h in hashes.items()): raise ValueError('Source checksum changed')
        if code_hashes() != identity['code']: raise ValueError('Code changed during run')
        if cfg['policy']==THREE_CLASS_POLICY and validation_gate(cfg)['sha256']!=identity['validation_gate']:
            raise ValueError('Validation gate changed during execution')
        if (digest(cfg['manifest']) != identity['manifest'] or digest(cfg['fasta']) != identity['reference']
                or digest(cfg['fasta']+'.fai') != identity['reference_fai']):
            raise ValueError('Manifest/reference changed during run')
        state.update(status='candidate_complete_not_training_approved', completed=time.time(),
                     outputs={str(p):digest(p) for p in dest.iterdir() if p.name.endswith(('.vcf.gz','.tbi','.json'))})
        save(root/'state.json',state)
        save(dest/'completion.json',state)
        return sid, state['status']
    except Exception as exc:
        state.update(status='failed',error=f'{type(exc).__name__}: {exc}',finished=time.time())
        save(root/'state.json',state); save(dest/'failure.json',state)
        raise


def write_candidate_manifest(root, samples):
    """Inventory successes only; a pilot never implies cohort completion."""
    target = root / 'candidate_manifest.tsv'
    temp = target.with_suffix('.tmp')
    with temp.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample_id', 'status', 'dna_consensus', 'truth_vcf', 'completion_report',
                         'candidate_vcf','training_label_vcf','policy'])
        for sample in samples:
            path = root / sample['sample_id'] / 'state.json'
            if not path.exists(): continue
            state = json.loads(path.read_text())
            if state.get('status') != 'candidate_complete_not_training_approved': continue
            dest = Path(state['output'])
            paths=state.get('final_artifacts',dict(consensus=str(dest/'refined.vcf.gz'),rescue=str(dest/'refined.rescue.vcf.gz')))
            policy=state['identity']['policy']
            if policy==THREE_CLASS_POLICY and 'final_artifacts' not in state:
                raise ValueError('Missing separated final artifacts')
            if policy==THREE_CLASS_POLICY and (set(paths)!={'consensus','rescue'} or any(
                    Path(paths[k])!=dest/f'three_class.{k}.vcf.gz' or paths[k] not in state['outputs'] for k in paths)):
                raise ValueError('Invalid separated final artifacts')
            writer.writerow([sample['sample_id'], state['status'], paths['consensus'],
                             '' if policy==THREE_CLASS_POLICY else paths['rescue'],dest/'completion.json',
                             paths['rescue'],'',policy])
    temp.replace(target)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, default=REPO/'examples/seq2neo/config/refined_native_v2_cohort.json')
    ap.add_argument('--plan', type=Path, help='Write a new preparation report; no generation')
    group=ap.add_mutually_exclusive_group()
    group.add_argument('--pilot', action='store_true'); group.add_argument('--all', action='store_true')
    ap.add_argument('--execute',action='store_true')
    ap.add_argument('--approve-pilot',action='store_true',help='Acknowledge reviewed pilot before the full cohort')
    args=ap.parse_args(); cfg=json.loads(args.config.read_text())
    plan=prepare(cfg)
    if args.plan:
        with args.plan.open('x') as handle: json.dump(plan,handle,indent=2); handle.write('\n')
    print(json.dumps({'status':plan['status'],'policy':plan['policy'],'training_approved':False,
                      'samples':len(plan['samples']),'resources':plan['resources']},indent=2),flush=True)
    if not args.execute: return
    if not (args.pilot or args.all): raise ValueError('Choose --pilot or --all for execution')
    if not plan['resources']['disk_ready']: raise ValueError('Insufficient estimated disk headroom')
    for executable in ('bcftools','prlimit'):
        if not shutil.which(executable): raise ValueError(f'Missing {executable}')
    root=Path(cfg['output_root']); root.mkdir(parents=True,exist_ok=True)
    with (root/'.cohort.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        identity={'policy':cfg['policy'],'manifest':plan['manifest_sha256'],'reference':digest(cfg['fasta']),
                  'reference_fai':digest(cfg['fasta']+'.fai'),'code':code_hashes()}
        if plan['validation_gate']:
            identity['validation_gate']=plan['validation_gate']['sha256']
        frozen=root/'run_identity.json'
        if frozen.exists() and json.loads(frozen.read_text()) != identity:
            raise ValueError('Changed code/input manifest/reference: use a new output namespace')
        if not frozen.exists(): save(frozen,identity)
        if args.all:
            if not args.approve_pilot: raise ValueError('Review pilot then supply --approve-pilot')
            for sample in plan['samples']:
                if sample['sample_id'] in cfg['pilot_samples']:
                    hashes={p:digest(p) for p in sample['files']}
                    if not completed(root/sample['sample_id'],identity,hashes):
                        raise ValueError('All pilot outputs must complete and pass integrity before --all')
        selected=[s for s in plan['samples'] if not args.pilot or s['sample_id'] in cfg['pilot_samples']]
        failures=[]
        with ThreadPoolExecutor(max_workers=cfg['workers']) as pool:
            futures={pool.submit(run_sample,s,cfg,identity):s['sample_id'] for s in selected}
            for future in as_completed(futures):
                try: print(*future.result(),flush=True)
                except Exception as exc: failures.append(futures[future]); print(f'FAILED {futures[future]}: {exc}',flush=True)
                write_candidate_manifest(root, plan['samples'])
        if failures: raise SystemExit(f'Failed samples retained for retry: {failures}')


if __name__=='__main__': main()
