"""Execute the production Nextflow channel operators without running aligners."""
import json
import gzip
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]


def run_nextflow(tmp_path, script_text, rows):
    script = tmp_path / "test.nf"
    script.write_text("import groovy.json.JsonSlurper\nimport groovy.json.JsonOutput\n" + script_text)
    executable = shutil.which("nextflow")
    if not executable:
        pytest.skip("Nextflow required for channel regression")
    return subprocess.run(
        [executable, "run", str(script), "-offline", "--rows", json.dumps(rows)],
        cwd=tmp_path, text=True, capture_output=True, timeout=90,
        env={**os.environ, "NXF_OFFLINE": "true", "NXF_ANSI_LOG": "false"},
    )


def run_pool(tmp_path, rows, branch="rna", missing=False):
    source = (ROOT / "subworkflows/local/bam_align/main.nf").read_text()
    counts = source.split("        reads_for_alignment.map { meta, reads ->", 1)[1]
    counts = "reads_for_alignment.map { meta, reads ->" + counts.split(
        "        reads_for_alignment = reads_for_alignment.map", 1
    )[0]
    pooling = source.split(f"        bam_mapped_{branch} =", 1)[1]
    pooling = f"bam_mapped_{branch} =" + pooling.split(
        f'        bam_mapped_{branch}.dump', 1
    )[0]
    result = run_nextflow(tmp_path,
        "workflow {\n"
        "def rows = new JsonSlurper().parseText(params.rows)\n"
        "reads_for_alignment = Channel.fromList(rows)\n"
        f"FASTQ_ALIGN_STAR = [out: [bam: Channel.fromList(rows{'.drop(1)' if missing else ''})]]\n"
        f"FASTQ_ALIGN = [out: [bam: Channel.fromList(rows{'.drop(1)' if missing else ''})]]\n"
        + counts + pooling + "\n"
        f"bam_mapped_{branch}.toList().view {{ items -> 'POOL_RESULT=' + JsonOutput.toJson("
        "items.collect { meta, bams -> [meta instanceof Map ? meta : meta.getGroupTarget(), bams] }) }\n}\n", rows)
    if missing:
        assert result.returncode != 0
        assert "Incomplete alignment pool" in result.stdout + result.stderr
        return
    assert result.returncode == 0, result.stdout + result.stderr
    return json.loads(next(line.split("=", 1)[1] for line in result.stdout.splitlines()
                           if line.startswith("POOL_RESULT=")))


def row(sample, status, lane, library=None, patient="P"):
    meta = dict(patient=patient, sample=sample, status=status, lane=lane,
                id=f"{sample}-{lane}", read_group=f"ID:{lane}", data_type="fastq",
                num_lanes=1, size=1)
    if library:
        meta.update(input_stage="raw_reads", library=library)
    return [meta, f"{patient}.{sample}.{lane}.bam"]


def test_hybrid_three_rna_libraries_emit_one_rt_pool(tmp_path):
    rows = [row("WES_LL_RT_1", 2, lib, lib) for lib in
            ("SRR9134699", "SRR9134697", "SRR9134727")]
    pools = run_pool(tmp_path, rows)
    assert len(pools) == 1, "RT disappears before preprocessing and RT/DN calling"
    assert sorted(pools[0][1]) == sorted(r[1] for r in rows)
    assert pools[0][0]["libraries"] == sorted(r[0]["library"] for r in rows)
    assert "library" not in pools[0][0] and "lane" not in pools[0][0]


@pytest.mark.parametrize("branch", ["dna", "rna"])
def test_seq2neo_legacy_triplet_preserves_all_roles(tmp_path, branch):
    rows = [row(sample, status, "LX") for sample, status in
            (("DN", 0), ("DT", 1), ("RT", 2))]
    for meta, _ in rows:
        meta.update(input_stage=[], library=[])
    pools = run_pool(tmp_path, rows, branch=branch)
    assert sorted((m["sample"], m["status"], b) for m, b in pools) == sorted(
        (m["sample"], m["status"], [b]) for m, b in rows)


@pytest.mark.parametrize("branch", ["dna", "rna"])
def test_legacy_multiple_lanes_and_patient_isolation(tmp_path, branch):
    rows = [row("RT", 2, lane, patient=patient)
            for patient in ("P", "Q") for lane in ("L1", "L2")]
    pools = run_pool(tmp_path, rows[::-1], branch=branch)
    assert len(pools) == 2
    for meta, bams in pools:
        assert len(bams) == 2
        assert all(b.startswith(meta["patient"] + ".") for b in bams)


def test_incomplete_pool_is_an_error(tmp_path):
    run_pool(tmp_path, [row("RT", 2, lib, lib) for lib in ("A", "B", "C")], missing=True)


def test_split_fastq_chunks_pool_with_library_provenance(tmp_path):
    rows = [row("RT", 2, "L1", "A") for _ in range(2)] + [row("RT", 2, "L2", "B")]
    rows[0][1], rows[1][1] = "A.chunk1.bam", "A.chunk2.bam"
    pools = run_pool(tmp_path, rows)
    assert len(pools) == 1 and len(pools[0][1]) == 3
    assert pools[0][0]["libraries"] == ["A", "B"]


def test_legacy_schema_empty_stage_is_not_a_hybrid_manifest(tmp_path):
    source = (ROOT / "subworkflows/local/utils_nfcore_rnadnavar_pipeline/main.nf").read_text()
    function = "def validateHybridManifest(rows) {" + source.split("def validateHybridManifest(rows) {", 1)[1]
    rows = [[dict(sample=s, status=i, input_stage=[], library=[])]
            for i, s in enumerate(("DN", "DT", "RT"))]
    result = run_nextflow(tmp_path, function + "\nworkflow {\n"
        "validateHybridManifest(new JsonSlurper().parseText(params.rows))\n"
        "println 'LEGACY_ACCEPTED'\n}\n", rows)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "LEGACY_ACCEPTED" in result.stdout


@pytest.mark.parametrize("staged", [False, True])
@pytest.mark.parametrize("realignment", [False, True])
def test_samplesheet_library_tags(tmp_path, staged, realignment):
    fastq = tmp_path / "reads.fastq.gz"
    with gzip.open(fastq, "wt") as handle:
        handle.write("@instrument:run:flowcell:1:tile:x:y\nA\n+\nI\n")
    rows = []
    for i, sample in enumerate(("DN", "DT", "RT")):
        meta = dict(patient="P", sample=sample, status=i, lane="LX",
                    input_stage="raw_reads" if staged else [],
                    library=f"{sample}_lib" if staged else [])
        rows.append([meta, str(fastq), str(fastq)] + [[] for _ in range(8)])
    module = ROOT / "subworkflows/local/samplesheet_to_channel/main.nf"
    alignment = (ROOT / "subworkflows/local/bam_align/main.nf").read_text()
    ingress = "input_sample_type = input_sample.branch{" + alignment.split(
        "input_sample_type = input_sample.branch{", 1)[1].split("// QC & TRIM", 1)[0]
    routing = "reads_for_alignment_status = reads_for_alignment.branch{" + alignment.split(
        "reads_for_alignment_status = reads_for_alignment.branch{", 1)[1].split("//  DNA mapping", 1)[0]
    tools = "'mutect2,strelka,deepsomatic,consensus,rescue,realignment'" if realignment else "null"
    result = run_nextflow(tmp_path,
        "params.step='mapping'\nparams.aligner='bwa-mem'\nparams.dbsnp='fixture'\n"
        f"params.tools={tools}\nparams.fasta='ref.fa'\nparams.seq_platform='ILLUMINA'\n"
        f"include {{ SAMPLESHEET_TO_CHANNEL }} from '{module}'\n"
        "workflow {\n"
        "def rows = new JsonSlurper().parseText(params.rows)\n"
        "SAMPLESHEET_TO_CHANNEL(Channel.fromList(rows).map { row -> "
        "[row[0], file(row[1]), file(row[2])] + row.drop(3) })\n"
        "SAMPLESHEET_TO_CHANNEL.out.input_sample.toList().view { items -> "
        "'TAGS=' + JsonOutput.toJson(items.collect { it[0] }) }\n"
        "input_sample = SAMPLESHEET_TO_CHANNEL.out.input_sample\n"
        + ingress + "\nreads_for_alignment = input_sample_type.fastq\n" + routing + "\n"
        "reads_for_alignment_status.dna.toList().view { 'DNA=' + JsonOutput.toJson(it.collect { it[0].sample }.sort()) }\n"
        "reads_for_alignment_status.rna.toList().view { 'RNA=' + JsonOutput.toJson(it.collect { it[0].sample }.sort()) }\n"
        "input_sample_type.caller_ready_bam.toList().view { 'BYPASS=' + it.size() }\n}\n", rows)
    assert result.returncode == 0, result.stdout + result.stderr
    metas = json.loads(next(line.split("=", 1)[1] for line in result.stdout.splitlines()
                            if line.startswith("TAGS=")))
    assert len(metas) == 3
    assert 'DNA=["DN","DT"]' in result.stdout
    assert 'RNA=["RT"]' in result.stdout
    assert 'BYPASS=0' in result.stdout
    for meta in metas:
        library = meta["sample"] + ("_lib" if staged else "")
        assert f"LB:{library}" in meta["read_group"]
        assert f'SM:{meta["sample"]}' in meta["read_group"]


def test_partially_staged_schema_manifest_still_fails(tmp_path):
    source = (ROOT / "subworkflows/local/utils_nfcore_rnadnavar_pipeline/main.nf").read_text()
    function = "def validateHybridManifest(rows) {" + source.split("def validateHybridManifest(rows) {", 1)[1]
    result = run_nextflow(tmp_path, function + "\nworkflow {\n"
        "validateHybridManifest(new JsonSlurper().parseText(params.rows))\n}\n",
        [[dict(input_stage="raw_reads")], [dict(input_stage=[])]])
    assert result.returncode != 0
    assert "partially staged" in result.stdout + result.stderr


@pytest.mark.parametrize("hybrid", [False, True])
def test_both_modalities_reach_caller_pair_channel(tmp_path, hybrid):
    rows = [row("DN", 0, "LX"), row("DT", 1, "LX")]
    if hybrid:
        rt_rows = [row("RT", 2, lib, lib) for lib in ("A", "B", "C")]
        pools = run_pool(tmp_path, rt_rows)
        rows += [[pools[0][0], pools[0][1]]]
    else:
        rows += [row("RT", 2, "LX")]
        pools = run_pool(tmp_path, rows)
        rows = pools
    source = (ROOT / "subworkflows/local/bam_variant_calling/main.nf").read_text()
    pairing = "def build_modality_pairs =" + source.split("def build_modality_pairs =", 1)[1].split(
        '            cram_variant_calling_pair.dump', 1)[0]
    result = run_nextflow(tmp_path, "workflow {\n"
        "def rows = new JsonSlurper().parseText(params.rows)\n"
        "def alignments = Channel.fromList(rows).map { meta, bam -> "
        "[meta.patient, meta + [id: meta.sample], bam, 'index'] }\n"
        "cram_variant_calling_normal_to_cross = alignments.filter { it[1].status == 0 }\n"
        "cram_variant_calling_pair_to_cross = alignments.filter { it[1].status > 0 }\n"
        + pairing + "\ncram_variant_calling_pair.toList().view { pairs -> "
        "'PAIRS=' + JsonOutput.toJson(pairs.collect { it[0] }) }\n}\n", rows)
    assert result.returncode == 0, result.stdout + result.stderr
    pairs = json.loads(next(line.split("=", 1)[1] for line in result.stdout.splitlines()
                            if line.startswith("PAIRS=")))
    assert sorted((p["id"], p["status"]) for p in pairs) == [("DT_vs_DN", 1), ("RT_vs_DN", 2)]


def test_hybrid_completion_requires_rna_and_rescue(tmp_path):
    pytest.importorskip("yaml", reason="SEQC2 runner requires PyYAML")
    spec = importlib.util.spec_from_file_location(
        "seqc2_runner", ROOT / "examples/seqc2/scripts/run_pipeline.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    config = runner.load_config(ROOT / "examples/seqc2/hybrid/config_pooling_fix.yaml", {})
    artifacts = config["completion_artifacts"]
    assert config["outdir"].endswith(".pooling_fix")
    assert config["state_file"] == "runs/pooling_fix_state.json"
    missing = []
    for pattern in artifacts:
        path = tmp_path / pattern.replace("**", "sample").replace("*", "fixture")
        if "RT_1_vs" in pattern or pattern.startswith("rescue/"):
            missing.append(path)
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
    assert not runner.is_complete(tmp_path, artifacts)
    for path in missing:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    assert runner.is_complete(tmp_path, artifacts)
