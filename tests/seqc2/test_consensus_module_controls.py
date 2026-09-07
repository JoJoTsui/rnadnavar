"""Optional local Nextflow execution of effective consensus/rescue controls (ticket 02).

Run with RUN_LOCAL_NEXTFLOW_TESTS=1 and Nextflow, bgzip, tabix and the repository
Python dependencies on PATH. Only two one-site synthetic VCFs are processed.
"""
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.skipif(
    os.environ.get("RUN_LOCAL_NEXTFLOW_TESTS") != "1",
    reason="opt-in tiny local Nextflow module execution",
)


@pytest.mark.parametrize("snv_threshold", [None, 0])
def test_explicit_zero_controls_reach_consensus_and_rescue(tmp_path, snv_threshold):
    for tool in ("nextflow", "bgzip", "tabix"):
        if not shutil.which(tool):
            pytest.skip(f"{tool} is unavailable")
    from cyvcf2 import VCF

    template = (
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr1,length=10000>\n'
        '##FILTER=<ID=PASS,Description="Passed">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n'
        'chr1\t1000\t.\tA\tG\t60\tPASS\t.\tGT:AD:DP\t0/0:100,0:100\t0/1:99,1:100\n'
    )
    for caller in ("mutect2", "deepsomatic"):
        source = tmp_path / f"sample.{caller}.vcf"
        source.write_text(template)
        subprocess.run(["bgzip", str(source)], check=True)
        subprocess.run(["tabix", "-p", "vcf", str(source) + ".gz"], check=True)

        for suffix in (".vcf.gz", ".vcf.gz.tbi"):
            shutil.copyfile(tmp_path / f"sample.{caller}{suffix}",
                            tmp_path / f"rna.{caller}{suffix}")
    rna_consensus = tmp_path / "rna.consensus.vcf"
    rna_consensus.write_text(
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chr1,length=10000>\n'
        '##FILTER=<ID=Somatic,Description="Somatic">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        'chr1\t1000\t.\tA\tG\t60\tSomatic\t.\n'
    )
    subprocess.run(["bgzip", str(rna_consensus)], check=True)
    subprocess.run(["tabix", "-p", "vcf", str(rna_consensus) + ".gz"], check=True)

    script = tmp_path / "main.nf"
    script.write_text(
        f"include {{ VCF_CONSENSUS }} from '{ROOT}/modules/local/vcf_consensus/main'\n"
        f"include {{ VCF_RESCUE }} from '{ROOT}/modules/local/vcf_rescue/main'\n"
        "workflow {\n"
        "  panel = ['mutect2', 'deepsomatic']\n"
        '  vcfs = panel.collect { file("${projectDir}/sample.${it}.vcf.gz") }\n'
        '  tbis = panel.collect { file("${projectDir}/sample.${it}.vcf.gz.tbi") }\n'
        '  rna_vcfs = panel.collect { file("${projectDir}/rna.${it}.vcf.gz") }\n'
        '  rna_tbis = panel.collect { file("${projectDir}/rna.${it}.vcf.gz.tbi") }\n'
        '  rna_consensus = file("${projectDir}/rna.consensus.vcf.gz")\n'
        '  rna_index = file("${projectDir}/rna.consensus.vcf.gz.tbi")\n'
        "  inputs = Channel.of(\n"
        "    [[id: 'default'], vcfs, tbis, panel, panel],\n"
        "    [[id: 'zero'], vcfs, tbis, panel, panel])\n"
        "  VCF_CONSENSUS(inputs)\n"
        "  rescue_inputs = VCF_CONSENSUS.out.vcf.map { meta, vcf, tbi ->\n"
        "    [meta, vcf, tbi, rna_consensus, rna_index, vcfs, tbis, panel, rna_vcfs, rna_tbis, panel]\n"
        "  }\n"
        "  VCF_RESCUE(rescue_inputs)\n"
        "}\n"
    )
    config = tmp_path / "nextflow.config"
    config.write_text(
        "process.executor = 'local'\n"
        "process.cpus = 1\n"
        "process.maxForks = 1\n"
        "process.ext.min_alt_support = { meta.id == 'zero' ? 0 : null }\n"
        "process.ext.rescue_min_dna_callers = { meta.id == 'zero' ? 0 : null }\n"
        "process.ext.rescue_min_rna_callers = { meta.id == 'zero' ? 0 : null }\n"
        'process.publishDir = [path: "${projectDir}/outputs", mode: \'copy\']\n'
        "docker.enabled = false\n"
        "singularity.enabled = false\n"
        "conda.enabled = false\n"
    )
    if snv_threshold is not None:
        with config.open("a") as handle:
            handle.write(f"process.ext.snv_thr = {snv_threshold}\n")
    environment = dict(os.environ)
    environment["PATH"] = os.pathsep.join([
        str(Path(sys.executable).parent), str(ROOT / "bin"), environment["PATH"]
    ])
    environment["NXF_OFFLINE"] = "true"
    completed = subprocess.run(
        ["nextflow", "-C", str(config), "run", str(script), "-ansi-log", "false"],
        cwd=tmp_path, env=environment, capture_output=True, text=True, timeout=120,
    )
    if snv_threshold == 0:
        assert completed.returncode != 0
        assert "snv_thr must be > 0" in completed.stdout + completed.stderr
        assert not list((tmp_path / "outputs").glob("*.consensus.vcf.gz"))
        return
    assert completed.returncode == 0, completed.stdout + completed.stderr
    with VCF(str(tmp_path / "outputs/zero.consensus.vcf.gz")) as reader:
        record = next(reader)
        assert record.INFO.get("N_SUPPORT_CALLERS") == 2
        assert record.FILTER == "Somatic"
    with VCF(str(tmp_path / "outputs/default.consensus.vcf.gz")) as reader:
        record = next(reader)
        assert record.INFO.get("N_SUPPORT_CALLERS") == 0
        assert record.FILTER != "Somatic"
    commands = "\n".join(p.read_text() for p in (tmp_path / "work").glob("*/*/.command.sh"))
    assert "--min_alt_support 0" in commands
    assert "--min_alt_support 3" in commands
    assert "--rescue_min_dna_callers 0" in commands
    assert "--rescue_min_dna_callers 1" in commands
    assert "--rescue_min_rna_callers 0" in commands
    assert "--rescue_min_rna_callers 1" in commands
    for sample in ("default", "zero"):
        assert (tmp_path / f"outputs/{sample}.rescued.vcf.gz.tbi").is_file()
    with VCF(str(tmp_path / "outputs/zero.rescued.vcf.gz")) as reader:
        assert next(reader).FILTER == "Somatic"
