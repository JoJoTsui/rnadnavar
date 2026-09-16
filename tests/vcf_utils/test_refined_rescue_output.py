import json
from pathlib import Path
import subprocess
import sys
from urllib.parse import unquote

import pysam
import pytest
import sqlite3

from apply_refined_rescue import transition, stage, panel


def test_panels_reject_unknown_or_repeated_identity(tmp_path):
    p = tmp_path / 'input.vcf'
    p.write_text('')
    with pytest.raises(ValueError):
        panel([f'mutect2={p}', f'mutect2={p}'])
    with pytest.raises(ValueError):
        panel([f'unknown={p}'])
    with pytest.raises(ValueError):
        panel([f'{c}={p}' for c in ('mutect2', 'strelka', 'deepsomatic')])


def test_previous_gated_output_cannot_be_gated_recursively(tmp_path):
    p = tmp_path / 'prior.vcf'
    p.write_text('##fileformat=VCFv4.2\n'
                 '##INFO=<ID=GATE_POLICY,Number=1,Type=String,Description="Policy">\n'
                 '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    with sqlite3.connect(':memory:') as db, pytest.raises(ValueError, match='recursively'):
        stage(db, p, 'dna')


def test_verification_and_protected_negative_precedence():
    assert transition("Artifact", "Somatic", True, "gate", None)[0] == "Artifact"
    assert transition("NoConsensus", "Somatic", True, "gate", "rejected")[0] == "NoConsensus"
    assert transition(None, "Germline", False, "gate", None)[0] == "Germline"
    assert transition("Somatic", "Somatic", False, "gate", "rejected")[0] == "Somatic"
    assert transition(None, "Somatic", True, "gate", "inconclusive")[0] == "NoConsensus"
    assert transition("NoConsensus", "RNAedit", True, "gate", None)[0] == "RNAedit"


def test_standalone_union_retains_evidence_and_negatives(tmp_path):
    header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
              '##INFO=<ID=GNOMAD_AF,Number=1,Type=Float,Description="AF">\n')
    for label in ("Somatic", "NoConsensus", "Artifact", "Germline"):
        header += f'##FILTER=<ID={label},Description="{label}">\n'
    columns = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    dna = tmp_path / "dna.vcf"
    dna.write_text(header + columns + ''.join(
        f'chr1\t{pos}\t.\tA\t{alt}\t50\t{label}\t.\n'
        for pos, alt, label in [(10, "G", "Somatic"), (20, "G", "NoConsensus"),
                                (50, "G", "Germline"), (60, "AT", "Somatic")]))
    rescue = tmp_path / "rescue.vcf"
    rescue.write_text(header + columns + ''.join(
        f'chr1\t{pos}\t.\tA\tG\t50\t{label}\t{info}\n'
        for pos, label, info in [(10, "Somatic", "."), (20, "Somatic", "."),
                                  (30, "Artifact", "."), (40, "Somatic", "GNOMAD_AF=0.1"),
                                  (50, "Somatic", ".")]))
    cmd = [sys.executable, str(Path(__file__).resolve().parents[2] / 'bin/apply_refined_rescue.py'),
           '--dna-consensus', str(dna), '--annotated-rescue', str(rescue),
           '--alignment-round', 'realignment', '--outdir', str(tmp_path / 'result')]
    for modality in ('dna', 'rna'):
        for caller in ('deepsomatic', 'mutect2', 'strelka'):
            path = tmp_path / f'{modality}.{caller}.vcf'
            path.write_text(header + '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
                            + columns.rstrip() + '\tFORMAT\tTUMOR\n'
                            + ''.join(f'chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tAD\t10,3\n' for pos in (20,40,50)))
            cmd += [f'--{modality}-vcf', f'{caller}={path}']
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    with pysam.VariantFile(tmp_path / 'result/refined.rescue.vcf.gz') as reader:
        rows = {r.pos: r for r in reader}
    assert {pos: list(r.filter) for pos,r in rows.items()} == {
        10:['Somatic'],20:['Somatic'],30:['Artifact'],40:['NoConsensus'],50:['Germline'],60:['Somatic']}
    for record in rows.values():
        assert record.info['UNIFIED_FILTER'] == list(record.filter)[0]
        original = unquote(record.info['GATE_SOURCE_RECORD']).split('\t')
        assert int(original[1]) == record.pos
    assert rows[20].info['RESCUE_PROMOTED'] == 'YES'
    assert rows[30].info['RESCUE_PROMOTED'] == 'NO'
    assert unquote(rows[20].info['GATE_DNA_RECORD']).split('\t')[6] == 'NoConsensus'
    report = json.loads((tmp_path / 'result/report.json').read_text())
    assert report['sources_unchanged']
    assert report['counts']['records_written'] == 6
    assert report['status'] == 'experimental_labels_not_training_ready'
    # Refuse to replace even this adapter's own previous output.
    assert subprocess.run(cmd, capture_output=True).returncode != 0
