import importlib.util
from pathlib import Path
import sys

import pysam

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "examples/seqc2/scripts"))
spec = importlib.util.spec_from_file_location("three_class_validation", ROOT / "examples/seqc2/scripts/validate_three_class_policy.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_queries_preserve_class_selection_and_original_labels(tmp_path):
    source = tmp_path / "source.vcf"
    header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
              '##INFO=<ID=CLASSIFICATION_RATIONALE,Number=1,Type=String,Description="rule">\n')
    labels = ["Somatic", "Germline", "Reference", "NoConsensus"]
    header += ''.join(f'##FILTER=<ID={x},Description="class">\n' for x in labels)
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    text = header + ''.join(f'chr1\t{i+1}\t.\tA\tT\t.\t{x}\tCLASSIFICATION_RATIONALE=decision:test|class:{x}\n'
                            for i,x in enumerate(labels))
    source.write_text(text)
    counts, reasons = module.class_queries(source, tmp_path)
    assert counts == dict.fromkeys(labels, 1)
    assert reasons == {"decision:test": 4}
    assert source.read_text() == text
    for i,label in enumerate(labels[:3]):
        with pysam.VariantFile(str(tmp_path / f"{label}.query.vcf.gz")) as f:
            records = list(f)
            assert len(records) == 1
            assert records[0].pos == i + 1
            assert set(records[0].filter) == {"PASS"}


def test_empty_classes_still_produce_readable_queries(tmp_path):
    source = tmp_path / "empty.vcf"
    source.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    assert module.class_queries(source, tmp_path) == ({}, {})
    for label in module.LABELS:
        with pysam.VariantFile(str(tmp_path / f"{label}.query.vcf.gz")) as f:
            assert list(f) == []


def test_bam_pilot_abstains_at_low_depth_and_detects_conflicts():
    from pilot_three_class_bam_evidence import assess
    assert assess("Reference", {"depth":20,"alt":0}, {"depth":20,"alt":0}) == "inconclusive_paired_evidence"
    assert assess("Reference", {"depth":60,"alt":0}, {"depth":60,"alt":0}) == "paired_zero_alt_corroborated_not_truth"
    assert assess("Reference", {"depth":60,"alt":0}, {"depth":20,"alt":3}) == "alt_read_conflict"
    assert assess("Germline", {"depth":40,"alt":20}, {"depth":40,"alt":0}) == "normal_read_corroborated_not_truth"
    assert assess("Germline", {"depth":60,"alt":0}, {"depth":40,"alt":20}) == "normal_reference_conflict"


def test_normal_gvcf_missing_evidence_is_not_reference(tmp_path):
    from check_three_class_normal_gvcf import outcome
    assert outcome(("chr1", 10, "A", "T"), []) == ("inconclusive", [])
    p = tmp_path / 'normal.vcf'
    p.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n'
                 '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="gq">\n'
                 '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp">\n'
                 '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\n'
                 'chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT:GQ:DP\t0/1:40:30\n'
                 'chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT:GQ:DP\t0/0:.:30\n')
    with pysam.VariantFile(str(p)) as f:
        rows = list(f)
    assert outcome(("chr1",10,"A","T"),rows[:1])[0] == 'normal_alt_corroboration'
    assert outcome(("chr1",20,"A","T"),rows[1:])[0] == 'inconclusive'


def test_rescue_plan_routes_both_rounds_and_never_calls_workflow(tmp_path):
    from validate_three_class_rescue import plan
    callers=('deepsomatic','mutect2','strelka')
    v=dict(status='complete_somatic_truth_screen_not_training_approved',
           sources_unchanged=True,code_unchanged=True,
           commands=[['python','consensus','--out_prefix',str(tmp_path/'dna')]],
           manifest=dict(truth='truth.vcf.gz',rescues=dict(first='first.vcf.gz',realignment='realign.vcf.gz'),
                         dna={c:'dna/'+c for c in callers},rna_first={c:'first/'+c for c in callers},
                         rna_realignment={c:'realign/'+c for c in callers}))
    jobs,truth,dna=plan(v,tmp_path/'out',Path('recommended.vcf.gz'))
    assert truth=='recommended.vcf.gz'
    assert [j['alignment_round'] for j in jobs]==['first','realignment']
    for job in jobs:
        cmd=job['command']
        assert '--experimental-three-class' in cmd
        assert cmd.count('--dna-vcf')==3 and cmd.count('--rna-vcf')==3
        assert 'nextflow' not in ' '.join(cmd)
        assert cmd[cmd.index('--annotated-rescue')+1]==v['manifest']['rescues'][job['alignment_round']]
