"""Synthetic alignments exercise measurements independently of release labels."""
import importlib.util
from pathlib import Path
from collections import Counter

import polars as pl
import pysam
import pytest

spec = importlib.util.spec_from_file_location("evidence", Path(__file__).parents[1] / "bin/export_variant_evidence.py")
e = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e)


def read(seq="AAAAA", cigar="5M", start=0, flag=0, mq=60, bq=30):
    r = pysam.AlignedSegment()
    r.query_name = "synthetic"
    r.query_sequence = seq
    r.query_qualities = [bq] * len(seq)
    r.reference_id = 0
    r.reference_start = start
    r.cigarstring = cigar
    r.flag = flag
    r.mapping_quality = mq
    return r


def obs(r, pos, ref, alt, star=False):
    return e.observation(r, pos, ref, alt, e.SETTINGS, star)


@pytest.mark.parametrize("seq,cigar,pos,ref,alt,state", [
    ("ACAAA", "5M", 1, "A", "C", "alt"),
    ("ACAAA", "5M", 1, "A", "G", "other"),
    ("AAAAA", "5M", 1, "A", "C", "ref"),
    ("AACAAA", "2M1I3M", 1, "A", "AC", "alt"),
    ("AAGAAA", "2M1I3M", 1, "A", "AC", "other"),
    ("AAAAA", "5M", 1, "A", "AC", "ref"),
    ("AAAA", "2M2D2M", 1, "AAA", "A", "alt"),
    ("AAAAAA", "6M", 1, "AAA", "A", "ref"),
    ("AAAAA", "2M1D3M", 1, "AAA", "A", "other"),
    ("GGACAAA", "2S5M", 1, "A", "C", "alt"),
    ("AACAA", "2M100N3M", 102, "A", "C", "alt"),
    ("AAGCAA", "2M1I3M", 2, "A", "C", "alt"),
])
def test_cigar_and_full_alleles(seq, cigar, pos, ref, alt, state):
    assert obs(read(seq, cigar), pos, ref, alt)[0] == state


@pytest.mark.parametrize("r,pos,ref,alt", [
    (read("AA", "2M"), 1, "A", "AC"),  # anchor alone is not REF or ALT
    (read("AA", "2M"), 1, "AAA", "A"),
    (read("AAAA", "2M2N2M"), 1, "AAA", "A"),  # splice is not deletion
    (read("AAAA", "2M2D2M"), 2, "A", "C"),  # deletion has no base
    (read("ANAAA"), 1, "A", "C"),
    (read(bq=19), 1, "A", "C"),
    (read(mq=19), 1, "A", "C"),
    (read(mq=255), 1, "A", "C"),
])
def test_unusable_reads(r, pos, ref, alt):
    assert obs(r, pos, ref, alt) is None


@pytest.mark.parametrize("flag", [4, 256, 512, 1024, 2048])
def test_excluded_flags(flag):
    assert obs(read(flag=flag), 1, "A", "C") is None


def test_star_mapping_quality():
    r = read(mq=255)
    assert obs(r, 1, "A", "C", True) is None
    r.set_tag("NH", 1)
    assert obs(r, 1, "A", "C", True)[2] == 60
    r.set_tag("NH", 2)
    assert obs(r, 1, "A", "C", True) is None


@pytest.mark.parametrize("flag,orientation", [(65, "f1r2"), (81, "f2r1"), (129, "f2r1"), (145, "f1r2"), (0, "unknown")])
def test_orientation(flag, orientation):
    assert obs(read(flag=flag), 1, "A", "C")[4] == orientation


def test_overlaps_count_both_and_denominator():
    acc = Counter()
    for r in [read(flag=65), read(flag=129), read("ACAAA"), read("AGAAA")]:
        e.add(acc, obs(r, 1, "A", "C"), e.SETTINGS)
    x = e.finish(acc, e.SETTINGS)
    assert (x["depth"], x["ref_count"], x["alt_count"], x["other_count"]) == (4, 2, 1, 1)
    assert x["vaf"] == 0.25 and x["vaf_denominator"] == 4
    assert x["ref_f1r2"] == x["ref_f2r1"] == 1


def test_missing_zero_and_cap_distinct():
    zero = e.finish(Counter(), e.SETTINGS)
    assert zero["depth"] == zero["alt_count"] == 0 and zero["vaf"] is None
    missing = e.unavailable("bam_absent")
    assert missing["depth"] is None and missing["alt_count"] is None
    acc = Counter()
    e.add(acc, obs(read(), 1, "A", "C"), e.SETTINGS)
    covered = e.finish(acc, e.SETTINGS)
    assert covered["status"] == "ok" and covered["vaf"] == 0.0
    assert e.finish(acc, {**e.SETTINGS, "max_depth": 0})["depth"] is None


def frame():
    return pl.DataFrame({"sample_id": ["PRJNA298376_4007"] * 2, "CHROM": ["chr1"] * 2,
                         "POS": [2, 2], "REF": ["A", "A"], "ALT": ["C", "G"],
                         "FILTER": ["Somatic", "Reference"]})


def test_exact_join_multiallelic_and_labels():
    f = frame()
    out = f.with_columns(pl.lit("train_pool").alias("pool"), pl.lit("test").alias("split"))
    e.reconcile(f, out, "train_pool")
    for bad in [out.head(1), pl.concat([out, out.head(1)]),
                out.with_columns(pl.lit("Germline").alias("FILTER")),
                out.with_columns((pl.col("POS") + 1).alias("POS")),
                out.with_columns(pl.lit("4007").alias("sample_id")),
                out.with_columns(pl.lit("reserved").alias("pool"))]:
        with pytest.raises(ValueError):
            e.reconcile(f, bad, "train_pool")


def test_reserved_precedes_chromosome_split():
    for chrom in ["chr1", "chr21", "chr22", "chrX", "chrY", "chrM", "new_contig"]:
        assert e.split("reserved", chrom) == "reserved"
    assert e.split("train_pool", "chr1") == "test"
    assert e.split("train_pool", "chr22") == "val"
    assert e.split("train_pool", "new_contig") == "train"


def test_content_bound_resume(tmp_path):
    p, receipt = tmp_path / "part.parquet", tmp_path / "part.json"
    assert not e.cache_valid(p, receipt, "x")
    p.write_bytes(b"original")
    e.atomic_json(receipt, {"identity": "x", "sha256": e.sha256(p)})
    assert e.cache_valid(p, receipt, "x")
    with pytest.raises(ValueError, match="identity"):
        e.cache_valid(p, receipt, "different_configuration_or_input")
    p.write_bytes(b"changed")
    with pytest.raises(ValueError, match="checksum"):
        e.cache_valid(p, receipt, "x")
    assert e.digest({"a": 1, "b": 2}) == e.digest({"b": 2, "a": 1})


def test_indexed_bam_integration(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 200 + "\n")
    pysam.faidx(str(fasta))
    path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(path), "wb", header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 200}]}) as b:
        for r in [read(), read("ACAAA"), read("AGAAA"), read("AACAAA", "2M1I3M"), read("AAAA", "2M2D2M")]:
            b.write(r)
    pysam.index(str(path))
    rows = frame().to_dicts()
    rows += [{**rows[0], "ALT": "AC"}, {**rows[0], "REF": "AAA", "ALT": "A"},
             {**rows[0], "POS": 100}, {**rows[0], "ALT": "C,G"},
             {**rows[0], "REF": "G"}, {**rows[0], "POS": 0}]
    with pysam.AlignmentFile(str(path), "rb") as b, pysam.FastaFile(str(fasta)) as f:
        out = e.measure(rows, b, f, e.SETTINGS)
        again = e.measure(rows, b, f, e.SETTINGS)
    assert out == again
    assert out[0]["alt_count"] == out[1]["alt_count"] == 1
    assert out[2]["alt_count"] == out[3]["alt_count"] == 1
    assert out[4]["status"] == "zero_usable_depth"
    assert out[5]["status"] == "unsupported_allele" and out[5]["alt_count"] is None
    assert out[6]["status"] == "reference_mismatch"
    assert out[7]["status"] == "invalid_coordinate"


def test_model_allowlist_has_no_labels():
    names = [f"bam_{t}_{k}" for t in e.MODALITIES for k in e.COUNTS + e.FLOATS]
    assert not set(names) & set(e.LABEL + ["pool", "split"])
    assert all(n in e.schema().names for n in names)


def test_absent_unreadable_and_missing_index(tmp_path):
    bam, info, error = e.open_bam(str(tmp_path / "absent.bam"), "sample", "DN")
    assert bam is None and error[0] == "bam_absent"
    path = tmp_path / "bad.bam"
    path.write_bytes(b"not a BAM")
    assert e.open_bam(str(path), "sample", "DN")[2][0] == "index_absent"
    Path(str(path) + ".bai").write_bytes(b"not an index")
    assert e.open_bam(str(path), "sample", "DN")[2][0] == "bam_unreadable"


def test_unavailable_keeps_all_keys():
    rows = frame().to_dicts()
    values = e.measure(rows, None, None, e.SETTINGS, input_error=("bam_absent", "registered input missing"))
    assert len(values) == len(rows)
    assert all(v["depth"] is None and v["status"] == "bam_absent" for v in values)


def test_multi_base_insertion_sequence_and_quality():
    r = read("AACTAAA", "2M2I3M")
    assert obs(r, 1, "A", "ACT")[0] == "alt"
    assert obs(r, 1, "A", "AGT")[0] == "other"
    r.query_qualities = [30, 30, 19, 30, 30, 30, 30]
    assert obs(r, 1, "A", "ACT") is None


def test_resume_mismatch_preserves_existing_provenance(tmp_path, monkeypatch):
    import json
    sid = 'PRJNA298376_4007'
    directory = tmp_path / 'development' / sid
    directory.mkdir(parents=True)
    path = directory / 'inputs.json'
    path.write_text(json.dumps({'identity': 'old_identity', 'bams': 'immutable_old_provenance'}))
    original = path.read_bytes()
    monkeypatch.setattr(e, 'sample_labels', lambda *args: frame())
    monkeypatch.setattr(e, 'open_bam', lambda *args: (None, {}, ('bam_absent', 'missing')))
    sample = {'sample_id': sid, 'pool': 'train_pool', **{field: 'missing.bam' for field in e.MODALITIES.values()}}
    with pytest.raises(ValueError, match='provenance preserved'):
        e.process_sample((str(tmp_path), sample, {'parquet': 'unused'},
                          {'pilot_rows': 0, 'deadline': 0, 'resume': True}))
    assert path.read_bytes() == original


def test_expired_budget_does_not_hash_queued_samples(tmp_path, monkeypatch):
    monkeypatch.setattr(e, 'open_bam', lambda *args: pytest.fail('Must not open or hash expired queued inputs'))
    result = e.process_sample((str(tmp_path), {'sample_id': 'sid', 'pool': 'train_pool'}, {}, {'deadline': 1}))
    assert result['rows'] == 0 and not result['complete']
    assert list(tmp_path.iterdir()) == []
