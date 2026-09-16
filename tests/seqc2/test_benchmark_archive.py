"""The manuscript archive is light, self-consistent and reproducible."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT/'docs/manuscript/frozen_native_gate_20260916'


def test_frozen_archive_integrity():
    subprocess.run([sys.executable, str(ARCHIVE/'verify.py')], check=True, capture_output=True)
    subprocess.run(['bash','-n',str(ARCHIVE/'reproduce.sh')], check=True)
    for path in ARCHIVE.rglob('*'):
        if path.is_file():
            assert path.stat().st_size < 1024 * 1024
            assert not path.name.endswith(('.vcf.gz','.bam','.cram','.sqlite'))


def test_export_refuses_existing_snapshot():
    result = subprocess.run([sys.executable, str(ROOT/'examples/seqc2/scripts/archive_frozen_benchmarks.py'),
                             '--archive', str(ARCHIVE)], capture_output=True, text=True)
    assert result.returncode != 0
    assert 'FileExistsError' in result.stderr
