import pytest
import os
from pathlib import Path
from gbrs.emase.emase_utils import bam2emase
from gbrs.emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix

THIS_DIR = Path(__file__).parent


@pytest.fixture
def bam_file():
    return str(THIS_DIR / 'r1.bam')


@pytest.fixture
def locusid_file():
    return str(THIS_DIR / 'IDs.txt')


@pytest.fixture
def haplotypes():
    return ['A', 'B']


def test_bam2emase(bam_file, haplotypes, locusid_file, tmp_path):
    output_file = tmp_path / 'alignments.transcriptome.h5'
    bam2emase(
        alignment_file=bam_file,
        haplotypes=haplotypes,
        locusid_file=locusid_file,
        output_file=str(output_file)
    )

    assert os.path.exists(output_file)
    apm = AlignmentPropertyMatrix(h5file=str(output_file))
    assert apm.num_loci > 0
    assert apm.num_haplotypes == len(haplotypes)
    assert apm.num_reads > 0
