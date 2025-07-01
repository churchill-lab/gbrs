import pytest
import os
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from gbrs.emase.commands import bam2emase_paired
from gbrs.emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix

THIS_DIR = Path(__file__).parent


@pytest.fixture
def bam_file_1(tmp_path):
    file_path = THIS_DIR / 'r1.bam'
    return file_path


@pytest.fixture
def bam_file_2(tmp_path):
    file_path = THIS_DIR / 'r2.bam'
    return file_path


@pytest.fixture
def locusid_file(tmp_path):
    file_path = THIS_DIR / 'IDs.txt'
    return file_path


def test_bam2emase_paired(tmp_path, bam_file_1, bam_file_2, locusid_file):
    output_file = tmp_path / "paired_emase_output.h5"
    bam_files = [bam_file_1, bam_file_2]
    haplotypes = ["A", "B"]

    bam2emase_paired(bam_files, haplotypes, locusid_file, output_file)

    assert os.path.exists(output_file)

    apm = AlignmentPropertyMatrix(h5file=output_file)
    # this class formats the h5 into a specific format for use downstream.
    # testing this class is done seperately, and is used here to obtain correct matrix formatting.

    test_matrix = np.array([0, 0, 1, 0, 2, 0, 3, 0, 4, 0])
    test_matrix = test_matrix.reshape((int(len(test_matrix) / 2), 2)).T

    test_sparse_mat = coo_matrix(
        (np.ones(test_matrix.shape[1]), test_matrix[:2]), shape=(6, 2)  # number of reads and number of loci
    ).tocsc()

    assert not (test_sparse_mat != apm.data[0]).A.all()
