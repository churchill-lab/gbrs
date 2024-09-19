import pytest
import os
from pathlib import Path
from scipy.sparse import coo_matrix
import numpy as np
from gbrs.gbrs.emase_utils import compress
from gbrs.emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix

THIS_DIR = Path(__file__).parent


@pytest.fixture
def emase_file_1(tmp_path):
    file_path = THIS_DIR / 'r1.h5'
    return file_path


@pytest.fixture
def emase_file_2(tmp_path):
    file_path = THIS_DIR / 'r2.h5'
    return file_path


def test_compress(tmp_path, emase_file_1, emase_file_2):
    output_file = tmp_path / "compressed_emase_output.h5"
    emase_files = [emase_file_1, emase_file_2]

    compress(emase_files, str(output_file))

    assert os.path.exists(output_file)

    apm = AlignmentPropertyMatrix(h5file=output_file)
    # this class formats the h5 into a specific format for use downstream.
    # testing this class is done separately, and is used here to obtain correct matrix formatting.

    test_matrix = np.array([0, 0, 1, 0])
    test_matrix = test_matrix.reshape((int(len(test_matrix) / 2), 2)).T

    test_sparse_mat = coo_matrix(
        (np.ones(test_matrix.shape[1]), test_matrix[:2]), shape=(4, 2)  # number of reads and number of loci
    ).tocsc()
    # compress collapes equvalent reads into a single read, so the number of reads is reduced

    assert (test_sparse_mat.A == apm.data[0].A).all()
