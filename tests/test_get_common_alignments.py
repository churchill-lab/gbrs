import pytest
import os
import h5py
from pathlib import Path
from gbrs.emase.emase_utils import get_common_alignments


THIS_DIR = Path(__file__).parent


@pytest.fixture
def emase_file_1(tmp_path):
    file_path = THIS_DIR / 'r1.h5'
    return file_path


@pytest.fixture
def emase_file_2(tmp_path):
    file_path = THIS_DIR / 'r2.h5'
    return file_path


def test_get_common_alignments(tmp_path, emase_file_1, emase_file_2):
    output_file = tmp_path / "common_alignments.h5"
    emase_files = [emase_file_1, emase_file_2]

    get_common_alignments(emase_files, output_file)

    assert os.path.exists(output_file)

    with h5py.File(output_file, 'r') as f:
        # https://stackoverflow.com/a/52299730
        assert (f['h1']['indices'][()] == [0, 1, 2, 3, 4, 1, 2, 3, 4]).all()
        assert (f['h1']['indptr'][()] == [0, 5, 9]).all()
        assert len(f['lname'][()]) == 2
        assert len(f['rname'][()]) == 6
