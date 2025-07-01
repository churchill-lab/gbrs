import pytest
import numpy as np
from pathlib import Path
from gbrs.emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix

THIS_DIR = Path(__file__).parent


@pytest.fixture
def empty_matrix():
    shape = (10, 5, 3)
    return AlignmentPropertyMatrix(shape=shape)


@pytest.fixture
def emase_file(tmp_path):
    file_path = THIS_DIR / 'r1.h5'
    return file_path


def test_initialization(empty_matrix):
    assert empty_matrix.num_loci == 10
    assert empty_matrix.num_haplotypes == 5
    assert empty_matrix.num_reads == 3


def test_copy(emase_file):
    matrix = AlignmentPropertyMatrix(h5file=emase_file)
    copied_matrix = matrix.copy()
    assert copied_matrix.shape == matrix.shape
    assert copied_matrix.num_loci == matrix.num_loci
    assert copied_matrix.num_haplotypes == matrix.num_haplotypes
    assert copied_matrix.num_reads == matrix.num_reads


def test_sum(emase_file):
    matrix = AlignmentPropertyMatrix(h5file=emase_file)
    # matrix.shape is num_loci, num_haplotypes, num_reads

    sum_locus = matrix.sum(axis=matrix.Axis.LOCUS)
    assert sum_locus.shape == (6, 2)

    sum_haplotype = matrix.sum(axis=matrix.Axis.HAPLOTYPE)
    assert sum_haplotype.shape == (6, 2)

    sum_read = matrix.sum(axis=matrix.Axis.READ)
    assert sum_read.shape == (2, 2)


def test_normalize_reads(emase_file):
    matrix = AlignmentPropertyMatrix(h5file=emase_file)
    matrix.normalize_reads(axis=matrix.Axis.LOCUS)
    assert matrix.data[0].shape == (6, 2)


def test_save_and_load(tmp_path, emase_file):
    matrix = AlignmentPropertyMatrix(h5file=emase_file)
    file_path = tmp_path / "r1.h5"
    matrix.save(file_path)

    loaded_matrix = AlignmentPropertyMatrix(h5file=file_path)
    assert loaded_matrix.shape == matrix.shape
    assert np.array_equal(loaded_matrix.data[0].toarray(), matrix.data[0].toarray())


def test_get_unique_reads(emase_file):
    matrix = AlignmentPropertyMatrix(h5file=emase_file)
    unique_reads = matrix.get_unique_reads()
    assert unique_reads.shape == matrix.shape


def test_count_unique_reads(emase_file):
    matrix = AlignmentPropertyMatrix(h5file=emase_file)
    unique_counts = matrix.count_unique_reads()
    assert unique_counts.shape == (2, 2)


def test_report_alignment_counts(tmp_path, emase_file):
    matrix = AlignmentPropertyMatrix(h5file=emase_file)
    file_path = tmp_path / "alignment_counts.txt"
    matrix.report_alignment_counts(file_path)
    assert file_path.exists()


def test_combine(emase_file):
    matrix = AlignmentPropertyMatrix(h5file=emase_file)
    combined_matrix = matrix.combine(matrix)
    assert combined_matrix.shape == (2, 2, 12)
