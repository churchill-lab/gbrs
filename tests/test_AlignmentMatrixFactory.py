import pytest
import os
from scipy.sparse import coo_matrix
import tables
import numpy as np
from pathlib import Path
from gbrs.emase.AlignmentMatrixFactory import AlignmentMatrixFactory
from gbrs.emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix

THIS_DIR = Path(__file__).parent


@pytest.fixture
def bam_file(tmp_path):
    file_path = THIS_DIR / 'r1.bam'
    return file_path


@pytest.fixture
def loci():
    return ['ENSMUST00000000001', 'ENSMUST00000000002']


@pytest.fixture
def haplotypes():
    return ['A', 'B']


def test_prepare(bam_file, loci, haplotypes, tmp_path):
    factory = AlignmentMatrixFactory(bam_file)
    factory.prepare(haplotypes, loci, outdir=tmp_path)

    assert factory.hname == haplotypes
    assert factory.lname == loci
    assert len(factory.rname) > 0
    assert len(factory.tmpfiles) == len(haplotypes)
    for hap in haplotypes:
        assert os.path.exists(factory.tmpfiles[hap])


def test_produce(bam_file, loci, haplotypes, tmp_path):
    factory = AlignmentMatrixFactory(bam_file)
    factory.prepare(haplotypes, loci, outdir=tmp_path)

    h5file = tmp_path / 'output.h5'
    factory.produce(h5file)

    assert os.path.exists(h5file)
    apm = AlignmentPropertyMatrix(h5file=h5file)
    # this class formats the h5 into a specific format for use downstream.
    # testing this class is done seperately, and is used here to obtain correct matrix formatting.

    test_matrix = np.array([0, 0, 1, 0, 2, 0, 3, 0, 4, 0])
    test_matrix = test_matrix.reshape((int(len(test_matrix) / 2), 2)).T

    test_sparse_mat = coo_matrix(
        (np.ones(test_matrix.shape[1]), test_matrix[:2]), shape=(6, 2)  # number of reads and number of loci
    ).tocsc()

    assert not (test_sparse_mat != apm.data[0]).A.all()
    # the assertion here is comparing two sparse matricies.
    # A warning is raised when using `==` to compare sparse matricies, so != is used per warning recommendation.

    with tables.open_file(h5file, 'r') as h5fh:
        assert h5fh.root._v_attrs.incidence_only
        assert h5fh.root._v_attrs.mtype == 'csc_matrix'
        assert h5fh.root._v_attrs.shape == (len(loci), len(haplotypes), len(factory.rname))
        for hid in range(len(haplotypes)):
            hap_group = getattr(h5fh.root, f'h{hid}')
            assert hap_group.indptr.shape[0] > 0
            assert hap_group.indices.shape[0] > 0
    # basic tests of the h5 file structure without class methods assigned by the AlignmentPropertyMatrix class.


def test_cleanup(bam_file, loci, haplotypes, tmp_path):
    factory = AlignmentMatrixFactory(bam_file)
    factory.prepare(haplotypes, loci, outdir=tmp_path)

    factory.cleanup()

    for hap in haplotypes:
        assert not os.path.exists(factory.tmpfiles[hap])
