import pytest
from gbrs.emase.EMfactory import EMfactory
from gbrs.emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix


@pytest.fixture
def mock_alignment_property_matrix():
    return AlignmentPropertyMatrix(h5file='tests/r1.compressed.h5', grpfile='tests/TranscriptGroupFile.txt')


@pytest.fixture
def em_factory(mock_alignment_property_matrix):
    return EMfactory(alignments=mock_alignment_property_matrix)


def test_prepare(em_factory):
    em_factory.prepare(pseudocount=0.1, lenfile='tests/TranscriptLengthFile.txt', read_length=100)
    assert em_factory.grp_conv_mat is not None
    assert em_factory.t2t_mat is not None
    assert em_factory.allelic_expression is not None


def test_reset(em_factory):
    em_factory.reset(pseudocount=0.1)
    assert em_factory.allelic_expression is not None


def test_get_allelic_expression(em_factory):
    em_factory.prepare(pseudocount=0.1, lenfile='tests/TranscriptLengthFile.txt', read_length=100)
    allelic_expression = em_factory.get_allelic_expression(at_group_level=False)
    assert allelic_expression is not None


def test_update_probability_at_read_level(em_factory):
    em_factory.prepare(pseudocount=0.1, lenfile='tests/TranscriptLengthFile.txt', read_length=100)
    em_factory.update_probability_at_read_level(model=4)
    assert em_factory.probability is not None


def test_update_allelic_expression(em_factory):
    em_factory.prepare(pseudocount=0.1, lenfile='tests/TranscriptLengthFile.txt', read_length=100)
    em_factory.update_allelic_expression(model=4)
    assert em_factory.allelic_expression is not None


def test_run(em_factory):
    em_factory.prepare(pseudocount=0.1, lenfile='tests/TranscriptLengthFile.txt', read_length=100)
    em_factory.run(model=4, tol=0.001, max_iters=10, verbose=False)
    assert em_factory.allelic_expression is not None


def test_report_read_counts(em_factory, tmp_path):
    em_factory.prepare(pseudocount=0.1, lenfile='tests/TranscriptLengthFile.txt', read_length=100)
    output_file = tmp_path / "read_counts.txt"
    em_factory.report_read_counts(filename=output_file)
    assert output_file.exists()


def test_report_depths(em_factory, tmp_path):
    em_factory.prepare(pseudocount=0.1, lenfile='tests/TranscriptLengthFile.txt', read_length=100)
    output_file = tmp_path / "depths.txt"
    em_factory.report_depths(filename=output_file)
    assert output_file.exists()


def test_export_posterior_probability(em_factory, tmp_path):
    em_factory.prepare(pseudocount=0.1, lenfile='tests/TranscriptLengthFile.txt', read_length=100)
    output_file = tmp_path / "posterior_probability.h5"
    em_factory.export_posterior_probability(filename=output_file)
    assert output_file.exists()

# to do: add tests for models 1,2,3. Currently only model 4 is tested.
