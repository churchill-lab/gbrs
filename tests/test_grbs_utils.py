# import pytest
import numpy as np
import os
from collections import OrderedDict

from gbrs.gbrs.gbrs_utils import (
    get_chromosome_info,
    get_founder_info,
    unit_vector,
    print_vecs,
    get_genotype_probability,
    ris_step,
    get_transition_prob
)


def test_get_chromosome_info(monkeypatch):
    def mock_loadtxt(*args, **kwargs):
        return np.array([(b'chr1', 100000), (b'chr2', 200000)], dtype=[('f0', 'S8'), ('f1', '<i4')])

    monkeypatch.setattr(np, 'loadtxt', mock_loadtxt)
    result = get_chromosome_info()
    expected = OrderedDict({'chr1': 100000, 'chr2': 200000})
    assert result == expected


def test_get_founder_info(monkeypatch):
    def mock_loadtxt(*args, **kwargs):
        return np.array([('A', '#FF0000'), ('B', '#00FF00')], dtype=[('f0', 'U1'), ('f1', 'U7')])

    monkeypatch.setattr(np, 'loadtxt', mock_loadtxt)
    result = get_founder_info()
    expected = OrderedDict({'A': '#FF0000', 'B': '#00FF00'})
    assert result == expected


def test_unit_vector():
    vector = np.array([1, 2, 3])
    result = unit_vector(vector)
    expected = vector / np.linalg.norm(vector)
    assert np.allclose(result, expected)


def test_print_vecs(capsys):
    vecs = np.array([[1, 2], [3, 4]])
    print_vecs(vecs)
    captured = capsys.readouterr()
    assert captured.out == "       1.0        2.0\n\n       3.0        4.0\n\n"


def test_get_genotype_probability():
    aln_profile = np.array([1, 2, 3])
    aln_specificity = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]
    result = get_genotype_probability(aln_profile, aln_specificity)
    assert result.shape == (6,)


def test_ris_step():
    result = ris_step('AA', 'AA', 0.1)
    assert isinstance(result, float)


def test_get_transition_prob(tmp_path, monkeypatch):
    marker_file = tmp_path / "marker.txt"
    marker_file.write_text("marker1\tchr1\t1\t0.1\nmarker2\tchr1\t2\t0.2\n")
    output_file = tmp_path / "ref.gene_pos.ordered.npz"

    def mock_savez_compressed(file, **kwargs):
        assert os.path.basename(file) == os.path.basename(output_file)

    monkeypatch.setattr(np, 'savez_compressed', mock_savez_compressed)
    get_transition_prob(str(marker_file), output_file=str(output_file))

# note: there are additonal functions that do not have testing routines.
# get_alignment_spec, recontruct, interpolate, plot, and export.
# These functions require inputs that are not easily mocked, and can be quite large.
