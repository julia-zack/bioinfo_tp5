"""Invariants of the contact maps and the geometry helpers, on hand-built
coordinates whose answer is known without running the code under test.

    venv/bin/python -m pytest tests
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from contacts import cutoff_map, neighbor_counts, occluders, shadow_map  # noqa: E402

RNG = np.random.default_rng(20261008)


def cloud(n=60, box=15.0):
    return RNG.uniform(0, box, size=(n, 3)), np.arange(1, n + 1)


# --- occlusion test on three atoms ------------------------------------------

def test_atom_on_the_segment_occludes():
    coords = np.array([[0, 0, 0], [3, 0.5, 0], [6, 0, 0]], float)
    assert list(occluders(coords, 0, 2, shadow_radius=1.0)) == [1]


def test_atom_off_to_the_side_does_not_occlude():
    coords = np.array([[0, 0, 0], [3, 2.0, 0], [6, 0, 0]], float)
    assert occluders(coords, 0, 2, shadow_radius=1.0).size == 0


@pytest.mark.parametrize("x", [-1.0, 7.0])
def test_atom_beyond_the_endpoints_does_not_occlude(x):
    # On the line i-j but outside the segment: t < 0 or t > 1.
    coords = np.array([[0, 0, 0], [x, 0, 0], [6, 0, 0]], float)
    assert occluders(coords, 0, 2, shadow_radius=1.0).size == 0


def test_perpendicular_distance_matches_hand_computation():
    # k at (2, 0.9, 0) is 0.9 from the x axis: occludes for S=1, not for S=0.8.
    coords = np.array([[0, 0, 0], [2, 0.9, 0], [6, 0, 0]], float)
    assert occluders(coords, 0, 2, 1.0).size == 1
    assert occluders(coords, 0, 2, 0.8).size == 0


# --- whole maps ---------------------------------------------------------------

def test_cutoff_map_symmetric_and_respects_sequence_separation():
    coords, resseq = cloud()
    m = cutoff_map(coords, resseq, cutoff=8.0, min_seq_sep=4)
    assert (m == m.T).all()
    for k in range(4):
        assert not np.diagonal(m, k).any()


def test_cutoff_map_matches_brute_force():
    coords, resseq = cloud(30)
    m = cutoff_map(coords, resseq, 8.0, 3)
    for i in range(30):
        for j in range(30):
            d = math.dist(coords[i], coords[j])
            assert m[i, j] == (d <= 8.0 and abs(int(resseq[i]) - int(resseq[j])) >= 3)


def test_shadow_with_zero_radius_equals_cutoff():
    coords, resseq = cloud()
    assert (shadow_map(coords, resseq, 8.0, 0.0, 4) == cutoff_map(coords, resseq, 8.0, 4)).all()


def test_shadow_is_symmetric_subset_and_monotone_in_S():
    coords, resseq = cloud()
    cut = cutoff_map(coords, resseq, 8.0, 4)
    previous = cut
    for s in (0.5, 1.0, 2.0, 3.0):
        sh = shadow_map(coords, resseq, 8.0, s, 4)
        assert (sh == sh.T).all()
        assert not (sh & ~cut).any()          # never invents a contact
        assert not (sh & ~previous).any()     # bigger S only removes
        previous = sh


def test_sequence_separation_uses_residue_numbers_not_positions():
    # Residues 10 and 20 are adjacent in the array but 10 apart in sequence.
    coords = np.array([[0, 0, 0], [3.8, 0, 0]], float)
    assert cutoff_map(coords, np.array([10, 20]), 8.0, 4)[0, 1]


# --- neighbour counts ---------------------------------------------------------

def test_neighbor_counts_exclude_self_and_split_by_class():
    coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [50, 0, 0]], float)
    counts = neighbor_counts(coords, "DLKS", 5.0, {"polar": set("DKS"), "nonpolar": set("L")})
    # Residue 0 (D) sees L and K, not itself and not the far S.
    assert counts["polar"][0] == 1 and counts["nonpolar"][0] == 1
    assert counts["polar"][3] == 0 and counts["nonpolar"][3] == 0
