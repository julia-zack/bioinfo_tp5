"""Loading invariants that need real files (downloaded to data/ on first use).

    venv/bin/python -m pytest tests
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from structures import amino_residues, load, mixed_altloc_residues  # noqa: E402


@pytest.mark.parametrize("pdb_id", ["1A6M", "3RGK", "1LNI", "1AKE"])
def test_no_residue_mixes_alternative_locations(pdb_id):
    # Biopython alone builds 1LNI Pro27 from two conformations (see structures.py).
    assert mixed_altloc_residues(load(pdb_id)) == []


@pytest.mark.parametrize("pdb_id", ["1A6M", "1LNI", "1LSE"])
def test_every_ca_is_a_carbon_bonded_to_n_and_c(pdb_id):
    for residue in amino_residues(load(pdb_id)[0]["A"]):
        ca = residue["CA"]
        assert ca.element == "C"
        assert ca - residue["N"] < 1.6 and ca - residue["C"] < 1.6, residue
