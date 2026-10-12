"""Pseudo-bond angles and dihedrals of a CA trace, and superposition helpers."""

import math

import numpy as np
from Bio.PDB import Superimposer
from Bio.PDB.vectors import calc_angle, calc_dihedral

# Consecutive CAs sit 3.8 A apart (trans peptide). Anything well above that
# means residues are missing in between, and an angle across the gap would
# mix two unrelated stretches of chain.
MAX_CA_CA = 4.2


def _consecutive(cas):
    return all(a - b <= MAX_CA_CA for a, b in zip(cas, cas[1:]))


def ca_angles(cas):
    """Angle CA(i)-CA(i+1)-CA(i+2) in degrees, keyed by the resseq of the
    middle CA. Windows that span a chain break are skipped."""
    out = {}
    for w in zip(cas, cas[1:], cas[2:]):
        if _consecutive(w):
            out[w[1].get_parent().id[1]] = math.degrees(calc_angle(*(a.get_vector() for a in w)))
    return out


def ca_dihedrals(cas):
    """Dihedral CA(i)..CA(i+3) in degrees, keyed by the resseq of CA(i+1)
    (the first atom of the central bond)."""
    out = {}
    for w in zip(cas, cas[1:], cas[2:], cas[3:]):
        if _consecutive(w):
            out[w[1].get_parent().id[1]] = math.degrees(calc_dihedral(*(a.get_vector() for a in w)))
    return out


def superimpose(ref_atoms, mobile_atoms, mobile_entity):
    """Fit mobile_atoms onto ref_atoms and move the whole mobile_entity with
    the resulting rotation/translation. Returns the Superimposer."""
    sup = Superimposer()
    sup.set_atoms(ref_atoms, mobile_atoms)
    sup.apply(list(mobile_entity.get_atoms()))
    return sup


def rmsd(atoms_a, atoms_b):
    """RMSD computed straight from coordinates, as an independent check on
    what Superimposer reports."""
    a = np.array([x.get_coord() for x in atoms_a])
    b = np.array([x.get_coord() for x in atoms_b])
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))
