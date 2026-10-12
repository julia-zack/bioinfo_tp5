"""Residue representatives, distance matrices, contact maps (cutoff and
shadow) and neighbour counts.

Everything here works on plain numpy coordinate arrays so it can be tested
with hand-built geometries, without any PDB file.
"""

import numpy as np

from structures import BACKBONE, amino_residues

POLAR = set("STKRDEQN")
NONPOLAR = set("ILVWYFM")


# --- one point per residue ------------------------------------------------

def side_chain_center(residue):
    """Mass-weighted center of the side-chain heavy atoms. Glycine has no side
    chain, so it falls back to its CA (same convention as `representative`)."""
    atoms = [a for a in residue if a.get_id() not in BACKBONE and a.element != "H"]
    if not atoms:
        return residue["CA"].get_coord()
    masses = np.array([a.mass for a in atoms])
    coords = np.array([a.get_coord() for a in atoms])
    return masses @ coords / masses.sum()


def representative(residue, kind):
    """Point that stands for a residue: 'CA', 'CB' or 'CM' (side-chain center
    of mass). Residues without CB (glycine, or a truncated side chain) use CA."""
    if kind == "CA":
        return residue["CA"].get_coord()
    if kind == "CB":
        return (residue["CB"] if "CB" in residue else residue["CA"]).get_coord()
    if kind == "CM":
        return side_chain_center(residue)
    raise ValueError(f"unknown representative '{kind}'")


def residue_points(chain, kind):
    """(residues, coords) for the standard residues of a chain that have a CA."""
    residues = [r for r in amino_residues(chain) if "CA" in r]
    coords = np.array([representative(r, kind) for r in residues], dtype=float)
    return residues, coords


# --- contact maps ----------------------------------------------------------

def distance_matrix(coords):
    diff = coords[:, None, :] - coords[None, :, :]
    return np.sqrt((diff ** 2).sum(-1))


def separation_mask(resseq, min_seq_sep):
    """True where |resseq_i - resseq_j| >= min_seq_sep. Uses residue numbers,
    not array positions, so a missing residue does not shift the count."""
    resseq = np.asarray(resseq)
    return np.abs(resseq[:, None] - resseq[None, :]) >= min_seq_sep


def cutoff_map(coords, resseq, cutoff, min_seq_sep):
    """Boolean contact map: distance <= cutoff and far enough in sequence."""
    return (distance_matrix(coords) <= cutoff) & separation_mask(resseq, min_seq_sep)


def occluders(coords, i, j, shadow_radius):
    """Indices k that occlude the i-j contact, following the TP's recipe.

    Project k on the line i->j: t = ((k - i) . v) / (v . v) with v = j - i.
    t in (0, 1) means the foot of the perpendicular p = i + t v falls between
    i and j. k occludes when that happens and |k - p| < S.
    """
    v = coords[j] - coords[i]
    rel = coords - coords[i]
    t = rel @ v / (v @ v)
    perpendicular = np.linalg.norm(rel - t[:, None] * v, axis=1)
    hit = (t > 0) & (t < 1) & (perpendicular < shadow_radius)
    hit[[i, j]] = False
    return np.flatnonzero(hit)


def shadow_map(coords, resseq, cutoff, shadow_radius, min_seq_sep):
    """Cutoff map with every occluded contact removed.

    Every atom is a candidate occluder, not only the ones within the cutoff of
    i: an atom whose projection falls inside the segment and is within S of it
    is necessarily within sqrt(d_ij^2 + S^2) of i anyway, and taking all of
    them keeps the test symmetric in i and j.
    """
    contacts = cutoff_map(coords, resseq, cutoff, min_seq_sep)
    shadow = contacts.copy()
    for i, j in zip(*np.nonzero(np.triu(contacts))):
        if occluders(coords, i, j, shadow_radius).size:
            shadow[i, j] = shadow[j, i] = False
    return shadow


def contact_pairs(contact_map, resseq):
    """Set of (resseq_i, resseq_j) with i < j, for comparing two structures
    whose arrays may not line up position by position."""
    return {(resseq[i], resseq[j]) for i, j in zip(*np.nonzero(np.triu(contact_map)))}


# --- neighbour counts --------------------------------------------------------

def neighbor_counts(coords, letters, radius, classes):
    """For each residue, how many *other* residues of each class lie within
    `radius` of it. `classes` maps a name to a set of one-letter codes."""
    within = distance_matrix(coords) <= radius
    np.fill_diagonal(within, False)
    letters = np.asarray(list(letters))
    return {name: within[:, np.isin(letters, list(members))].sum(1)
            for name, members in classes.items()}
