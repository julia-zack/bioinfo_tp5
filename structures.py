"""Download, cache and parse PDB/mmCIF files, plus the residue predicates the
exercises keep needing.

Files land in data/ the first time they are requested and are read from
there afterwards, so the exercises run offline once everything is cached.
"""

import urllib.request
import warnings
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.SeqUtils import seq1

DATA = Path(__file__).parent / "data"
OUT = DATA / "out"
RCSB = "https://files.rcsb.org/download/{}.{}"

BACKBONE = {"N", "CA", "C", "O", "OXT"}


def fetch(pdb_id, fmt="pdb"):
    """Path to data/<ID>.<fmt>, downloading it from RCSB if it is not there."""
    pdb_id = pdb_id.upper()
    path = DATA / f"{pdb_id}.{fmt}"
    if not path.exists():
        DATA.mkdir(exist_ok=True)
        urllib.request.urlretrieve(RCSB.format(pdb_id, fmt), path)
    return path


def load(pdb_id, fmt="pdb"):
    """Parse a structure. PERMISSIVE is the parser default; the warnings it
    raises about small format irregularities are silenced on purpose."""
    path = fetch(pdb_id, fmt)
    parser = PDBParser(PERMISSIVE=1, QUIET=True) if fmt == "pdb" else MMCIFParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = parser.get_structure(pdb_id.upper(), path)
    select_consistent_altlocs(structure)
    return structure


def select_consistent_altlocs(structure):
    """Make every residue use a single alternative location for all its atoms.

    Biopython picks the altloc of each disordered atom independently (highest
    occupancy), which can assemble a residue from two different conformations:
    in 1LNI Pro27, N and C come from B (occupancy 0.51) and CA from A (a 0.49
    tie), giving a 0.86 A CA-C "bond". Here each residue gets the altloc with
    the highest total occupancy over its atoms; ties go to the first label.
    """
    for residue in structure.get_residues():
        atoms = [a for a in residue if a.is_disordered()]
        if not atoms:
            continue
        occupancy = {}
        for atom in atoms:
            for altloc, child in atom.child_dict.items():
                occupancy[altloc] = occupancy.get(altloc, 0) + child.get_occupancy()
        best = max(sorted(occupancy), key=occupancy.get)
        for atom in atoms:
            if atom.disordered_has_id(best):
                atom.disordered_select(best)


def mixed_altloc_residues(structure):
    """Residues whose currently selected atoms come from more than one altloc."""
    return [r for r in structure.get_residues()
            if len({a.get_altloc() for a in r if a.get_altloc() != " "}) > 1]


def out_path(name):
    OUT.mkdir(parents=True, exist_ok=True)
    return OUT / name


# The hetero field of residue.id is " " for standard residues, "W" for water
# and "H_<resname>" for any other hetero group.

def is_amino(residue):
    return residue.id[0] == " "


def is_water(residue):
    return residue.id[0] == "W"


def is_ligand(residue):
    return residue.id[0].startswith("H_")


def amino_residues(chain):
    return [r for r in chain if is_amino(r)]


def sequence(chain):
    """One-letter sequence of the standard residues of a chain, in file order."""
    return "".join(seq1(r.get_resname()) for r in amino_residues(chain))


def ca_atoms(chain):
    """CA of every standard residue that has one. Filtering by residue type
    matters: a calcium ion is also an atom called 'CA'."""
    return [r["CA"] for r in amino_residues(chain) if "CA" in r]


def secondary_structure(pdb_id, chain_id="A"):
    """Residue number -> 'hélice' / 'hebra β' from the HELIX and SHEET records."""
    out = {}
    for line in fetch(pdb_id).read_text().splitlines():
        if line.startswith("HELIX") and line[19] == chain_id:
            out.update({n: "hélice" for n in range(int(line[21:25]), int(line[33:37]) + 1)})
        elif line.startswith("SHEET") and line[21] == chain_id:
            out.update({n: "hebra β" for n in range(int(line[22:26]), int(line[33:37]) + 1)})
    return out

