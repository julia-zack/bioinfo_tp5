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
        return parser.get_structure(pdb_id.upper(), path)


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
