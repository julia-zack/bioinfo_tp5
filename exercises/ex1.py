"""Ejercicio 1: cargar la misma estructura en PDB y en mmCIF y explorar el header."""

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from structures import fetch, load

PDB_ID = "1A6M"  # oximioglobina de cachalote, 1.0 A

HEADER_KEYS = ["name", "head", "deposition_date", "release_date", "structure_method",
               "resolution", "author", "has_missing_residues", "missing_residues"]

CIF_KEYS = ["_struct.title", "_exptl.method", "_refine.ls_d_res_high",
            "_entity.pdbx_description", "_entity_src_nat.pdbx_organism_scientific",
            "_symmetry.space_group_name_H-M", "_cell.length_a", "_cell.length_b",
            "_cell.length_c"]


def counts(structure):
    return {"modelos": len(structure),
            "cadenas": sum(1 for _ in structure.get_chains()),
            "residuos": sum(1 for _ in structure.get_residues()),
            "átomos": sum(1 for _ in structure.get_atoms())}


def run():
    pdb = load(PDB_ID, "pdb")
    cif = load(PDB_ID, "cif")

    print(f"== header del PDB ({PDB_ID}.pdb)")
    for key in HEADER_KEYS:
        print(f"  {key:22} {pdb.header.get(key)}")
    for mol, info in pdb.header["compound"].items():
        print(f"  compound[{mol}]       {info}")
    print(f"  source[1]              {pdb.header['source']['1'].get('organism_scientific')}")

    print(f"\n== header del mmCIF ({PDB_ID}.cif) vía MMCIFParser")
    for key in HEADER_KEYS[:6]:
        print(f"  {key:22} {cif.header.get(key)}")

    mmcif = MMCIF2Dict(fetch(PDB_ID, "cif"))
    print(f"\n== MMCIF2Dict: {len(mmcif)} claves; algunas que el header del PDB no expone")
    for key in CIF_KEYS:
        print(f"  {key:42} {mmcif.get(key)}")

    print("\n== ¿los dos formatos dan el mismo objeto?")
    a, b = counts(pdb), counts(cif)
    for key in a:
        print(f"  {key:9} pdb={a[key]:6}  cif={b[key]:6}  {'ok' if a[key] == b[key] else 'DIFIEREN'}")
