"""Ejercicio 2: hemo, conteos, secuencias y un PDB con solo los CA."""

import re
from collections import Counter

from Bio.PDB import PDBIO, Select
from Bio.SeqUtils import seq1

from structures import (amino_residues, ca_atoms, fetch, is_amino, is_ligand, is_water, load,
                        out_path, sequence)

MB = "1A6M"   # oximioglobina, una cadena
HB = "1HHO"   # oxihemoglobina, dos cadenas (alfa y beta): para que "cada cadena" signifique algo


def formula(pdb_id, resname):
    """Element counts of a hetero group, read from the FORMUL record of the
    PDB file (Biopython does not parse it). E.g. 'C34 H32 FE N4 O4'."""
    for line in fetch(pdb_id).read_text().splitlines():
        if line.startswith("FORMUL") and line[12:15].strip() == resname:
            return {re.match(r"[A-Z]+", t)[0]: int(re.sub(r"[A-Z]+", "", t) or 1)
                    for t in line[19:].split()}
    raise ValueError(f"no FORMUL record for {resname} in {pdb_id}")


def a():
    structure = load(MB)

    # Dos maneras independientes de encontrar el hemo, que tienen que coincidir:
    # por su código en el PDB, y por contenido (la única molécula no proteica con hierro).
    by_name = [r for r in structure.get_residues() if r.get_resname() == "HEM"]
    by_iron = [r for r in structure.get_residues() if is_ligand(r) and "FE" in r]
    print(f"por nombre (HEM):          {[r.full_id for r in by_name]}")
    print(f"por contenido de hierro:   {[r.full_id for r in by_iron]}")
    assert by_name == by_iron, "el nombre y el contenido de hierro no señalan el mismo residuo"

    [heme] = by_name
    print(f"\nHEM  id={heme.id}  full_id={heme.full_id}")

    # Tercer control: los átomos del residuo tienen que coincidir con la fórmula
    # declarada en el archivo, sin contar hidrógenos (los rayos X casi no los resuelven).
    expected = {el: n for el, n in formula(MB, "HEM").items() if el != "H"}
    found = Counter(atom.element for atom in heme if atom.element != "H")
    print(f"fórmula del archivo (sin H): {expected}  -> {sum(expected.values())} átomos")
    print(f"átomos en el residuo (sin H): {dict(found)}  -> {sum(found.values())} átomos")
    assert found == Counter(expected), "los átomos del residuo no coinciden con la fórmula"

    print("\nResiduos que no son proteína ni agua:")
    for name, n in Counter(r.get_resname() for r in structure.get_residues() if is_ligand(r)).items():
        print(f"  {name}: {n}")


def count(structure):
    model = structure[0]
    residues = list(model.get_residues())
    atoms = list(model.get_atoms())
    by_kind = Counter("proteína" if is_amino(a.get_parent()) else
                      "agua" if is_water(a.get_parent()) else "hetero" for a in atoms)
    return {
        "cadenas": len(model),
        "residuos (todos)": len(residues),
        "residuos de proteína": sum(map(is_amino, residues)),
        "átomos totales": len(atoms),
        "átomos de proteína": by_kind["proteína"],
        # En el formato PDB las aguas también son HETATM; se muestran ambos números.
        "HETATM (incluye agua)": by_kind["hetero"] + by_kind["agua"],
        "HETATM (sin agua)": by_kind["hetero"],
        "átomos de agua": by_kind["agua"],
    }


def raw_atom_counts(pdb_id):
    """Count atoms straight from the text of the file, without Biopython.
    An atom with alternative positions (altloc) appears in several lines;
    it is counted once, by (chain, residue number + insertion code, atom name)."""
    seen = {"ATOM": set(), "HETATM": set(), "HOH": set()}
    lines = Counter()
    for line in fetch(pdb_id).read_text().splitlines():
        record = line[:6].strip()
        if record not in ("ATOM", "HETATM"):
            continue
        lines[record] += 1
        atom = (line[21], line[22:27], line[12:16].strip())
        seen[record].add(atom)
        if line[17:20] == "HOH":
            seen["HOH"].add(atom)
    return {k: len(v) for k, v in seen.items()}, lines


def b():
    for pdb_id in (MB, HB):
        print(f"== {pdb_id} (modelo 0)")
        result = count(load(pdb_id))
        for key, value in result.items():
            print(f"  {key:24} {value}")
        assert result["átomos totales"] == (result["átomos de proteína"]
                                            + result["HETATM (incluye agua)"])

        # Control independiente: los mismos números contados sobre el texto del archivo.
        raw, lines = raw_atom_counts(pdb_id)
        print(f"  texto crudo: {lines['ATOM']} líneas ATOM y {lines['HETATM']} HETATM; "
              f"sin repetir altlocs: {raw['ATOM']} / {raw['HETATM']} / {raw['HOH']} de agua")
        assert (raw["ATOM"], raw["HETATM"], raw["HOH"]) == (
            result["átomos de proteína"], result["HETATM (incluye agua)"], result["átomos de agua"])


def seqres(pdb_id):
    """Sequence of each chain as declared in the SEQRES records: what the
    authors put in the experiment, whether or not it was seen in the density."""
    chains = {}
    for line in fetch(pdb_id).read_text().splitlines():
        if line.startswith("SEQRES"):
            chains.setdefault(line[11], []).extend(line[19:].split())
    return {chain: "".join(seq1(r) for r in names) for chain, names in chains.items()}


def c():
    for pdb_id in (MB, HB):
        declared = seqres(pdb_id)
        for chain in load(pdb_id)[0]:
            seq = sequence(chain)
            numbers = [r.id[1] for r in amino_residues(chain)]
            gaps = [(x, y) for x, y in zip(numbers, numbers[1:]) if y != x + 1]
            print(f">{pdb_id}_{chain.id} len={len(seq)} huecos_en_numeración={gaps or 'ninguno'}")
            print(seq)
            # Control independiente: la secuencia de las coordenadas contra la declarada en SEQRES.
            # Si difieren, faltan residuos sin resolver o hay aminoácidos modificados.
            print(f"  igual a SEQRES ({len(declared[chain.id])} aa): {seq == declared[chain.id]}")
            assert seq == declared[chain.id]


class OnlyProteinCA(Select):
    """Keeps the CA of standard residues. Checking the residue too keeps a
    calcium ion (also named 'CA') out."""

    def accept_residue(self, residue):
        return is_amino(residue)

    def accept_atom(self, atom):
        return atom.get_id() == "CA"


def d():
    structure = load(MB)
    path = out_path(f"{MB}_CA.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path), OnlyProteinCA())

    # Releer lo escrito y comparar contra lo esperado, en vez de confiar en el Select.
    written = load_written(path)
    expected = sum(len(ca_atoms(ch)) for ch in structure[0])
    names = {a.get_id() for a in written.get_atoms()}
    print(f"escrito {path.relative_to(path.parents[2])}")
    n_written = sum(1 for _ in written.get_atoms())
    print(f"átomos releídos={n_written} esperados={expected} nombres={names}")
    assert n_written == expected and names == {"CA"}


def load_written(path):
    from Bio.PDB import PDBParser
    return PDBParser(QUIET=True).get_structure("ca", path)
