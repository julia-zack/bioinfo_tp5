"""Ejercicio 5: superponer dos mioglobinas homólogas.

Referencia: 1A6M (cachalote). Móvil: 3RGK (humana, mutante K45R).
Las dos secuencias no tienen el mismo largo ni la misma numeración, así que
los CA equivalentes no se pueden emparejar por número de residuo: se
emparejan por un alineamiento de secuencia. PairwiseAligner se usa acá como
herramienta (el TP no pide implementar el alineamiento).
"""

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBIO

from geometry import rmsd, superimpose
from structures import amino_residues, load, out_path, sequence

REF, MOBILE, CHAIN = "1A6M", "3RGK", "A"


def aligner():
    al = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"))
    al.open_gap_score, al.extend_gap_score = -10, -0.5
    return al


def equivalent_cas(ref_chain, mob_chain):
    """CA pairs of residues aligned to each other (match or mismatch)."""
    ref_res = [r for r in amino_residues(ref_chain) if "CA" in r]
    mob_res = [r for r in amino_residues(mob_chain) if "CA" in r]
    aln = aligner().align(sequence(ref_chain), sequence(mob_chain))[0]
    pairs = []
    for (r0, r1), (m0, m1) in zip(*aln.aligned):
        pairs += zip(ref_res[r0:r1], mob_res[m0:m1])
    return aln, [(r["CA"], m["CA"]) for r, m in pairs]


def run():
    ref, mob = load(REF), load(MOBILE)
    aln, pairs = equivalent_cas(ref[0][CHAIN], mob[0][CHAIN])
    ref_ca, mob_ca = [p[0] for p in pairs], [p[1] for p in pairs]
    identical = sum(r.get_parent().get_resname() == m.get_parent().get_resname() for r, m in pairs)

    print(aln)
    print(f"pares de CA equivalentes: {len(pairs)}  identidad: {identical / len(pairs):.1%}")
    print(f"RMSD antes de superponer (coordenadas crudas): {rmsd(ref_ca, mob_ca):7.2f} A")

    sup = superimpose(ref_ca, mob_ca, mob)
    rotation, translation = sup.rotran
    print(f"RMSD reportado por Superimposer:              {sup.rms:7.3f} A")
    # Chequeo independiente: recalcular desde las coordenadas ya movidas.
    print(f"RMSD recalculado tras aplicar la matriz:      {rmsd(ref_ca, mob_ca):7.3f} A")
    print(f"rotación:\n{rotation.round(4)}\ntraslación: {translation.round(3)}")

    path = out_path(f"{MOBILE}_on_{REF}.pdb")
    io = PDBIO()
    io.set_structure(mob)
    io.save(str(path))
    print(f"estructura móvil superpuesta guardada en data/out/{path.name}")
