"""Ejercicio 4: geometría del O2 unido al hierro del hemo.

"N" en el diedro N-Fe-O-O es ambiguo. Se reportan los cuatro nitrógenos
pirrólicos del hemo (NA..ND): el diedro con cada uno da la orientación
azimutal del O2 sobre el plano del hemo. El NE2 de la histidina proximal no
sirve para esto: NE2-Fe-O es casi colineal, y un diedro cuyo ángulo central
es ~180 grados está indefinido (el plano NE2-Fe-O no existe).
"""

import math

from Bio.PDB.vectors import calc_angle, calc_dihedral

from structures import fetch, load

# (pdb, cadena). La His proximal de cada sitio no se escribe a mano: se lee
# de los registros LINK del archivo (ver proximal_his).
SITES = [("1A6M", "A"),    # oximioglobina de cachalote
         ("1HHO", "A"),    # oxihemoglobina humana, alfa
         ("1HHO", "B")]    # oxihemoglobina humana, beta

PYRROLE_N = ["NA", "NB", "NC", "ND"]


def residue_named(chain, name):
    [residue] = [r for r in chain if r.get_resname() == name]
    return residue


def degrees(f, *atoms):
    return math.degrees(f(*(a.get_vector() for a in atoms)))


def links(pdb_id):
    """Bonds declared by the authors in the LINK records, which Biopython
    does not parse. Each one as ((atom, resname, chain, resseq), same, distance)."""
    def side(line, start):
        return (line[start:start + 4].strip(), line[start + 5:start + 8].strip(),
                line[start + 9], int(line[start + 10:start + 14]))
    return [(side(line, 12), side(line, 42), float(line[73:78]))
            for line in fetch(pdb_id).read_text().splitlines() if line.startswith("LINK")]


def linked_to_iron(pdb_id, chain_id, resname, atom=None):
    """(resseq, atom name, distance) of every atom of residue type `resname`
    that LINK declares bonded to the heme iron of the given chain."""
    out = []
    for a, b, dist in links(pdb_id):
        for fe, other in ((a, b), (b, a)):
            if fe[:3] == ("FE", "HEM", chain_id) and other[1] == resname and other[2] == chain_id:
                if atom is None or other[0] == atom:
                    out.append((other[3], other[0], dist))
    return out


def proximal_his(pdb_id, chain):
    """The histidine bonded to the iron, read from LINK and cross-checked by
    geometry: its NE2 must be the closest histidine NE2 to the iron."""
    [(his, _, _)] = linked_to_iron(pdb_id, chain.id, "HIS", "NE2")
    fe = residue_named(chain, "HEM")["FE"]
    closest = min((r for r in chain if r.get_resname() == "HIS" and "NE2" in r),
                  key=lambda r: r["NE2"] - fe)
    assert closest.id[1] == his, f"LINK dice His{his}, la geometría His{closest.id[1]}"
    return his


def run():
    for pdb_id, chain_id in SITES:
        chain = load(pdb_id)[0][chain_id]
        heme, oxy = residue_named(chain, "HEM"), residue_named(chain, "OXY")
        fe = heme["FE"]
        # El O unido es el más cercano al Fe; los nombres O1/O2 del archivo no lo garantizan.
        o_bound, o_far = sorted(oxy, key=lambda o: o - fe)
        his = proximal_his(pdb_id, chain)
        ne2 = chain[his]["NE2"]

        # Control: la distancia Fe-O calculada tiene que coincidir con la que
        # los autores declararon en LINK para ese mismo oxígeno (LINK tiene 2 decimales).
        [(_, _, link_dist)] = linked_to_iron(pdb_id, chain_id, "OXY", o_bound.get_id())
        assert abs((fe - o_bound) - link_dist) < 0.01, (fe - o_bound, link_dist)

        print(f"== {pdb_id} cadena {chain_id}  (O unido = {o_bound.get_id()}, His proximal = His{his} según LINK)")
        print(f"  d(Fe-O)            {fe - o_bound:6.2f} A   (LINK: {link_dist:.2f})")
        print(f"  d(O-O)             {o_bound - o_far:6.2f} A")
        print(f"  ang(Fe-O-O)        {degrees(calc_angle, fe, o_bound, o_far):6.1f} grados")
        print(f"  ang(NE2-Fe-O)      {degrees(calc_angle, ne2, fe, o_bound):6.1f} grados   (His{his}, ~180 = colineal)")
        dihedrals = [degrees(calc_dihedral, heme[n], fe, o_bound, o_far) for n in PYRROLE_N]
        for n, value in zip(PYRROLE_N, dihedrals):
            print(f"  diedro({n}-Fe-O-O)  {value:6.1f} grados")

        # Control: NA..ND rodean al Fe separados ~90 grados entre sí, así que los
        # cuatro diedros tienen que diferir en ~90 grados de uno al siguiente.
        steps = [abs((b - a + 180) % 360 - 180) for a, b in zip(dihedrals, dihedrals[1:] + dihedrals[:1])]
        print(f"  saltos entre diedros consecutivos: {[round(x) for x in steps]}")
        assert all(abs(x - 90) < 20 for x in steps), steps
