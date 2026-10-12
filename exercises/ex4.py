"""Ejercicio 4: geometría del O2 unido al hierro del hemo.

"N" en el diedro N-Fe-O-O es ambiguo. Se reportan los cuatro nitrógenos
pirrólicos del hemo (NA..ND): el diedro con cada uno da la orientación
azimutal del O2 sobre el plano del hemo. El NE2 de la histidina proximal no
sirve para esto: NE2-Fe-O es casi colineal, y un diedro cuyo ángulo central
es ~180 grados está indefinido (el plano NE2-Fe-O no existe).
"""

import math

from Bio.PDB.vectors import calc_angle, calc_dihedral

from structures import load

# (pdb, cadena, nro de la His proximal)
SITES = [("1A6M", "A", 93),    # oximioglobina de cachalote
         ("1HHO", "A", 87),    # oxihemoglobina humana, alfa
         ("1HHO", "B", 92)]    # oxihemoglobina humana, beta

PYRROLE_N = ["NA", "NB", "NC", "ND"]


def residue_named(chain, name):
    [residue] = [r for r in chain if r.get_resname() == name]
    return residue


def degrees(f, *atoms):
    return math.degrees(f(*(a.get_vector() for a in atoms)))


def run():
    for pdb_id, chain_id, his in SITES:
        chain = load(pdb_id)[0][chain_id]
        heme, oxy = residue_named(chain, "HEM"), residue_named(chain, "OXY")
        fe = heme["FE"]
        # El O unido es el más cercano al Fe; los nombres O1/O2 del archivo no lo garantizan.
        o_bound, o_far = sorted(oxy, key=lambda o: o - fe)
        ne2 = chain[his]["NE2"]

        print(f"== {pdb_id} cadena {chain_id}  (O unido = {o_bound.get_id()})")
        print(f"  d(Fe-O)            {fe - o_bound:6.2f} A")
        print(f"  d(O-O)             {o_bound - o_far:6.2f} A")
        print(f"  ang(Fe-O-O)        {degrees(calc_angle, fe, o_bound, o_far):6.1f} grados")
        print(f"  ang(NE2-Fe-O)      {degrees(calc_angle, ne2, fe, o_bound):6.1f} grados   (His{his}, ~180 = colineal)")
        for n in PYRROLE_N:
            print(f"  diedro({n}-Fe-O-O)  {degrees(calc_dihedral, heme[n], fe, o_bound, o_far):6.1f} grados")
