"""Ejercicio 9: mapa de contactos "shadow" (Noel, Whitford y Onuchic 2012)
sobre los CA de 1AKE, con la recta de oclusión que describe la consigna.

Diferencia con el método original: en el trabajo, el átomo k proyecta una
sombra cónica (la luz sale del centro de i) y j queda tapado si cae dentro de
ese cono. La consigna lo simplifica a un cilindro de radio S alrededor del
segmento i-j: k ocluye si su proyección cae entre i y j y está a menos de S
del segmento. Es lo que se implementa en contacts.occluders.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from contacts import cutoff_map, shadow_map
from pairs import MIN_SEQ_SEP
from plots import kept_removed_map
from structures import ca_atoms, load

PDB_ID, CHAIN = "1AKE", "A"
SHADOW_RADII = [0.5, 1.0, 2.0, 3.0]   # A
CUTOFFS = [8.0, 10.0, 12.0]           # A
PLOT_CUTOFF = 10.0                    # con 8 A la sombra casi no saca nada (ver tabla)
FIGS = Path(__file__).parent


def ca_coords():
    cas = ca_atoms(load(PDB_ID)[0][CHAIN])
    return (np.array([a.get_coord() for a in cas]),
            np.array([a.get_parent().id[1] for a in cas]))


def run():
    coords, resseq = ca_coords()

    print(f"{PDB_ID} cadena {CHAIN}, {len(resseq)} CA, |i-j| >= {MIN_SEQ_SEP}")
    print("nº de contactos (pares)        " + "".join(f"S={s:<6g}" for s in SHADOW_RADII))
    for cutoff in CUTOFFS:
        n_cut = cutoff_map(coords, resseq, cutoff, MIN_SEQ_SEP).sum() // 2
        n_sh = [shadow_map(coords, resseq, cutoff, s, MIN_SEQ_SEP).sum() // 2 for s in SHADOW_RADII]
        print(f"  cutoff {cutoff:4g} A: {n_cut:5}  ->  " + "".join(f"{n:<8}" for n in n_sh))

    cutoff = PLOT_CUTOFF
    cut = cutoff_map(coords, resseq, cutoff, MIN_SEQ_SEP)
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    for ax, s in zip(axes.flat, SHADOW_RADII):
        sh = shadow_map(coords, resseq, cutoff, s, MIN_SEQ_SEP)
        # La sombra solo puede sacar contactos: si apareciera uno que no está en el corte simple, es un error.
        assert not (sh & ~cut).any(), "el mapa shadow tiene contactos que el de corte no tiene"
        kept_removed_map(ax, sh, cut, resseq,
                         f"S = {s:g} Å, C = {cutoff:g} Å: quedan {sh.sum() // 2} de {cut.sum() // 2}",
                         "quedan con la sombra", "eliminados por la sombra")
    fig.suptitle(f"{PDB_ID} (CA): contactos del corte por distancia que quedan (azul) "
                 f"o elimina (gris claro) el método de la sombra")
    fig.tight_layout()
    fig.savefig(FIGS / "ex9_shadow_vs_cutoff.png", dpi=150)
    print("figura: exercises/ex9_shadow_vs_cutoff.png")
