"""Comparison of two conformations of the same protein on their CA trace:
contacts by cutoff, contacts exclusive to each, contact maps, and CA
angles/dihedrals. Exercise 7 runs it on adenylate kinase, exercise 8 on SERCA.
"""

import json

import matplotlib.pyplot as plt
import numpy as np

from contacts import contact_pairs, cutoff_map, distance_matrix
from geometry import ca_angles, ca_dihedrals
from plots import dihedral_panel, difference_map, secondary_structure_bands
from structures import OUT, ca_atoms, load, secondary_structure

# Parámetros de contacto (CA): explícitos y guardados junto con los resultados.
CUTOFF = 8.0          # A
MIN_SEQ_SEP = 4       # |i-j| >= 4: excluye los vecinos de secuencia y la vuelta de hélice i,i+3
# Contactos exclusivos de A: están en el mapa de A y no en el de B.
# Contactos "solo de A" según Whitford et al. 2007 (Models and Methods y Fig. 6,
# donde A es la forma cerrada y el conjunto se llama QLigand): exclusivos de A
# y, además, con la distancia CA-CA en B más de un 50% mayor que en A. Descarta
# los pares que solo cruzaron el umbral del cutoff sin separarse de verdad.
PAPER_MIN_RATIO = 1.5


def common_cas(chain_a, chain_b):
    """CA atoms of the residues present in both chains, matched by residue
    number (same protein, same numbering), plus those numbers."""
    a = {ca.get_parent().id[1]: ca for ca in ca_atoms(chain_a)}
    b = {ca.get_parent().id[1]: ca for ca in ca_atoms(chain_b)}
    resseq = sorted(a.keys() & b.keys())
    return [a[n] for n in resseq], [b[n] for n in resseq], np.array(resseq)


def exclusive(map_a, map_b):
    """Contacts present in A and absent in B."""
    return map_a & ~map_b


def paper_exclusive(map_a, dist_a, map_b, dist_b, min_ratio=PAPER_MIN_RATIO):
    """Exclusive to A and more than (min_ratio - 1) * 100 % further apart in B."""
    return exclusive(map_a, map_b) & (dist_b > min_ratio * dist_a)


def compare(name_a, name_b, chain_id="A", figs_dir=None, prefix="", domains=None, notes=(), labels=None):
    """Run the whole comparison and return a dict with every result. `labels`
    maps a PDB code to a short description shown next to it in printouts and
    figures (file names and JSON keys keep the bare code)."""
    lab_a, lab_b = (f"{n} ({labels[n]})" if labels and n in labels else n for n in (name_a, name_b))
    chain_a, chain_b = load(name_a)[0][chain_id], load(name_b)[0][chain_id]
    cas_a, cas_b, resseq = common_cas(chain_a, chain_b)
    xyz_a = np.array([a.get_coord() for a in cas_a])
    xyz_b = np.array([a.get_coord() for a in cas_b])
    dist_a, dist_b = distance_matrix(xyz_a), distance_matrix(xyz_b)
    map_a = cutoff_map(xyz_a, resseq, CUTOFF, MIN_SEQ_SEP)
    map_b = cutoff_map(xyz_b, resseq, CUTOFF, MIN_SEQ_SEP)
    only_a, only_b = exclusive(map_a, map_b), exclusive(map_b, map_a)
    paper_a = paper_exclusive(map_a, dist_a, map_b, dist_b)
    paper_b = paper_exclusive(map_b, dist_b, map_a, dist_a)

    pairs_a, pairs_b = contact_pairs(map_a, resseq), contact_pairs(map_b, resseq)
    print(f"residuos en común ({chain_id}): {len(resseq)}  "
          f"[{lab_a}: {len(ca_atoms(chain_a))}, {lab_b}: {len(ca_atoms(chain_b))}]")
    print(f"contactos CA (d <= {CUTOFF} A, |i-j| >= {MIN_SEQ_SEP}):")
    print(f"  {lab_a}: {len(pairs_a)}   {lab_b}: {len(pairs_b)}   compartidos: {len(pairs_a & pairs_b)}")
    print(f"  exclusivos de {lab_a}: {only_a.sum() // 2}   exclusivos de {lab_b}: {only_b.sum() // 2}")
    print(f"  según el trabajo (exclusivos y > {PAPER_MIN_RATIO:g} veces más lejos en la otra): "
          f"solo {lab_a}: {paper_a.sum() // 2}   solo {lab_b}: {paper_b.sum() // 2}")

    angles_a, angles_b = ca_angles(cas_a), ca_angles(cas_b)
    dihedrals_a, dihedrals_b = ca_dihedrals(cas_a), ca_dihedrals(cas_b)
    print(f"ángulos CA: {len(angles_a)} / {len(angles_b)}   diedros CA: {len(dihedrals_a)} / {len(dihedrals_b)}")
    common_d = sorted(dihedrals_a.keys() & dihedrals_b.keys())
    delta = np.array([_angle_diff(dihedrals_a[k], dihedrals_b[k]) for k in common_d])
    top = np.argsort(-np.abs(delta))[:8]
    print("diedros que más cambian (residuo: Δ grados): "
          + ", ".join(f"{common_d[k]}: {delta[k]:+.0f}" for k in top))

    result = {
        "parametros": {"cutoff": CUTOFF, "min_seq_sep": MIN_SEQ_SEP,
                       "paper_min_ratio": PAPER_MIN_RATIO, "cadena": chain_id},
        "resseq": resseq.tolist(),
        "contactos": {name_a: sorted(map(list, pairs_a)), name_b: sorted(map(list, pairs_b))},
        "exclusivos": {name_a: sorted(map(list, contact_pairs(only_a, resseq))),
                       name_b: sorted(map(list, contact_pairs(only_b, resseq)))},
        "solo_segun_trabajo": {name_a: sorted(map(list, contact_pairs(paper_a, resseq))),
                               name_b: sorted(map(list, contact_pairs(paper_b, resseq)))},
        "angulos": {name_a: angles_a, name_b: angles_b},
        "diedros": {name_a: dihedrals_a, name_b: dihedrals_b},
    }
    path = OUT / f"{prefix}{name_a}_{name_b}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, default=int, indent=1))
    print(f"contactos, ángulos y diedros guardados en data/out/{path.name}")

    if figs_dir is not None:
        _figures(figs_dir, f"{prefix}{name_a}_{name_b}", name_a, name_b, resseq, map_a, map_b,
                 angles_a, angles_b, dihedrals_a, dihedrals_b, chain_id, domains, notes, (lab_a, lab_b))
    return result


def _angle_diff(x, y):
    """Signed difference x - y wrapped into (-180, 180]."""
    return (x - y + 180) % 360 - 180


def _figures(figs_dir, stem, name_a, name_b, resseq, map_a, map_b,
             angles_a, angles_b, dih_a, dih_b, chain_id="A", domains=None, notes=(), shown=None):
    """One file per plot: contact map, CA angles, and CA dihedrals of each structure.
    `shown` are the names to display (code plus description)."""
    shown_a, shown_b = shown or (name_a, name_b)
    def save(fig, suffix):
        path = figs_dir / f"{stem}_{suffix}.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"figura: exercises/{path.name}")

    fig, ax = plt.subplots(figsize=(10, 9.5), layout="constrained")
    difference_map(ax, map_a, map_b, resseq, shown_a, shown_b,
                   f"mapa de contactos CA (≤{CUTOFF:g} Å, |i−j|≥{MIN_SEQ_SEP}): {shown_a} vs {shown_b}",
                   domains=domains, notes=notes)
    save(fig, "mapa_contactos")

    fig, ax = plt.subplots(figsize=(14, 4.5), layout="constrained")
    secondary_structure_bands(ax, secondary_structure(name_a, chain_id))
    for name, d in ((shown_a, angles_a), (shown_b, angles_b)):
        ax.plot(list(d), list(d.values()), lw=0.8, label=name)
    ax.set(xlabel="residuo (CA central)", ylabel="ángulo CA(i−1)-CA(i)-CA(i+1) [°]",
           title=f"ángulos CA (fondo: hélice rojo, hebra β verde, según HELIX/SHEET de {name_a})")
    ax.legend(fontsize=8)
    save(fig, "angulos")

    for name, label, dih in ((name_a, shown_a, dih_a), (name_b, shown_b, dih_b)):
        fig, ax = plt.subplots(figsize=(14, 5), layout="constrained")
        dihedral_panel(ax, dih, secondary_structure(name, chain_id),
                       f"diedros CA de {label} (eje circular cortado en −40°: ±180° en el medio)")
        save(fig, f"diedros_{name}")
