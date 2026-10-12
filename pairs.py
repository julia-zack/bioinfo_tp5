"""Comparison of two conformations of the same protein on their CA trace:
contacts by cutoff, contacts exclusive to each, contact maps, and CA
angles/dihedrals. Exercise 7 runs it on adenylate kinase, exercise 8 on SERCA.
"""

import json

import matplotlib.pyplot as plt
import numpy as np

from contacts import contact_pairs, cutoff_map, distance_matrix
from geometry import ca_angles, ca_dihedrals
from plots import difference_map, split_map
from structures import OUT, ca_atoms, load

# Parámetros de contacto (CA): explícitos y guardados junto con los resultados.
CUTOFF = 8.0          # A
MIN_SEQ_SEP = 4       # |i-j| >= 4: excluye los vecinos de secuencia y la vuelta de hélice i,i+3
# Un contacto "exclusivo de A" está en el mapa de A, no en el de B, y además
# la distancia en B supera a la de A en al menos EXCLUSIVE_MIN_DELTA. Con 0 es
# la diferencia de conjuntos pura; un valor > 0 descarta los pares que en B
# quedan apenas por encima del corte (no se separaron, solo cruzaron el umbral).
EXCLUSIVE_MIN_DELTA = 0.0


def common_cas(chain_a, chain_b):
    """CA atoms of the residues present in both chains, matched by residue
    number (same protein, same numbering), plus those numbers."""
    a = {ca.get_parent().id[1]: ca for ca in ca_atoms(chain_a)}
    b = {ca.get_parent().id[1]: ca for ca in ca_atoms(chain_b)}
    resseq = sorted(a.keys() & b.keys())
    return [a[n] for n in resseq], [b[n] for n in resseq], np.array(resseq)


def exclusive(map_a, dist_a, map_b, dist_b, min_delta):
    return map_a & ~map_b & (dist_b - dist_a >= min_delta)


def compare(name_a, name_b, chain_id="A", figs_dir=None, prefix=""):
    """Run the whole comparison and return a dict with every result."""
    chain_a, chain_b = load(name_a)[0][chain_id], load(name_b)[0][chain_id]
    cas_a, cas_b, resseq = common_cas(chain_a, chain_b)
    xyz_a = np.array([a.get_coord() for a in cas_a])
    xyz_b = np.array([a.get_coord() for a in cas_b])
    dist_a, dist_b = distance_matrix(xyz_a), distance_matrix(xyz_b)
    map_a = cutoff_map(xyz_a, resseq, CUTOFF, MIN_SEQ_SEP)
    map_b = cutoff_map(xyz_b, resseq, CUTOFF, MIN_SEQ_SEP)
    only_a = exclusive(map_a, dist_a, map_b, dist_b, EXCLUSIVE_MIN_DELTA)
    only_b = exclusive(map_b, dist_b, map_a, dist_a, EXCLUSIVE_MIN_DELTA)

    pairs_a, pairs_b = contact_pairs(map_a, resseq), contact_pairs(map_b, resseq)
    print(f"residuos en común ({chain_id}): {len(resseq)}  "
          f"[{name_a}: {len(ca_atoms(chain_a))},{name_b}: {len(ca_atoms(chain_b))}]")
    print(f"contactos CA (d <= {CUTOFF} A, |i-j| >= {MIN_SEQ_SEP}):")
    print(f"  {name_a}: {len(pairs_a)}   {name_b}: {len(pairs_b)}   compartidos: {len(pairs_a & pairs_b)}")
    print(f"  solo {name_a}: {only_a.sum() // 2}   solo {name_b}: {only_b.sum() // 2}"
          f"   (Δd mínimo = {EXCLUSIVE_MIN_DELTA} A)")

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
                       "exclusive_min_delta": EXCLUSIVE_MIN_DELTA, "cadena": chain_id},
        "resseq": resseq.tolist(),
        "contactos": {name_a: sorted(map(list, pairs_a)), name_b: sorted(map(list, pairs_b))},
        "exclusivos": {name_a: sorted(map(list, contact_pairs(only_a, resseq))),
                       name_b: sorted(map(list, contact_pairs(only_b, resseq)))},
        "angulos": {name_a: angles_a, name_b: angles_b},
        "diedros": {name_a: dihedrals_a, name_b: dihedrals_b},
    }
    path = OUT / f"{prefix}{name_a}_{name_b}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, default=int, indent=1))
    print(f"contactos, ángulos y diedros guardados en data/out/{path.name}")

    if figs_dir is not None:
        _figures(figs_dir / f"{prefix}{name_a}_{name_b}.png", name_a, name_b, resseq,
                 map_a, map_b, angles_a, angles_b, dihedrals_a, dihedrals_b)
    return result


def _angle_diff(x, y):
    """Signed difference x - y wrapped into (-180, 180]."""
    return (x - y + 180) % 360 - 180


def _figures(path, name_a, name_b, resseq, map_a, map_b, angles_a, angles_b, dih_a, dih_b):
    fig = plt.figure(figsize=(14, 9))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.4, 1])
    split_map(fig.add_subplot(grid[0, 0]), map_a, map_b, resseq, name_a, name_b,
              f"mapas de contacto CA (≤{CUTOFF:g} Å, |i−j|≥{MIN_SEQ_SEP})")
    difference_map(fig.add_subplot(grid[0, 1]), map_a, map_b, resseq, name_a, name_b,
                   "contactos compartidos y exclusivos")
    ax = fig.add_subplot(grid[1, 0])
    for name, d in ((name_a, angles_a), (name_b, angles_b)):
        ax.plot(list(d), list(d.values()), lw=0.8, label=name)
    ax.set(xlabel="residuo (CA central)", ylabel="ángulo CA(i−1)-CA(i)-CA(i+1) [°]")
    ax.legend(fontsize=8)
    ax = fig.add_subplot(grid[1, 1])
    for name, d in ((name_a, dih_a), (name_b, dih_b)):
        ax.plot(list(d), list(d.values()), ".", ms=2.5, label=name)
    ax.set(xlabel="residuo (CA i+1)", ylabel="diedro CA(i)…CA(i+3) [°]", ylim=(-180, 180))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"figura: exercises/{path.name}")
