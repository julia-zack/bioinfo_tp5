"""Ejercicio 8: lo mismo que el 7 para la Ca2+-ATPasa del retículo (SERCA1a
de conejo), los tres pares posibles entre:
  1SU4  E1·2Ca2+ (dos Ca unidos)
  4XOU  Ca2-E1·Mg-AMPPCP (análogo de ATP)
  3B9B  E2·BeF3- (análogo de fosfoenzima E2-P)
Misma proteína y misma numeración, así que los CA se emparejan por número."""

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt

from pairs import CUTOFF, PAPER_MIN_RATIO, compare
from structures import ca_atoms, load, rel

STRUCTURES = ["1SU4", "4XOU", "3B9B"]
# Apodos que se muestran junto al código en gráficos y salidas.
LABELS = {"1SU4": "con Ca", "4XOU": "Ca + ATP", "3B9B": "sin Ca"}

# Dominios de la SERCA1a según la bibliografía (no figuran en el archivo):
# A (actuador), P (fosforilación), N (unión de nucleótido) y M (membrana, 10 hélices).
DOMAINS = {"A": list(range(1, 41)) + list(range(124, 244)),
           "P": list(range(330, 360)) + list(range(605, 738)),
           "N": list(range(360, 605)),
           "M": list(range(41, 124)) + list(range(244, 330)) + list(range(738, 995))}


def run():
    for a, b in combinations(STRUCTURES, 2):
        print(f"\n######## {a} ({LABELS[a]}) vs {b} ({LABELS[b]})")
        result = compare(a, b, figs_dir=Path(__file__).parent, prefix="ex8_", domains=DOMAINS, labels=LABELS)
        contact_changes(a, b, result)


def domain_of(n):
    return next(d for d, residues in DOMAINS.items() if n in residues)


def contact_changes(name_a, name_b, result):
    """Like ex7_cambio_contactos, for each structure of the pair: one arrow per
    exclusive contact, from its CA-CA distance in the structure that has it
    (dot) to its distance in the other one (arrow head). With hundreds of
    contacts per set, the ones that meet the paper's criterion and the ones
    that do not go to separate figures, without residue labels, and the
    colour only says whether the pair is within one domain or between two."""
    ca = {p: {a.get_parent().id[1]: a for a in ca_atoms(load(p)[0]["A"])} for p in (name_a, name_b)}
    colors = {"mismo dominio": "0.55", "dominios distintos": "tab:orange"}
    for own, other in ((name_a, name_b), (name_b, name_a)):
        paper = {tuple(x) for x in result["solo_segun_trabajo"][own]}
        rows = [(ca[own][i] - ca[own][j], ca[other][i] - ca[other][j], i, j)
                for i, j in map(tuple, result["exclusivos"][own])]
        rows.sort(key=lambda r: r[1] / r[0])                      # least change at the bottom
        for meets, suffix, what in ((True, "cumplen", "cumplen"), (False, "no_cumplen", "NO cumplen")):
            subset = [r for r in rows if ((r[2], r[3]) in paper) == meets]
            # Every exclusive contact either passes the 1.5 x tick or does not; check the split.
            assert all((r[1] > PAPER_MIN_RATIO * r[0]) == meets for r in subset)
            _arrow_figure(subset, own, other, what, suffix, colors, name_a, name_b)


def _arrow_figure(rows, own, other, what, suffix, colors, name_a, name_b):
    fig, ax = plt.subplots(figsize=(10, 3 + 0.065 * len(rows)), layout="constrained")
    counts = {}
    for y, (d_own, d_other, i, j) in enumerate(rows):
        key = "mismo dominio" if domain_of(i) == domain_of(j) else "dominios distintos"
        counts[key] = counts.get(key, 0) + 1
        ax.annotate("", xy=(d_other, y), xytext=(d_own, y),
                    arrowprops=dict(arrowstyle="-|>", color=colors[key], lw=1, mutation_scale=6,
                                    shrinkA=0, shrinkB=0))
        ax.scatter(d_own, y, color=colors[key], s=6, zorder=3)
        ax.plot([PAPER_MIN_RATIO * d_own] * 2, [y - 0.4, y + 0.4], color="k", lw=0.6)   # 1.5 x threshold
    for key, color in colors.items():
        ax.scatter([], [], color=color, label=f"{key} ({counts.get(key, 0)})")
    ax.axvline(CUTOFF, color="0.4", ls="--", lw=0.8)
    top = max(r[1] for r in rows)
    ax.set(xlim=(3, max(16, top + 2)), ylim=(-1, len(rows)), yticks=[],
           xlabel=f"distancia CA–CA (Å):  ● en {own} ({LABELS[own]})  →  punta en {other} ({LABELS[other]})",
           title=f"Par {name_a} – {name_b}: los {len(rows)} contactos exclusivos de {own} ({LABELS[own]}) "
                 f"que {what} el criterio del trabajo\n"
                 f"(rayita negra = 1,5 × la distancia en {own}; línea punteada = límite de contacto, 8 Å)")
    ax.legend(title="residuos del par en…", fontsize=8, loc="lower right")
    path = Path(__file__).parent / f"ex8_{name_a}_{name_b}_cambio_contactos_{own}_{suffix}.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"figura: {rel(path)}")
