"""Ejercicio 7: adenilato quinasa de E. coli, cerrada (1AKE, con Ap5A) vs
abierta (4AKE), solo carbonos alfa, cadena A de cada una."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from geometry import ca_dihedrals
from pairs import compare
from plots import SS_COLORS, dihedral_kind, dihedral_panel
from structures import ca_atoms, load, secondary_structure

CLOSED, OPEN = "1AKE", "4AKE"
LABELS = {CLOSED: "cerrada", OPEN: "abierta"}
# Dominios de la adenilato quinasa según la bibliografía (no figuran en el archivo).
DOMAINS = {"NMP": range(30, 60), "LID": range(122, 160)}
DOMAINS["CORE"] = [n for n in range(1, 215) if not any(n in r for r in DOMAINS.values())]
# Groups of contacts that only the closed form has (criterion of the paper),
# pointed at in the contact map: ((smaller residue, larger residue), label, text position).
NOTES = [((55, 163), "NMP (53–58) con 157–170:\nse tocan solo en la cerrada", (130, 60)),
         ((15, 132), "CORE (14–18) con LID (131–133):\nsolo en la cerrada", (60, 20)),
         ((137, 202), "LID (137–138) con CORE (201–202)", (175, 120))]


def run():
    result = compare(CLOSED, OPEN, figs_dir=Path(__file__).parent, prefix="ex7_", domains=DOMAINS, notes=NOTES, labels=LABELS)
    contact_changes()
    return result


def dihedral_explainer():
    """Figure that explains how to read the CA dihedral plot, on 1AKE data."""
    cas = ca_atoms(load(CLOSED)[0]["A"])
    dihedrals = ca_dihedrals(cas)
    ss = secondary_structure(CLOSED)
    kind = {k: dihedral_kind(k, ss) for k in dihedrals}
    colors = SS_COLORS

    fig = plt.figure(figsize=(16, 11), layout="constrained")
    grid = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.4])

    # 1. what a dihedral is: looking down the CA2-CA3 bond
    ax = fig.add_subplot(grid[0, 0])
    ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(plt.Circle((0, 0), 0.12, color="k"))
    ax.text(0.17, -0.05, "tubo CA2→CA3\n(visto de frente)", fontsize=8, va="top")
    ax.annotate("", xy=(0, 1), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", lw=2.5, color="k"))
    ax.text(0, 1.1, "CA 1", ha="center", fontsize=10)
    for angle, label, color in ((49, "CA 4 en una\nhélice: +49°", "tab:red"),
                                (-153, "CA 4 en una\nhebra β: −153°", "tab:green")):
        rad = np.radians(90 - angle)    # measured clockwise from CA1, as seen down the bond
        ax.annotate("", xy=(np.cos(rad), np.sin(rad)), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", lw=2, color=color, ls="--"))
        ax.text(1.25 * np.cos(rad), 1.25 * np.sin(rad), label, ha="center", va="center", color=color, fontsize=9)
    ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.4)
    ax.set_title("1. Mirando a lo largo de la unión CA2→CA3,\nel diedro es el giro de CA4 respecto de CA1", fontsize=10)

    # 2. the dihedral lives on a circle
    ax = fig.add_subplot(grid[0, 1], projection="polar")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.set_thetalim(0, 2 * np.pi)
    for t in ("otro", "hebra β", "hélice"):
        ks = [k for k in dihedrals if kind[k] == t]
        ax.scatter(np.radians([dihedrals[k] % 360 for k in ks]), [1.0 if t != "otro" else 0.75] * len(ks),
                   color=colors[t], s=14, alpha=0.8)
    ax.set_ylim(0, 1.2); ax.set_yticks([])
    ax.set_xticks(np.radians([0, 49, 90, 180, 270]), ["0°", "+49°", "+90°", "+180° = −180°", "−90°"])
    ax.set_title("2. El diedro vive en un círculo:\n+170° y −170° quedan uno al lado del otro", fontsize=10)

    # 3. the same values unrolled, as in the exercise figure
    ax = fig.add_subplot(grid[0, 2])
    ax.axhspan(140, 180, color="tab:green", alpha=0.1)
    ax.axhspan(-180, -140, color="tab:green", alpha=0.1)
    for t in ("otro", "hebra β", "hélice"):
        ks = [k for k in dihedrals if kind[k] == t]
        ax.scatter(ks, [dihedrals[k] for k in ks], color=colors[t], s=12, label=f"{t} ({len(ks)})")
    ax.annotate("las hebras quedan arriba Y abajo:\nson el mismo grupo, partido",
                xy=(104, -160), xytext=(115, -60), fontsize=8, color="tab:green",
                arrowprops=dict(arrowstyle="->", color="tab:green"))
    ax.annotate("", xy=(195, 172), xytext=(150, -45), arrowprops=dict(arrowstyle="->", color="tab:green"))
    ax.set(ylim=(-180, 180), yticks=[-180, -90, 0, 49, 90, 180], xlabel="residuo", ylabel="diedro CA (°)")
    ax.set_title("3. El mismo círculo \"desenrollado\" sobre una recta (el gráfico del ejercicio)", fontsize=10)
    ax.legend(fontsize=8, loc="upper left", title="según HELIX/SHEET del archivo", title_fontsize=8)
    # 4. the same values with the circle cut where there are no dihedrals
    dihedral_panel(fig.add_subplot(grid[1, :]), dihedrals, ss,
                   "4. El collar cortado en −40°, donde no hay ningún diedro: las hebras quedan enteras alrededor "
                   "de ±180°, y la línea une los residuos en orden, como en el gráfico de ángulos")
    fig.savefig(Path(__file__).parent / "ex7_diedros_explicacion.png", dpi=150)
    plt.close(fig)
    path = Path(__file__).parent / "ex7_diedros_explicacion.png"
    print("figura: exercises/ex7_diedros_explicacion.png")


def contact_changes():
    """One arrow per contact exclusive to the closed form: from its CA-CA
    distance in 1AKE (dot) to its distance in 4AKE (arrow head), sorted by how
    much it changes, with a tick at 1.5 x the closed distance (the paper's
    threshold) and coloured by the domains of the pair."""
    import json
    from pairs import PAPER_MIN_RATIO, CUTOFF
    result = json.loads((Path(__file__).parents[1] / "data" / "out" / f"ex7_{CLOSED}_{OPEN}.json").read_text())
    ca = {p: {a.get_parent().id[1]: a for a in ca_atoms(load(p)[0]["A"])} for p in (CLOSED, OPEN)}
    paper = {tuple(x) for x in result["solo_segun_trabajo"][CLOSED]}
    dom = lambda n: next(d for d, rng in DOMAINS.items() if n in rng)
    rows = []
    for i, j in map(tuple, result["exclusivos"][CLOSED]):
        rows.append((ca[CLOSED][i] - ca[CLOSED][j], ca[OPEN][i] - ca[OPEN][j], i, j))
    rows.sort(key=lambda r: r[1] / r[0])                      # least change at the bottom
    palette = {"CORE": "0.55", "NMP": "tab:olive", "LID": "tab:cyan", "CORE – LID": "tab:orange",
               "CORE – NMP": "tab:green", "LID – NMP": "tab:purple"}

    fig, ax = plt.subplots(figsize=(10, 18), layout="constrained")
    seen = set()
    for y, (dc, do, i, j) in enumerate(rows):
        key = " – ".join(sorted({dom(i), dom(j)}))
        color = palette[key]
        strong = (i, j) in paper
        ax.annotate("", xy=(do, y), xytext=(dc, y),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=1.8 if strong else 1, alpha=1 if strong else 0.6,
                                    mutation_scale=9, shrinkA=0, shrinkB=0))
        ax.scatter(dc, y, color=color, s=14, zorder=3, label=None if key in seen else key)
        seen.add(key)
        ax.plot([PAPER_MIN_RATIO * dc] * 2, [y - 0.35, y + 0.35], color="k", lw=1)   # 1.5 x threshold
        ax.text(31.5, y, f"{i}–{j}", va="center", fontsize=7, color="k" if strong else "0.5",
                fontweight="bold" if strong else "normal")
    ax.axvline(CUTOFF, color="0.4", ls="--", lw=0.8)
    ax.text(CUTOFF, len(rows) + 0.3, " límite de contacto (8 Å)", fontsize=8, color="0.3")
    n_strong = sum((i, j) in paper for _, _, i, j in rows)
    ax.axhline(len(rows) - n_strong - 0.5, color="k", lw=0.8)
    ax.text(30, len(rows) - n_strong - 0.2, f"↑ {n_strong} cumplen la regla del trabajo\n"
            f"↓ {len(rows) - n_strong} no la cumplen", ha="right", va="bottom", fontsize=9)
    ax.set(xlim=(4, 34), ylim=(-1, len(rows) + 1), yticks=[],
           xlabel=f"distancia CA–CA (Å):  ● en {CLOSED} ({LABELS[CLOSED]})  →  punta en {OPEN} ({LABELS[OPEN]})",
           title=f"Cómo cambia cada contacto exclusivo de {CLOSED} ({LABELS[CLOSED]}) al abrirse la enzima\n"
                 f"(rayita negra = 1,5 × la distancia en la cerrada: la flecha tiene que pasarla para cumplir la regla)")
    ax.legend(title="dominios del par", fontsize=8, loc="center right", bbox_to_anchor=(0.9, 0.3))
    path = Path(__file__).parent / "ex7_cambio_contactos.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"figura: exercises/{path.name}")
