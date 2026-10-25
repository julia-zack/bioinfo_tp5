"""Ejercicio 6: SASA, contactos, vecinos polares/no polares y pKa de Asp.

Proteínas: las dos de la tabla del punto v, cadena A de cada una.
  1LNI  ribonucleasa Sa (1.0 A; tiene una cadena B que se ignora)
  1LSE  lisozima de huevo (1.7 A; trae hidrógenos, 1LNI no)

La SASA se calcula sobre la cadena sola, sin aguas, ligandos ni hidrógenos:
el agua cristalográfica y la otra cadena del cristal taparían superficie que
en solución está expuesta, y 1LSE trae H que 1LNI no tiene, así que se sacan
para medir las dos igual (medido: dejarlos sube la SASA de 1LSE solo un 1%).
"""

from functools import cache

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1
from scipy.stats import pearsonr, spearmanr

from contacts import NONPOLAR, POLAR, cutoff_map, neighbor_counts, residue_points
from structures import OUT, is_amino, load, rel

PROTEINS = {"1LNI": "A", "1LSE": "A"}
KINDS = ["CA", "CB", "CM"]
CONTACT_CUTOFF = 8.0     # A
CONTACT_MIN_SEP = 3      # |i-j| >= 3: los vecinos i+-1, i+-2 están siempre "en contacto"
RADII = [6.0, 8.0, 10.0, 12.0]
CLASSES = {"polar": POLAR, "nopolar": NONPOLAR}
ASP_REF_PKA = 4.0
PLOT_RADIUS = 8.0        # A, radio de los gráficos de 6v (el de 6ii y del notebook)
PLOT_KIND = "CB"         # punto de los gráficos con rectas de 6v
# Misma convención en todo el TP (y en los notebooks): color = grupo de vecinos,
# forma = proteína.
GROUP_COLORS = {"polar": "tab:blue", "nopolar": "tab:orange"}
PROTEIN_MARKERS = {"1LNI": "o", "1LSE": "s"}

# (pdb, residuo, pKa experimental) — tabla de la consigna
PKA = [("1LNI", 17, 3.72), ("1LNI", 25, 4.87), ("1LNI", 33, 2.39), ("1LNI", 79, 7.37),
       ("1LNI", 84, 3.01), ("1LNI", 93, 3.09), ("1LSE", 18, 2.66), ("1LSE", 48, 1.6),
       ("1LSE", 52, 3.68), ("1LSE", 66, 0.9), ("1LSE", 87, 2.07), ("1LSE", 101, 4.09),
       ("1LSE", 119, 3.2)]

# SASA máxima por residuo (Tien et al. 2013, valores teóricos), para la SASA relativa.
MAX_ASA = dict(A=129, R=274, N=195, D=193, C=167, Q=225, E=223, G=104, H=224, I=197,
               L=201, K=236, M=224, F=240, P=159, S=155, T=172, W=285, Y=263, V=174)

FIGS = __import__("pathlib").Path(__file__).parent


# Aspárticos cuya ausencia cambia mucho alguna de las rectas de 6v: (pdb, residuo, grupo afectado)
INFLUENTIAL = [("1LNI", 79, "nopolar"), ("1LSE", 66, "polar"), ("1LSE", 52, "polar")]


def trend_lines(rows, kind=PLOT_KIND, radius=PLOT_RADIUS):
    """Polar and nonpolar neighbours vs pKa on the same axes, each with its
    least-squares line (the "X"): with all 13 Asp, and once without each of the
    residues that move a line the most, to show how fragile the trends are.
    In the "without" panels the lines for all 13 are drawn dotted for reference
    and the removed residue is marked with a cross."""
    groups = (("polar", "polares", GROUP_COLORS["polar"]), ("nopolar", "no polares", GROUP_COLORS["nopolar"]))
    names = {"polar": "polares", "nopolar": "no polares"}
    col = lambda cls: f"{cls}_{kind}_{radius:g}"

    def fit(data, cls):
        return np.polyfit(data[col(cls)], data.pka_exp, 1)

    full = {cls: fit(rows, cls) for cls, _, _ in groups}
    panels = [(rows, f"los {len(rows)} Asp", None, None)]
    for pdb_id, resseq, cls in INFLUENTIAL:
        drop = (rows.pdb == pdb_id) & (rows.resseq == resseq)
        panels.append((rows[~drop], f"sin el Asp{resseq} de {pdb_id} (mirar la recta de los {names[cls]})",
                       rows[drop], cls))

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharey=True, sharex=True)
    print(f"\nrectas de tendencia ({kind}, {radius:g} Å): pendiente = cambio de pKa por vecino")
    for ax, (data, title, removed, changed) in zip(axes.flat, panels):
        for cls, label, color in groups:
            x, y = data[col(cls)], data.pka_exp
            slope, intercept = fit(data, cls)
            for pdb_id, marker in PROTEIN_MARKERS.items():
                sel = data.pdb == pdb_id
                ax.scatter(x[sel], y[sel], color=color, marker=marker, s=30, alpha=0.7,
                           label=f"vecinos {label} ({pdb_id})")
            for xi, yi, n in zip(x, y, data.resseq):
                ax.annotate(str(n), (xi, yi), textcoords="offset points", xytext=(4, 3),
                            fontsize=7, color=color)
            xs = np.array([0, rows[col(cls)].max()])
            ax.plot(xs, intercept + slope * xs, color=color, lw=2)
            if removed is not None:
                if cls == changed:          # reference line only for the group whose line changes
                    full_slope, full_intercept = full[cls]
                    ax.plot(xs, full_intercept + full_slope * xs, color=color, lw=1, ls=":")
                ax.scatter(removed[col(cls)], removed.pka_exp, marker="x", s=80, color=color, lw=2)
            print(f"  {title.split(' (')[0]:22} {label:11} pendiente {slope:+.2f}")
        ax.axhline(ASP_REF_PKA, color="0.6", ls="--", lw=0.8)
        ax.set_title(title, fontsize=10)
    for ax in axes[1]:
        ax.set_xlabel(f"nº de vecinos (desde {kind}, radio {radius:g} Å)")
    for ax in axes[:, 0]:
        ax.set_ylabel("pKa medido (línea gris: referencia 4)")
    axes[0][0].legend(fontsize=7)
    fig.suptitle("Tendencia de polares y no polares (cada Asp aparece dos veces: una por grupo; círculo 1LNI, "
                 "cuadrado 1LSE).\nEn los paneles \"sin\": punteado = cómo era con los 13 la recta que cambia; cruz = el Asp que se sacó",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGS / "ex6_vecinos_vs_pka_rectas.png", dpi=150)
    plt.close(fig)
    print(f"figura: {rel(FIGS / 'ex6_vecinos_vs_pka_rectas.png')}")


def bare_chain(pdb_id, chain_id):
    """Copy of the chain with only standard residues and no hydrogens."""
    chain = load(pdb_id)[0][chain_id].copy()
    for residue in list(chain):
        if not is_amino(residue):
            chain.detach_child(residue.id)
            continue
        for atom in list(residue):
            if atom.element == "H":
                residue.detach_child(atom.id)
    return chain


def residue_table(pdb_id, chain_id):
    """One row per residue: SASA, RSA, contact counts and neighbour counts."""
    chain = bare_chain(pdb_id, chain_id)
    ShrakeRupley().compute(chain, level="R")
    residues, _ = residue_points(chain, "CA")
    df = pd.DataFrame({
        "pdb": pdb_id,
        "resseq": [r.id[1] for r in residues],
        "aa": [seq1(r.get_resname()) for r in residues],
        "sasa": [r.sasa for r in residues],
    })
    df["rsa"] = df.sasa / df.aa.map(MAX_ASA)
    resseq = df.resseq.to_numpy()
    for kind in KINDS:
        _, coords = residue_points(chain, kind)
        df[f"contactos_{kind}"] = cutoff_map(coords, resseq, CONTACT_CUTOFF, CONTACT_MIN_SEP).sum(1)
        for radius in RADII:
            for cls, n in neighbor_counts(coords, df.aa, radius, CLASSES).items():
                df[f"{cls}_{kind}_{radius:g}"] = n
    return df


@cache
def _all_residues():
    return pd.concat([residue_table(p, c) for p, c in PROTEINS.items()], ignore_index=True)


def all_residues():
    """Cached: SASA is the slow part and every sub-item needs the same table."""
    return _all_residues().copy()


def i():
    df = all_residues()
    for pdb_id, sub in df.groupby("pdb"):
        print(f"== {pdb_id}: {len(sub)} residuos, SASA total {sub.sasa.sum():.0f} A^2")
        print(sub.nlargest(5, "sasa")[["resseq", "aa", "sasa", "rsa"]].to_string(index=False))
        print("  más enterrados:")
        print(sub.nsmallest(5, "sasa")[["resseq", "aa", "sasa", "rsa"]].to_string(index=False))
    path = OUT / "ex6_residuos.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"\ntabla completa por residuo en data/out/{path.name}")


def ii():
    df = all_residues()
    print(f"contacto = distancia <= {CONTACT_CUTOFF} A y |i-j| >= {CONTACT_MIN_SEP}")
    print(df.groupby("pdb")[[f"contactos_{k}" for k in KINDS]].agg(["mean", "max"]).round(1))


def iii():
    df = all_residues()
    rows = []
    for pdb_id, sub in df.groupby("pdb"):
        for kind in KINDS:
            for target in ("sasa", "rsa"):
                x, y = sub[f"contactos_{kind}"], sub[target]
                rows.append({"pdb": pdb_id, "contactos": kind, "vs": target,
                             "spearman": spearmanr(x, y).statistic,
                             "pearson": pearsonr(x, y).statistic})
    print(pd.DataFrame(rows).round(2).to_string(index=False))

    fig, axes = plt.subplots(len(PROTEINS), len(KINDS), figsize=(12, 7), sharey=True)
    for row, (pdb_id, sub) in zip(axes, df.groupby("pdb")):
        for ax, kind in zip(row, KINDS):
            ax.scatter(sub[f"contactos_{kind}"], sub.rsa, s=12, alpha=0.6)
            rho = spearmanr(sub[f"contactos_{kind}"], sub.rsa).statistic
            ax.set_title(f"{pdb_id} — contactos por {kind} (ρ={rho:.2f})", fontsize=10)
            ax.set_xlabel(f"nº de contactos (≤{CONTACT_CUTOFF:g} Å, |i−j|≥{CONTACT_MIN_SEP})")
        row[0].set_ylabel("SASA relativa")
    fig.tight_layout()
    fig.savefig(FIGS / "ex6_contactos_vs_sasa.png", dpi=150)
    print(f"figura: {rel(FIGS / 'ex6_contactos_vs_sasa.png')}")


def iv():
    df = all_residues()
    cols = [f"{cls}_{kind}_{r:g}" for kind in KINDS for r in RADII for cls in CLASSES]
    print("promedio de vecinos por residuo, según representante y radio:")
    summary = df.groupby("pdb")[cols].mean().T.round(1)
    summary.index = pd.MultiIndex.from_tuples([(c, k, float(r)) for c, k, r in
                                               (x.split("_") for x in summary.index)],
                                              names=["clase", "repr", "radio"])
    print(summary.unstack("clase").sort_index().to_string())

    # Cómo crecen los vecinos con el radio. Un panel por punto de medición: la
    # línea es el promedio de los aminoácidos y la banda, la mitad central
    # (percentiles 25 a 75), para ver cuánto varían entre sí. El último panel
    # compara el crecimiento con el del volumen de la esfera.
    kind_styles = {"CA": ("-", "o"), "CB": ("--", "s"), "CM": (":", "^")}
    radii = np.array(RADII)
    fig, axes = plt.subplots(len(PROTEINS), len(KINDS) + 1, figsize=(16, 8), sharex=True)
    for (pdb_id, sub), row in zip(df.groupby("pdb"), axes):
        ax_rel = row[-1]
        for ax, (kind, (ls, marker)) in zip(row, kind_styles.items()):
            for cls, label in (("polar", "polares"), ("nopolar", "no polares")):
                values = [sub[f"{cls}_{kind}_{r:g}"] for r in RADII]
                ax.plot(radii, [v.mean() for v in values], ls=ls, marker=marker,
                        color=GROUP_COLORS[cls], label=f"{label} (promedio)")
                ax.fill_between(radii, [v.quantile(0.25) for v in values], [v.quantile(0.75) for v in values],
                                color=GROUP_COLORS[cls], alpha=0.2, label=f"{label} (mitad central)")
            ax.set_title(f"{pdb_id}, desde {kind}: vecinos por aminoácido ({len(sub)} aa)", fontsize=9)
            total = np.array([(sub[f"polar_{kind}_{r:g}"] + sub[f"nopolar_{kind}_{r:g}"]).mean()
                              for r in RADII])
            ax_rel.plot(radii, total / total[0], ls=ls, marker=marker, color="0.25", label=f"desde {kind}")
        ax_rel.plot(radii, (radii / radii[0]) ** 3, color="tab:red", ls="--", lw=1,
                    label="si crecieran con el volumen, (r/6)³")
        ax_rel.set_title(f"{pdb_id}: crecimiento del promedio (relativo a 6 Å)", fontsize=9)
        ax_rel.set_ylabel("veces lo que hay a 6 Å")
        row[0].set_ylabel("nº de vecinos")
        for ax in row[1:len(KINDS)]:
            ax.sharey(row[0])
    for ax in axes[-1]:
        ax.set_xlabel("radio de la esfera (Å)")
        ax.set_xticks(radii)
    axes[0][0].legend(fontsize=7)
    axes[0][-1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIGS / "ex6_vecinos_vs_radio.png", dpi=150)
    plt.close(fig)
    print(f"figura: {rel(FIGS / 'ex6_vecinos_vs_radio.png')}")


def pka_table():
    df = all_residues().set_index(["pdb", "resseq"])
    rows = df.loc[[(p, n) for p, n, _ in PKA]].reset_index()
    assert (rows.aa == "D").all()
    rows["pka_exp"] = [x for _, _, x in PKA]
    rows["dpka"] = rows.pka_exp - ASP_REF_PKA
    return rows


def v():
    rows = pka_table()
    print(rows[["pdb", "resseq", "pka_exp", "dpka", "rsa", "polar_CB_8", "nopolar_CB_8",
                "polar_CM_8", "nopolar_CM_8"]].round(2).to_string(index=False))

    corr = []
    for kind in KINDS:
        for r in RADII:
            for cls in CLASSES:
                res = spearmanr(rows[f"{cls}_{kind}_{r:g}"], rows.dpka)
                corr.append({"repr": kind, "radio": r, "clase": cls,
                             "rho": res.statistic, "p": res.pvalue})
    corr = pd.DataFrame(corr)
    print(f"\nSpearman entre nº de vecinos y ΔpKa (n={len(rows)}):")
    print(corr.pivot_table(index=["repr", "radio"], columns="clase", values=["rho", "p"])
          .round(2).to_string())
    rsa = spearmanr(rows.rsa, rows.dpka)
    print(f"\nreferencia: Spearman(SASA relativa, ΔpKa) = {rsa.statistic:.2f} (p={rsa.pvalue:.2f})")

    trend_lines(rows)

    # Gráfico de los datos (no de las correlaciones): un punto por aspártico.
    # Un solo radio, fijado de antemano (el mismo de 6ii y del notebook), para
    # no elegir el que mejor correlaciona.
    radius = PLOT_RADIUS
    fig, axes = plt.subplots(len(CLASSES), len(KINDS), figsize=(12, 7), sharey=True)
    for row, (cls, label) in zip(axes, (("polar", "polares"), ("nopolar", "no polares"))):
        for ax, kind in zip(row, KINDS):
            for pdb_id, sub in rows.groupby("pdb"):
                x = sub[f"{cls}_{kind}_{radius:g}"]
                ax.scatter(x, sub.pka_exp, color=GROUP_COLORS[cls], marker=PROTEIN_MARKERS[pdb_id],
                           s=30, label=pdb_id)
                for xi, yi, n in zip(x, sub.pka_exp, sub.resseq):
                    ax.annotate(str(n), (xi, yi), textcoords="offset points", xytext=(4, 3), fontsize=7)
            ax.axhline(ASP_REF_PKA, color="0.6", ls="--", lw=0.8)
            ax.set_title(f"vecinos {label}, medidos desde {kind}", fontsize=10)
            ax.set_xlabel(f"nº de vecinos {label} (radio {radius:g} Å)")
        row[0].set_ylabel("pKa medido (línea: referencia 4)")
    axes[0][0].legend(fontsize=8, loc="upper left")   # círculo 1LNI, cuadrado 1LSE
    fig.suptitle(f"Asp de la tabla (n={len(rows)}): pKa medido vs número de vecinos", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "ex6_vecinos_vs_pka.png", dpi=150)
    plt.close(fig)
    print(f"figura: {rel(FIGS / 'ex6_vecinos_vs_pka.png')}")
