"""Ejercicio 6: SASA, contactos, vecinos polares/no polares y pKa de Asp.

Proteínas: las dos de la tabla del punto v, cadena A de cada una.
  1LNI  ribonucleasa Sa (1.0 A; tiene una cadena B que se ignora)
  1LSE  lisozima de huevo (1.7 A; trae hidrógenos, 1LNI no)

La SASA se calcula sobre la cadena sola, sin aguas, ligandos ni hidrógenos:
el agua cristalográfica y la otra cadena del cristal taparían superficie que
en solución está expuesta, y los H inflarían la SASA de 1LSE respecto a 1LNI.
"""

from functools import cache

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1
from scipy.stats import pearsonr, spearmanr

from contacts import NONPOLAR, POLAR, cutoff_map, neighbor_counts, residue_points
from structures import OUT, is_amino, load

PROTEINS = {"1LNI": "A", "1LSE": "A"}
KINDS = ["CA", "CB", "CM"]
CONTACT_CUTOFF = 8.0     # A
CONTACT_MIN_SEP = 3      # |i-j| >= 3: los vecinos i+-1, i+-2 están siempre "en contacto"
RADII = [6.0, 8.0, 10.0, 12.0]
CLASSES = {"polar": POLAR, "nopolar": NONPOLAR}
ASP_REF_PKA = 4.0

# (pdb, residuo, pKa experimental) — tabla de la consigna
PKA = [("1LNI", 17, 3.72), ("1LNI", 25, 4.87), ("1LNI", 33, 2.39), ("1LNI", 79, 7.37),
       ("1LNI", 84, 3.01), ("1LNI", 93, 3.09), ("1LSE", 18, 2.66), ("1LSE", 48, 1.6),
       ("1LSE", 52, 3.68), ("1LSE", 66, 0.9), ("1LSE", 87, 2.07), ("1LSE", 101, 4.09),
       ("1LSE", 119, 3.2)]

# SASA máxima por residuo (Tien et al. 2013, valores teóricos), para la SASA relativa.
MAX_ASA = dict(A=129, R=274, N=195, D=193, C=167, Q=225, E=223, G=104, H=224, I=197,
               L=201, K=236, M=224, F=240, P=159, S=155, T=172, W=285, Y=263, V=174)

FIGS = __import__("pathlib").Path(__file__).parent


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
    print("figura: exercises/ex6_contactos_vs_sasa.png")


def iv():
    df = all_residues()
    cols = [f"{cls}_{kind}_{r:g}" for kind in KINDS for r in RADII for cls in CLASSES]
    print("promedio de vecinos por residuo, según representante y radio:")
    summary = df.groupby("pdb")[cols].mean().T.round(1)
    summary.index = pd.MultiIndex.from_tuples([(c, k, float(r)) for c, k, r in
                                               (x.split("_") for x in summary.index)],
                                              names=["clase", "repr", "radio"])
    print(summary.unstack("clase").sort_index().to_string())


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

    pivot = corr.pivot_table(index="clase", columns=["repr", "radio"], values="rho")
    fig, ax = plt.subplots(figsize=(11, 2.8))
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(pivot.shape[1]), [f"{k}\n{r:g} Å" for k, r in pivot.columns], fontsize=8)
    ax.set_yticks(range(pivot.shape[0]), pivot.index)
    for (y, x), val in np.ndenumerate(pivot.values):
        ax.text(x, y, f"{val:.2f}", ha="center", va="center", fontsize=8)
    ax.set_title(f"ρ de Spearman: nº de vecinos vs ΔpKa de Asp (n={len(rows)}, 1LNI+1LSE)")
    fig.colorbar(im, ax=ax, label="ρ")
    fig.tight_layout()
    fig.savefig(FIGS / "ex6_vecinos_vs_pka.png", dpi=150)
    print("figura: exercises/ex6_vecinos_vs_pka.png")
