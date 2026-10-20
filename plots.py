"""Contact-map figures shared by exercises 7, 8 and 9."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap


def _extent(resseq):
    return (resseq[0] - 0.5, resseq[-1] + 0.5, resseq[-1] + 0.5, resseq[0] - 0.5)


def split_map(ax, lower, upper, resseq, lower_label, upper_label, title):
    """Two maps in one square: `lower` below the diagonal, `upper` above."""
    n = len(resseq)
    tri = np.tril(np.ones((n, n), bool), -1)
    img = np.where(tri, lower, upper).astype(float)   # both maps are symmetric
    ax.imshow(np.where(img, 1.0, np.nan), cmap=ListedColormap(["black"]),
              extent=_extent(resseq), interpolation="nearest")
    ax.plot([resseq[0], resseq[-1]], [resseq[0], resseq[-1]], color="0.7", lw=0.5)
    ax.set_xlabel("residuo")
    ax.set_ylabel("residuo")
    ax.set_title(f"{title}\nbajo la diagonal: {lower_label} · sobre la diagonal: {upper_label}",
                 fontsize=10)


def difference_map(ax, map_a, map_b, resseq, label_a, label_b, title, domains=None, notes=()):
    """Contacts of two structures in one map, shown once (below the diagonal;
    the map is symmetric): shared in grey, A-only in red, B-only in blue.

    domains: name -> residues, drawn as dashed boundaries and labelled on the
    top and right axes. notes: ((x, y), text, (text_x, text_y)) arrows, with
    x the smaller and y the larger residue of the pair; the empty upper
    triangle is room for their text."""
    n = len(resseq)
    code = np.full(map_a.shape, np.nan)
    code[map_a & map_b] = 0
    code[map_a & ~map_b] = 1
    code[~map_a & map_b] = 2
    code[np.triu_indices(n)] = np.nan
    colors = ["0.7", "tab:red", "tab:blue"]
    if n <= 300:
        ax.imshow(code, cmap=ListedColormap(colors), vmin=0, vmax=2, extent=_extent(resseq),
                  interpolation="nearest")
    else:
        # with many residues a cell is smaller than a pixel: draw each contact as a dot
        for value, color in enumerate(colors):
            rows, cols = np.nonzero(code == value)
            ax.scatter(np.asarray(resseq)[cols], np.asarray(resseq)[rows], s=2.5, marker="s",
                       color=color, lw=0, zorder=2 + value)
        ax.set_xlim(resseq[0] - 0.5, resseq[-1] + 0.5)
        ax.set_ylim(resseq[-1] + 0.5, resseq[0] - 0.5)
        ax.set_aspect("equal")
    ax.plot([resseq[0], resseq[-1]], [resseq[0], resseq[-1]], color="0.6", lw=0.6)
    ax.set_xlabel("residuo (el de número menor del par)")
    ax.set_ylabel("residuo (el de número mayor del par)")
    ax.set_title(title, fontsize=10)
    if domains:
        runs = []                                    # contiguous (start, end, name) pieces
        for name, members in domains.items():
            members = sorted(members)
            start = members[0]
            for a, b in zip(members, members[1:] + [None]):
                if b != a + 1:
                    runs.append((start, a, name))
                    start = b
        for start, end, name in runs:
            if start > resseq[0]:
                ax.axvline(start - 0.5, color="0.5", lw=0.6, ls="--")
                ax.axhline(start - 0.5, color="0.5", lw=0.6, ls="--")
        centers = [(start + end) / 2 for start, end, _ in runs]
        names = [name for _, _, name in runs]
        top = ax.secondary_xaxis("top")
        top.set_xticks(centers, names, fontsize=8)
        right = ax.secondary_yaxis("right")
        right.set_yticks(centers, names, fontsize=8)
    for (x, y), text, (tx, ty) in notes:
        ax.annotate(text, xy=(x, y), xytext=(tx, ty), fontsize=8, ha="center",
                    arrowprops=dict(arrowstyle="->", color="0.2", lw=0.8),
                    bbox=dict(boxstyle="round", fc="white", ec="0.6", lw=0.5))
    handles = [plt.Line2D([], [], marker="s", ls="", color=c, label=l) for c, l in
               [("0.7", "en las dos"), ("tab:red", f"solo {label_a}"), ("tab:blue", f"solo {label_b}")]]
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.9)


def kept_removed_map(ax, kept, base, resseq, title, kept_label, removed_label):
    """Half contact map of a filtered map against the one it was filtered from:
    contacts that are kept in blue, those the filter removed in light grey.
    `kept` must be a subset of `base`."""
    n = len(resseq)
    code = np.full(base.shape, np.nan)
    code[base & ~kept] = 0
    code[kept] = 1
    code[np.triu_indices(n)] = np.nan
    ax.imshow(code, cmap=ListedColormap(["0.85", "tab:blue"]), vmin=0, vmax=1,
              extent=_extent(resseq), interpolation="nearest")
    ax.plot([resseq[0], resseq[-1]], [resseq[0], resseq[-1]], color="0.6", lw=0.6)
    ax.set_xlabel("residuo (el de número menor del par)")
    ax.set_ylabel("residuo (el de número mayor del par)")
    ax.set_title(title, fontsize=10)
    handles = [plt.Line2D([], [], marker="s", ls="", color=c, label=l) for c, l in
               [("tab:blue", kept_label), ("0.85", removed_label)]]
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.9)


SS_COLORS = {"hélice": "tab:red", "hebra β": "tab:green", "otro": "0.55"}


def dihedral_kind(k, ss):
    """The CA dihedral keyed by residue k spans CA k-1..k+2; it counts as helix
    or strand only if the four residues are."""
    return next((t for t in ("hélice", "hebra β") if all(ss.get(n) == t for n in range(k - 1, k + 3))), "otro")


def secondary_structure_bands(ax, ss):
    """Shade the residues the file declares as helix (red) or strand (green)."""
    for n, kind in ss.items():
        ax.axvspan(n - 0.5, n + 0.5, color=SS_COLORS[kind], alpha=0.12, lw=0)


DIHEDRAL_CUT = -40.0   # deg: no CA dihedral of 1AKE/4AKE falls between -60 and -20


def to_window(values, cut=DIHEDRAL_CUT):
    """Map angles to [cut, cut + 360): the circle is cut at `cut` instead of at ±180."""
    return (np.asarray(values, dtype=float) - cut) % 360 + cut


def wrapped_dihedral_line(ax, dihedrals, cut=DIHEDRAL_CUT, **kwargs):
    """Plot CA dihedrals along the chain as an ordinary line, on an axis that
    cuts the circle at `cut` instead of at ±180. With the cut placed where no
    dihedral falls (between -60 and -20 deg for 1AKE/4AKE), the line never needs
    to wrap, and strands (near ±180) stay together in the middle of the plot.
    Above +180 the axis continues with -180, -160, ... up to the cut."""
    ks = sorted(dihedrals)
    y = to_window([dihedrals[k] for k in ks], cut)
    x = np.array(ks, dtype=float)
    gaps = np.flatnonzero(np.diff(x) != 1) + 1         # don't join across chain breaks
    label = kwargs.pop("label", None)
    for xs, ys in zip(np.split(x, gaps), np.split(y, gaps)):
        ax.plot(xs, ys, **kwargs)
    ax.plot([], [], label=label, **kwargs)
    ticks = np.arange(cut, cut + 361, 40)
    ax.set_yticks(ticks, [f"{(t + 180) % 360 - 180:+.0f}°" for t in ticks])
    ax.set_ylim(cut, cut + 360)
    ax.axhline(180, color="0.5", lw=0.6, ls=":")
    ax.text(1.0, 180, " ±180°", transform=ax.get_yaxis_transform(), va="center", fontsize=7, color="0.4")


def dihedral_panel(ax, dihedrals, ss, title):
    """One structure's CA dihedrals along the chain: points coloured by the
    secondary structure declared in the file, joined in order by a grey line,
    on the circular axis cut at DIHEDRAL_CUT, over helix/strand bands."""
    secondary_structure_bands(ax, ss)
    wrapped_dihedral_line(ax, dihedrals, color="0.55", lw=0.7, zorder=2)
    kinds = {k: dihedral_kind(k, ss) for k in dihedrals}
    for t in ("otro", "hebra β", "hélice"):
        ks = [k for k in dihedrals if kinds[k] == t]
        ax.scatter(ks, to_window([dihedrals[k] for k in ks]), color=SS_COLORS[t], s=12, zorder=3,
                   label=f"{t} ({len(ks)})")
    ax.set(xlabel="residuo (CA i+1)", ylabel="diedro CA(i)…CA(i+3)", title=title)
    ax.legend(fontsize=8, loc="upper left", title="color: sus 4 residuos en\nel tramo (HELIX/SHEET)", title_fontsize=8)
