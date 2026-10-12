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


def difference_map(ax, map_a, map_b, resseq, label_a, label_b, title):
    """Shared contacts in grey, A-only in red, B-only in blue."""
    code = np.full(map_a.shape, np.nan)
    code[map_a & map_b] = 0
    code[map_a & ~map_b] = 1
    code[~map_a & map_b] = 2
    cmap = ListedColormap(["0.75", "tab:red", "tab:blue"])
    ax.imshow(code, cmap=cmap, vmin=0, vmax=2, extent=_extent(resseq), interpolation="nearest")
    ax.set_xlabel("residuo")
    ax.set_ylabel("residuo")
    ax.set_title(title, fontsize=10)
    handles = [plt.Line2D([], [], marker="s", ls="", color=c, label=l) for c, l in
               [("0.75", "ambas"), ("tab:red", f"solo {label_a}"), ("tab:blue", f"solo {label_b}")]]
    ax.legend(handles=handles, loc="upper right", fontsize=7, framealpha=0.9)
