"""Draw the equivalent round-beam cross-section used by the map."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch


OUT = Path(__file__).resolve().parent / "figures" / "fig_space_charge_envelope"
plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "font.size": 11,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.bf": "Times New Roman:bold",
        "mathtext.fallback": None,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

fig, ax = plt.subplots(figsize=(4.0, 2.5))
cx = -0.20
force_color = "#233c48"
ax.add_patch(Circle((cx, 0), 0.85, facecolor="#e1f1f0", edgecolor="#087e8b", lw=1.7))
ax.add_patch(Circle((cx, 0), 0.65, fill=False, edgecolor="#c44e52", lw=1.5, ls="--"))
ax.plot(cx, 0, "k+", ms=8, mew=1.2)
for radius, angle in [(0.30, 18), (0.38, 115), (0.46, 210), (0.54, 300)]:
    direction = np.array([np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle))])
    start = np.array([cx, 0]) + radius * direction
    end = start + 0.50 * radius * direction
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>",
                                 mutation_scale=7, lw=1.25, color=force_color))
    ax.plot(*start, "o", color=force_color, ms=2.9)

ax.plot([0.83, 1.09], [0.55, 0.55], color="#087e8b", lw=1.7)
ax.text(1.20, 0.55, r"$R(I)$", va="center")
ax.plot([0.83, 1.09], [0.18, 0.18], color="#c44e52", lw=1.5, ls="--")
ax.text(1.20, 0.18, r"$R_{\mathrm{id}}$", va="center")
ax.add_patch(FancyArrowPatch((0.83, -0.21), (1.09, -0.21), arrowstyle="-|>",
                             mutation_scale=9, lw=1.25, color=force_color))
ax.text(1.20, -0.21, r"$\mathbf{F}_{\mathrm{sc}}$", va="center")
ax.plot(0.96, -0.59, "k+", ms=8, mew=1.2)
ax.text(1.20, -0.59, "Centroid", va="center")
ax.set_xlim(-1.13, 2.05)
ax.set_ylim(-0.99, 0.99)
ax.set_aspect("equal")
ax.axis("off")
fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(OUT.with_suffix(".png"), dpi=240, bbox_inches="tight")
print(f"Saved {OUT.with_suffix('.pdf')}")
