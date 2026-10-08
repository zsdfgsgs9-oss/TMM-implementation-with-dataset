"""Plot the retained space-charge map against the round-beam envelope check."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SOURCE = (
    ROOT
    / "RSI_revision_simulation_package_20260901"
    / "06_source_code"
    / "revision_simulation"
    / "results"
    / "space_charge_kv_comparison.csv"
)
OUT = ROOT / "figures" / "fig_space_charge_reviewer4"

with SOURCE.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))

assert len(rows) == 9 and float(rows[0]["I_A"]) == 0
current_na = np.array([float(row["I_A"]) * 1e9 for row in rows[1:]])
kv_final = np.array([float(row["KV_diameter_final_nm"]) for row in rows[1:]])
tmm_final = np.array([2 * float(row["TMM_R_total_final_nm"]) for row in rows[1:]])
kv_peak = np.array([float(row["KV_peak_delta_diameter_nm"]) for row in rows[1:]])
tmm_peak = np.array([float(row["native_peak_delta_diameter_nm"]) for row in rows[1:]])
r0_m = float(rows[0]["R0_nm"]) * 1e-9
nominal_source_radius = 18e-9
source_area = np.pi * nominal_source_radius**2

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "font.size": 8.6,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.bf": "Times New Roman:bold",
        "mathtext.fallback": None,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)
fig, axes = plt.subplots(1, 2, figsize=(7.05, 3.15), constrained_layout=True)
colors = {"map": "#087e8b", "kv": "#c44e52"}

axes[0].plot(current_na, tmm_final, "o-", color=colors["map"], lw=1.6,
             markersize=4.3, label="Envelope map")
axes[0].plot(current_na, kv_final, "s--", color=colors["kv"], lw=1.4,
             markersize=4, label="K-V-type envelope")
axes[0].set_ylabel("Image-plane diameter (nm)")
axes[0].set_ylim(31.6, 37.1)
axes[0].legend(frameon=False, loc="lower left", fontsize=9)

axes[1].plot(current_na, tmm_peak, "o-", color=colors["map"], lw=1.6,
             markersize=4.3, label="Envelope map")
axes[1].plot(current_na, kv_peak, "s--", color=colors["kv"], lw=1.4,
             markersize=4, label="K-V-type envelope")
axes[1].set_yscale("log")
axes[1].set_ylabel("Maximum diameter increment (nm)")
axes[1].set_ylim(8e-4, 3e2)

ticks = [1e-3, 1e-2, 1e-1, 1, 10, 100]
tick_labels = ["0.001", "0.01", "0.1", "1", "10", "100"]
for ax, label in zip(axes, ["(a)", "(b)"]):
    ax.set_xscale("log")
    ax.set_xlim(8e-4, 1.3e2)
    ax.set_xticks(ticks, tick_labels)
    ax.set_xlabel("Beam current at emission, $I$ (nA)")
    ax.grid(True, which="major", color="#d7dde0", lw=0.55)
    ax.tick_params(direction="in", top=False, right=True, length=3)
    ax.text(0.025, 0.965, label, transform=ax.transAxes,
            ha="left", va="top", fontweight="bold")

def current_to_density(current):
    return np.asarray(current) * 1e-9 / source_area

def density_to_current(density):
    return np.asarray(density) * source_area * 1e9

top = axes[0].secondary_xaxis(
    "top", functions=(current_to_density, density_to_current)
)
top.set_xscale("log")
top.set_xticks([1e3, 1e5, 1e6, 1e7, 1e8])
top.set_xticklabels(["$10^3$", "$10^5$", "$10^6$", "$10^7$", "$10^8$"])
top.set_xlabel("Emission-plane $J_0$ (A/m$^2$)", labelpad=4)
top.tick_params(direction="in", length=3)

fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(OUT.with_suffix(".png"), dpi=250, bbox_inches="tight")
print(f"Saved {OUT.with_suffix('.pdf')}")
print(f"Sampled R0 = {r0_m * 1e9:.6f} nm")
print(f"Nominal 18 nm source: 1 nA = {current_to_density(1):.0f} A/m^2")
