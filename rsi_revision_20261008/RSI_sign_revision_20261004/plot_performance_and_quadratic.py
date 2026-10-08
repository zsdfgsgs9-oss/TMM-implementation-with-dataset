"""Redraw measured performance and COMSOL-only quadratic-field sensitivity."""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "RSI最终版报告与汇总数据_精简版_20261004"
FIGURES = ROOT / "figures"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--performance-only", action="store_true")
args = parser.parse_args()


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


font = font_manager.findfont("Times New Roman", fallback_to_default=False)
print("Font:", font)
plt.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def save(fig, name):
    for extension in ("pdf", "png"):
        fig.savefig(FIGURES / f"{name}.{extension}", dpi=300,
                    facecolor="white", bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)


timing = read_csv(SOURCE / "06_五千万粒子实测" /
                  "tmm_timing_summary_with_50million.csv")
n = np.array([int(row["N"]) for row in timing])
core = np.array([float(row["core_median_s"]) for row in timing])
minimum = np.array([float(row["core_min_s"]) for row in timing])
maximum = np.array([float(row["core_max_s"]) for row in timing])
assert np.array_equal(n, [1000, 5000, 10000, 50000, 100000,
                          1000000, 10000000, 50000000])
assert all(int(row["formal_repeats"]) == 5 and int(row["warmups"]) == 2
           for row in timing)
matched_time = float(core[n == 50000][0])
solve_time = 1210.97878
ratio = solve_time / matched_time
print(f"Matched solver/core ratio: {ratio:.4f}")

# These five points are retained author measurements, not matched new solves.
historical_n = np.array([1000, 5000, 10000, 50000, 100000])
historical_time = np.array([182, 227, 290, 613, 1048])
fig, axes = plt.subplots(1, 2, figsize=(7.15, 3.65),
                         gridspec_kw={"width_ratios": [1.55, 1]})
ax, comparison = axes
ax.set_title("(a) Runtime scaling", loc="left", pad=10)
ax.errorbar(n, core, yerr=np.vstack((core-minimum, maximum-core)),
            fmt="o-", color="#187a67", markersize=4, capsize=2,
            linewidth=1.3, label="TMM core: median (min--max)")
ax.plot(historical_n, historical_time, "s--", color="#777777",
        markersize=4, linewidth=1.1, label="COMSOL: earlier configuration")
ax.plot(50000, solve_time, "*", color="#b7303f", markersize=9,
        label="COMSOL: Case 06 input")
ax.set(xscale="log", yscale="log", xlabel="Particle count, n",
       ylabel="Transport time (s)", xlim=(700, 80000000), ylim=(0.002, 100000))
ax.set_xticks([1000, 10000, 100000, 1000000, 10000000, 50000000],
              labels=[r"$10^3$", r"$10^4$", r"$10^5$", r"$10^6$",
                      r"$10^7$", r"$5\times10^7$"])
ax.grid(which="major", color="0.88", linewidth=0.5)
ax.text(0.02, 0.98, "COMSOL memory limit above ~100k particles\n"
        "(earlier configuration)",
        transform=ax.transAxes, va="top", fontsize=8.4)
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.09, 0.005),
           frameon=False, borderaxespad=0, handlelength=2.3)
comparison.set_title("(b) Same input: 50,000 particles", loc="left", pad=10)
comparison.bar([0, 1], [matched_time, solve_time], width=0.52,
               color=["#187a67", "#b7303f"])
comparison.errorbar(0, matched_time,
                    yerr=[[matched_time-minimum[n == 50000][0]],
                          [maximum[n == 50000][0]-matched_time]],
                    fmt="none", ecolor="black", capsize=3, linewidth=1)
comparison.set(yscale="log", ylabel="Transport time (s)", ylim=(0.03, 12000),
               xticks=[0, 1], xticklabels=["TMM\ntransport core", "COMSOL\nparticle solve"])
comparison.text(0, matched_time*1.8, "0.145 s", ha="center", fontsize=9)
comparison.text(1, solve_time*1.8, "1211 s", ha="center", fontsize=9)
comparison.grid(axis="y", which="major", color="0.88", linewidth=0.5)
comparison.set_axisbelow(True)
comparison.text(0.5, -0.43, rf"Solver/core time ratio: ${ratio/1000:.2f}\times10^3$",
                transform=comparison.transAxes, ha="center", va="top", fontsize=9)
fig.subplots_adjust(left=0.09, right=0.985, top=0.86, bottom=0.36, wspace=0.37)
save(fig, "fig_performance_final_20261004")
if args.performance_only:
    print("Performance PDF and PNG updated; quadratic figure unchanged.")
    raise SystemExit(0)

rows = read_csv(SOURCE / "02_双场口径rtol1e7" /
                "dual_origin_truncation_summary.csv")
intervals = read_csv(SOURCE / "02_双场口径rtol1e7" / "paired_bootstrap_CI95.csv")
ci_by_case = {(row["field_reference"], row["case_name"]): row for row in intervals}
for origin in ("entrance", "midpoint"):
    selected = [row for row in rows if row["field_reference"] == origin]
    baseline = float(selected[0]["NR_full_FW99_nm"])
    for row in selected:
        delta = float(row["NR_full_FW99_nm"]) - baseline
        assert abs(delta-float(row["H_minus_H0_FW99_nm"])) < 1e-12

fig, axes = plt.subplots(1, 2, figsize=(7.15, 2.75), sharey=True)
for ax, origin, title, color in zip(
        axes, ("midpoint", "entrance"),
        ("(a) Midpoint-referenced linear baseline", "(b) Entrance-referenced baseline"),
        ("#187a67", "#326cac")):
    selected = [row for row in rows if row["field_reference"] == origin]
    h = np.array([float(row["H_Tpm2"]) for row in selected])
    delta = np.array([float(row["H_minus_H0_FW99_nm"])*1000 for row in selected])
    low = np.array([float(ci_by_case[(origin, row["case_name"])]["H_H0_CI95_low_nm"])*1000
                    for row in selected])
    high = np.array([float(ci_by_case[(origin, row["case_name"])]["H_H0_CI95_high_nm"])*1000
                     for row in selected])
    ax.axhline(0, color="0.5", linestyle="--", linewidth=0.8)
    ax.errorbar(h, delta, yerr=np.vstack((delta-low, high-delta)), fmt="o-",
                color=color, markersize=4, capsize=3, linewidth=1.2)
    ax.set_title(title, loc="left", fontsize=10, pad=9)
    ax.set(xlabel=r"Quadratic coefficient, H (T/m$^2$)", xticks=h,
           xlim=(-0.3, 5.1), ylim=(-1.05, 1.1))
    ax.grid(axis="y", color="0.9", linewidth=0.5)
axes[0].set_ylabel(r"FW99 change, $\delta_H$ ($10^{-3}$ nm)")
fig.text(0.5, 0.02, "Same 50,000 particle identities; 1000 paired bootstrap resamples; "
         "95% sampling intervals", ha="center", fontsize=9)
fig.subplots_adjust(left=0.1, right=0.985, top=0.87, bottom=0.25, wspace=0.19)
save(fig, "fig_quadratic_field_sensitivity_20261004")
print("Four figure files saved; delta_H arithmetic verified for both field baselines.")
