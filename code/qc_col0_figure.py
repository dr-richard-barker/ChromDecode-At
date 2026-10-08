"""fig40: power curve and OSD-120 module bound, pooled vs Col-0-only OSD-218.

Reads Tables S43, S65 and S67 (qc_col0_power_bound.py). No values are typed in.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
T = f"{ROOT}/supplementary_tables"

pw = pd.read_csv(f"{T}/qc_power_summary_col0.csv")
bd = pd.read_csv(f"{T}/qc_osd120_module_bound_col0.csv")
fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1.25, 1]})

ax = axes[0]
style = {("pooled", 0): ("#9e9e9e", "-", "o"), ("pooled", 2): ("#9e9e9e", "--", "s"),
         ("col0", 0): ("#0072B2", "-", "o"), ("col0", 2): ("#0072B2", "--", "s")}
for pool, sub in pw.groupby("pool"):
    key = "col0" if "Col-0 only" in pool else "pooled"
    short = "Col-0 only (8 v 8)" if key == "col0" else "pooled Col-0 + WS (as published)"
    for nuis, s in sub.groupby("nuisance"):
        c, ls, mk = style[(key, nuis)]
        s = s.sort_values("k")
        ax.plot(s.k, 100 * s.power, ls=ls, marker=mk, color=c,
                label=f"{short}{' + 2 nuisance covariates' if nuis else ''}")
ax.axhline(80, color="k", lw=0.8, ls=":")
ax.text(16, 82, "80% power", ha="right", fontsize=8)
ax.annotate("k = 8 uses every Col-0 sample\n(identical draws, not a sampling estimate)",
            xy=(8, 100), xytext=(8.6, 40), fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.8))
ax.set_xticks([2, 3, 4, 6, 8, 12, 16])
ax.set_xlabel("Replicates per group (k), subsampled from OSD-218 roots")
ax.set_ylabel("Power: P(Polycomb padj < 0.05), %")
ax.set_ylim(-3, 108)
ax.set_title("A  Empirical power, 20 draws per k", loc="left", fontsize=11)
ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)

ax = axes[1]
order = ["redox_module", "pc_down_module_pooled_as_published", "pc_down_module_col0_only", "pc_down_module_ws_only"]
names = {"redox_module": "redox module", "pc_down_module_pooled_as_published": "Polycomb-down\n(pooled; published)",
         "pc_down_module_col0_only": "Polycomb-down\n(Col-0 only)", "pc_down_module_ws_only": "Polycomb-down\n(WS only)"}
b = bd.set_index("module").loc[order]
x = np.arange(len(order))
ax.bar(x - 0.18, b.median_strat_lfc_Alight, 0.36, color="#E69F00", label="flight effect within Alight")
ax.bar(x + 0.18, b.median_strat_lfc_dark, 0.36, color="#56B4E9", label="flight effect within dark")
ax.axhline(0, color="k", lw=0.8)
for xi, r in zip(x, b.itertuples()):
    ax.text(xi, max(r.median_strat_lfc_Alight, r.median_strat_lfc_dark, 0) + 0.05,
            f"n = {r.n_genes}\nsig. interactions: {r.inter_padj_lt05}\nmedian |int.| {r.inter_median_abs_lfc:.2f}",
            ha="center", va="bottom", fontsize=7.5)
ax.set_xticks(x, [names[m] for m in order], fontsize=8.5)
ax.set_ylabel("Median log2 fold-change of module genes")
ax.set_ylim(-0.7, 1.05)
ax.set_title("B  OSD-120 Col-0 WT: flight effect per light stratum", loc="left", fontsize=11)
ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
fig.tight_layout()
for ext in ("svg", "png"):
    fig.savefig(f"{ROOT}/figures/fig40_power_bound_col0.{ext}", dpi=200, bbox_inches="tight")
print("saved fig40")
