"""fig38: what the 2026-10 QC re-runs changed (original vs corrected contrasts).

Reads only the QC tables written by qc_genotype_rerun.py and qc_sweep_rerun.py
(Tables S48, S51, S53) plus the original sweep tables (S12, S51).
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
FIG = f"{ROOT}/figures"

g = pd.read_csv(f"{T}/qc_genotype_rerun.csv").set_index(["accession", "analysis"])
s5 = pd.read_csv(f"{T}/sweep_enrichment5_corrected.csv")
d314 = pd.read_csv(f"{T}/qc_osd314_design_rerun.csv").set_index("contrast")
s36 = pd.read_csv(f"{T}/sweep_enrichment36_long.csv")


def pc5(acc, kind):
    r = s5[(s5.accession == acc) & (s5.kind == kind) &
           (s5.state_group == "Polycomb_repressed")].iloc[0]
    return r.odds_ratio, r.padj


def stars(p):
    return "ns" if not np.isfinite(p) or p >= 0.05 else \
        "***" if p < 1e-3 else "**" if p < 1e-2 else "*"


# panel A rows: (label, original (odds, padj), corrected [(label, odds, padj)])
rows = [
    ("OSD-314\nseedlings", pc5("OSD-314", "down"),
     [("0g vs 1g\n+ light", d314.loc["0g_vs_1g_light", "polycomb5_down_odds"],
       d314.loc["0g_vs_1g_light", "polycomb5_down_padj"]),
      ("0.3g vs 1g\n+ light", d314.loc["03g_vs_1g_light", "polycomb5_down_odds"],
       d314.loc["03g_vs_1g_light", "polycomb5_down_padj"])]),
    ("OSD-218\nroots", tuple(g.loc[("OSD-218", "v8_as_published"),
                                   ["pc_odds_down", "pc_padj_down"]]),
     [("Col-0\nonly", *g.loc[("OSD-218", "col0_only"), ["pc_odds_down", "pc_padj_down"]]),
      ("WS\nonly", *g.loc[("OSD-218", "ws_only"), ["pc_odds_down", "pc_padj_down"]])]),
    ("OSD-406\nsuborbital roots", tuple(g.loc[("OSD-406", "v8_as_published"),
                                             ["pc_odds_down", "pc_padj_down"]]),
     [("+ genotype\n& rocket", *g.loc[("OSD-406", "all_genotype_rocket_covariates"),
                                    ["pc_odds_down", "pc_padj_down"]])]),
]

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8),
                         gridspec_kw={"width_ratios": [2.3, 1]})
ax = axes[0]
x = 0
ticks, labels = [], []
for label, (o0, p0), corr in rows:
    xs = [x] + [x + 0.9 * (i + 1) for i in range(len(corr))]
    vals = [o0] + [c[1] for c in corr]
    ps = [p0] + [c[2] for c in corr]
    cols = ["#9e9e9e"] + ["#0072B2"] * len(corr)
    ax.bar(xs, vals, width=0.8, color=cols)
    for xi, v, p, name in zip(xs, vals, ps, ["as\npublished"] + [c[0] for c in corr]):
        ax.text(xi, v + 0.08, f"{v:.2f}\n{stars(p)}", ha="center", va="bottom", fontsize=9)
        ax.text(xi, -0.18, name, ha="center", va="top", fontsize=8.5, rotation=0,
                wrap=True)
    ticks.append(np.mean(xs))
    labels.append(label)
    x = xs[-1] + 1.6
ax.axhline(1, ls="--", lw=0.8, color="k")
ax.set_xticks(ticks)
ax.set_xticklabels(labels, fontsize=10)
ax.tick_params(axis="x", pad=40, length=0)
ax.set_ylim(0, 4.9)
ax.set_ylabel("Polycomb-state odds ratio\n(down-regulated genes)")
ax.set_title("A  Down-gene Polycomb enrichment: published vs corrected contrast",
             fontsize=11, loc="left")

ax = axes[1]
o0 = s36[(s36.accession == "OSD-314") & (s36.kind == "up") & (s36.state == 11)].iloc[0]
vals = [o0.odds_ratio, d314.loc["0g_vs_1g_light", "S11_up_odds"],
        d314.loc["03g_vs_1g_light", "S11_up_odds"]]
ps = [o0.padj, d314.loc["0g_vs_1g_light", "S11_up_padj"],
      d314.loc["03g_vs_1g_light", "S11_up_padj"]]
names = ["as published\n(0g+0.3g pooled)", "0g vs 1g\n+ light", "0.3g vs 1g\n+ light"]
ax.bar(range(3), vals, color=["#9e9e9e", "#0072B2", "#0072B2"], width=0.7)
for i, (v, p) in enumerate(zip(vals, ps)):
    ax.text(i, v + 0.3, f"{v:.1f}\n{stars(p)}", ha="center", va="bottom", fontsize=9)
ax.axhline(1, ls="--", lw=0.8, color="k")
ax.set_xticks(range(3))
ax.set_xticklabels(names, fontsize=8.5)
ax.set_ylim(0, vals[0] * 1.25)
ax.set_ylabel("PCSD S11 odds ratio (up-regulated genes)")
ax.set_title("B  OSD-314 up genes in S11", fontsize=11, loc="left")

fig.tight_layout()
for ext in ("svg", "png"):
    fig.savefig(f"{FIG}/fig38_qc_corrections.{ext}", dpi=150, bbox_inches="tight")
print("saved fig38_qc_corrections")
