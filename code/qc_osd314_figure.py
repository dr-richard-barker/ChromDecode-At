"""fig39: OSD-314 diagnostics on the published vs corrected contrasts.

Reads Tables S58, S60, S63 (written by qc_osd314_diagnostics.py) and, for the
running-enrichment curves, the three OSD-314 DE tables (S55, S56, S57) plus the
gene-state map (D1). No values are typed in.
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

LAB = {"published": "published (0g + 0.3g pooled)", "03g_light": "0.3g vs 1g + light",
       "0g_light": "0g vs 1g + light"}
COL = {"published": "#9e9e9e", "03g_light": "#0072B2", "0g_light": "#56B4E9"}
ORDER = ["published", "03g_light", "0g_light"]

null = pd.read_csv(f"{T}/qc_osd314_null_calibration_corrected.csv")
grid = pd.read_csv(f"{T}/qc_osd314_threshold_grid_corrected.csv")
go = pd.read_csv(f"{T}/qc_osd314_go_replication_corrected.csv")
tss = pd.read_csv(f"{ROOT}/data/derived/gene_state_assignments.csv")
poly = set(tss.loc[tss.pcsd_state.isin([11, 12, 13, 15]), "agi"])
DE = {"published": f"{T}/OSD-314_DE_sweep_miscoded.csv",
      "0g_light": f"{T}/OSD-314_0g_vs_1g_light_DE.csv",
      "03g_light": f"{T}/OSD-314_03g_vs_1g_light_DE.csv"}

fig, axes = plt.subplots(2, 2, figsize=(13, 10.5))

# A: observed odds vs null 99th percentile, per state
ax = axes[0, 0]
states = ["S11", "S12", "S13", "S15"]
w = 0.26
for k, c in enumerate(ORDER):
    sub = null[null.contrast == c].set_index("state").loc[states]
    x = np.arange(len(states)) + (k - 1) * w
    ax.bar(x, sub.observed_odds, width=w, color=COL[c], label=LAB[c])
    ax.scatter(x, sub.null_p99, marker="_", s=180, color="k", zorder=3,
               label="null 99th percentile" if k == 0 else None)
    for xi, (o, p) in zip(x, zip(sub.observed_odds, sub.empirical_p)):
        ax.text(xi, o * 1.08, f"{o:.1f}", ha="center", va="bottom", fontsize=8)
ax.set_yscale("log")
ax.set_ylim(0.5, 60)
ax.set_xticks(range(len(states)), states)
ax.set_xlabel("PCSD Polycomb state")
ax.set_ylabel("Odds ratio, up-regulated genes (log scale)")
ax.set_title("A  Null calibration (10,000 random gene sets of equal size)", loc="left", fontsize=11)
ax.legend(frameon=False, fontsize=8, loc="upper right")

# B: threshold dose-response for S11 and S15 at padj < 0.05
ax = axes[0, 1]
for c in ORDER:
    sub = grid[(grid.contrast == c) & (grid.padj == 0.05)].sort_values("lfc")
    ax.plot(sub.lfc, sub.S11, "-o", color=COL[c], label=f"{LAB[c]}: S11")
    ax.plot(sub.lfc, sub.S15, "--s", color=COL[c], label=f"{LAB[c]}: S15")
ax.axhline(1, color="k", lw=0.8, ls=":")
ax.set_xticks([0.5, 1, 2])
ax.set_xlabel("|log2 fold-change| threshold (padj < 0.05)")
ax.set_ylabel("Odds ratio, up-regulated genes")
ax.set_title("B  Dose-response with DE stringency", loc="left", fontsize=11)
ax.legend(frameon=False, fontsize=7.5, ncol=1, loc="upper left")

# C: running enrichment of Polycomb genes along the log2FC ranking
ax = axes[1, 0]
for c in ORDER:
    de = pd.read_csv(DE[c], index_col=0)
    ranked = de["log2fc"].sort_values(ascending=False)
    hit = ranked.index.isin(poly)
    run = np.cumsum(np.where(hit, 1 / hit.sum(), -1 / (len(hit) - hit.sum())))
    ax.plot(np.linspace(0, 1, len(run)), run, color=COL[c], lw=1.4, label=LAB[c])
ax.axhline(0, color="k", lw=0.6)
ax.set_xlabel("Gene rank by log2 fold-change (fraction; up-regulated at left)")
ax.set_ylabel("Running enrichment sum")
ax.set_title(f"C  Threshold-free enrichment of {len(poly & set(ranked.index)):,} Polycomb genes",
             loc="left", fontsize=11)
ax.legend(frameon=False, fontsize=8)

# D: GO replication of the OSD-37 redox module
ax = axes[1, 1]
terms = go[go.contrast == "published"]["name"].tolist()
short = [t if len(t) <= 32 else t[:30] + "…" for t in terms]
y = np.arange(len(terms))
for k, c in enumerate(ORDER):
    sub = go[go.contrast == c].set_index("name").loc[terms]
    yy = y + (k - 1) * 0.26
    ax.barh(yy, sub.OSD314_fold, height=0.26, color=COL[c],
            label=f"{LAB[c]} (n = {int(sub.OSD314_polycomb_down_n.iloc[0])} Polycomb-down)")
    for yi, f, q in zip(yy, sub.OSD314_fold, sub.OSD314_padj_within_contrast):
        if q < 0.05:
            ax.text(f * 1.05, yi, "*", va="center", fontsize=9)
ax.set_xscale("log")
ax.set_yticks(y, short, fontsize=8)
ax.invert_yaxis()
ax.axvline(1, color="k", lw=0.8, ls=":")
ax.set_xlabel("Fold enrichment in OSD-314 Polycomb-down genes (log scale)")
ax.set_title("D  Replication of the top 10 OSD-37 GO terms", loc="left", fontsize=11)
ax.legend(frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.45, -0.12), ncol=1,
          title="* BH-adjusted p < 0.05 within contrast", title_fontsize=7.5)

fig.tight_layout()
for ext in ("svg", "png"):
    fig.savefig(f"{FIG}/fig39_osd314_corrected_diagnostics.{ext}", dpi=200, bbox_inches="tight")
print("saved fig39")
