"""Aggregate OSDR sweep results: cross-dataset heatmap + summary tables."""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys
import glob
import os

sys.path.insert(0, f"{ROOT}/code")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

OUT = f"{ROOT}/results"
SWEEP = f"{OUT}/sweep"
FIG = f"{ROOT}/figures"

summary = pd.read_csv(f"{SWEEP}/sweep_summary.csv")

# ---- collect enrichment tables (DE-based up/down + expression-only) --------
rows = []
for _, r in summary.iterrows():
    acc = r["accession"]
    for kind, fn in [("up", f"{SWEEP}/{acc}_enrichment_up.csv"),
                     ("down", f"{SWEEP}/{acc}_enrichment_down.csv"),
                     ("expr_only", f"{SWEEP}/{acc}_expr_only_enrichment.csv")]:
        if os.path.exists(fn):
            e = pd.read_csv(fn)
            for _, er in e.iterrows():
                rows.append({"accession": acc, "kind": kind,
                             "state_group": er["st"],
                             "odds_ratio": er["odds_ratio"],
                             "padj": er["padj"],
                             "n_set": er["n_in_set"]})
long = pd.DataFrame(rows)

# ---- cross-dataset heatmap: up-regulation enrichment -----------------------
up = long[long["kind"] == "up"].pivot_table(
    index="accession", columns="state_group", values="odds_ratio")
up_fdr = long[long["kind"] == "up"].pivot_table(
    index="accession", columns="state_group", values="padj")
order = [c for c in ["Accessible_promoter", "Active_transcribed",
                     "Polycomb_repressed", "Intergenic_quiet",
                     "Heterochromatin_TE"] if c in up.columns]
up = up[order].reindex([a for a in summary["accession"] if a in up.index])
# FDR must share the odds table's row/column order: annotations index both
# positionally (2026-10 QC fix: stars were landing on the wrong cells).
up_fdr = up_fdr.reindex(index=up.index, columns=up.columns)

annot = up.copy().astype(object)
for i in range(up.shape[0]):
    for j in range(up.shape[1]):
        o = up.iloc[i, j]
        q = up_fdr.iloc[i, j]
        if np.isnan(o):
            annot.iloc[i, j] = ""
            continue
        star = ""
        if not np.isnan(q):
            star = "***" if q < 1e-3 else ("**" if q < 1e-2 else
                                           ("*" if q < 0.05 else ""))
        annot.iloc[i, j] = f"{o:.1f}{star}"

fig, ax = plt.subplots(figsize=(8, max(3, 0.45 * len(up) + 2)))
sns.heatmap(up, annot=annot, fmt="", cmap="RdBu_r", center=1,
            vmin=0, vmax=3, ax=ax,
            cbar_kws={"label": "Odds ratio (up-regulated DEGs)"})
ax.set_title("Chromatin-state enrichment of up-regulated DEGs across OSDR datasets")
ax.set_ylabel("")
fig.savefig(f"{FIG}/fig10_sweep_heatmap_up.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig10_sweep_heatmap_up.png", dpi=150, bbox_inches="tight")

# ---- expression-only heatmap ----------------------------------------------
eo = long[long["kind"] == "expr_only"].pivot_table(
    index="accession", columns="state_group", values="odds_ratio")
if len(eo):
    eo = eo[order].reindex([a for a in summary["accession"] if a in eo.index])
    fig, ax = plt.subplots(figsize=(8, max(3, 0.45 * len(eo) + 2)))
    sns.heatmap(eo, cmap="RdBu_r", center=1, vmin=0, vmax=3, ax=ax,
                cbar_kws={"label": "Odds ratio (top-quartile expressed genes)"})
    ax.set_title("Chromatin-state composition of highly expressed genes\n"
                 "(datasets without auto-detectable contrasts)")
    ax.set_ylabel("")
    fig.savefig(f"{FIG}/fig11_sweep_heatmap_expr_only.svg", bbox_inches="tight")
    fig.savefig(f"{FIG}/fig11_sweep_heatmap_expr_only.png", dpi=150,
                bbox_inches="tight")

# ---- per-dataset enrichment bar charts for DE datasets --------------------
import figures as fg
for acc in summary[summary["status"] == "ok"]["accession"]:
    for kind in ["up", "down"]:
        fn = f"{SWEEP}/{acc}_enrichment_{kind}.csv"
        if os.path.exists(fn):
            e = pd.read_csv(fn)
            fg.fig_enrichment(e, "st", {g: g for g in e["st"]},
                              f"{FIG}/fig12_{acc}_enrichment_{kind}",
                              title=f"{acc}: {kind}-regulated DEGs")

# ---- final summary table ----------------------------------------------------
long.to_csv(f"{OUT}/sweep_enrichment_long.csv", index=False)
print("Aggregation done.")
print(summary[summary["status"] == "ok"].to_string(index=False))
