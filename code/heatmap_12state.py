"""Cross-dataset heatmap at AraENCODE 12-state resolution.

Re-runs state enrichment for each DE dataset against the AraENCODE 12-state
TSS assignments (instead of the 5 collapsed groups) and regenerates the
cross-dataset heatmap.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys
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

import chromdecode as cd

OUT = f"{ROOT}/results"
SWEEP = f"{OUT}/sweep"
FIG = f"{ROOT}/figures"

tss = pd.read_csv(f"{OUT}/gene_state_assignments.csv")
st = tss.dropna(subset=["araencode_state"]).rename(
    columns={"araencode_state": "st"})[["agi", "st"]]
bg = st["agi"].unique()

summary = pd.read_csv(f"{SWEEP}/sweep_summary.csv")
rows = []
for acc in summary[summary["status"] == "ok"]["accession"]:
    de_fn = f"{SWEEP}/{acc}_DE.csv"
    if not os.path.exists(de_fn):
        continue
    de = pd.read_csv(de_fn, index_col=0)
    for kind in ["up", "down"]:
        genes = de[de["direction"] == kind].index
        if len(genes) == 0:
            continue
        e = cd.state_enrichment(st, genes, background=bg, state_col="st")
        e.to_csv(f"{SWEEP}/{acc}_enrichment12_{kind}.csv", index=False)
        for _, er in e.iterrows():
            rows.append({"accession": acc, "kind": kind,
                         "state": er["st"], "odds_ratio": er["odds_ratio"],
                         "padj": er["padj"], "n_set": er["n_in_set"]})
        print(f"{acc} {kind}: top={e.iloc[0]['st']} "
              f"odds={e.iloc[0]['odds_ratio']:.2f} FDR={e.iloc[0]['padj']:.1e}")
long = pd.DataFrame(rows)
long.to_csv(f"{OUT}/sweep_enrichment12_long.csv", index=False)

# ---- heatmap: up-regulated, 12-state resolution ----------------------------
up = long[long["kind"] == "up"].pivot_table(
    index="accession", columns="state", values="odds_ratio")
up_fdr = long[long["kind"] == "up"].pivot_table(
    index="accession", columns="state", values="padj")
state_order = [s for s in ["Active", "Transcription1", "Transcription2",
                           "Bivalent", "Repression", "Flanking",
                           "Quiescent", "Heterochromatin"]
               if s in up.columns]
up = up[state_order].reindex(sorted(up.index))
up_fdr = up_fdr[state_order].reindex(up.index)

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

fig, ax = plt.subplots(figsize=(10, max(3, 0.5 * len(up) + 1.5)))
sns.heatmap(up, annot=annot, fmt="", cmap="RdBu_r", center=1,
            vmin=0, vmax=3, ax=ax, linewidths=0.5,
            cbar_kws={"label": "Odds ratio (up-regulated DEGs)"})
ax.set_title("AraENCODE 12-state enrichment of up-regulated DEGs\n"
             "across OSDR Arabidopsis datasets")
ax.set_ylabel("")
ax.set_xlabel("")
fig.savefig(f"{FIG}/fig13_sweep_heatmap_12state_up.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig13_sweep_heatmap_12state_up.png", dpi=150,
            bbox_inches="tight")

# ---- down-regulated heatmap ------------------------------------------------
dn = long[long["kind"] == "down"].pivot_table(
    index="accession", columns="state", values="odds_ratio")
dn_fdr = long[long["kind"] == "down"].pivot_table(
    index="accession", columns="state", values="padj")
dn = dn[state_order].reindex(sorted(dn.index))
dn_fdr = dn_fdr[state_order].reindex(dn.index)
if len(dn):
    annot = dn.copy().astype(object)
    for i in range(dn.shape[0]):
        for j in range(dn.shape[1]):
            o = dn.iloc[i, j]
            q = dn_fdr.iloc[i, j]
            if np.isnan(o):
                annot.iloc[i, j] = ""
                continue
            star = ""
            if not np.isnan(q):
                star = "***" if q < 1e-3 else ("**" if q < 1e-2 else
                                               ("*" if q < 0.05 else ""))
            annot.iloc[i, j] = f"{o:.1f}{star}"
    fig, ax = plt.subplots(figsize=(10, max(3, 0.5 * len(dn) + 1.5)))
    sns.heatmap(dn, annot=annot, fmt="", cmap="RdBu_r", center=1,
                vmin=0, vmax=3, ax=ax, linewidths=0.5,
                cbar_kws={"label": "Odds ratio (down-regulated DEGs)"})
    ax.set_title("AraENCODE 12-state enrichment of down-regulated DEGs\n"
                 "across OSDR Arabidopsis datasets")
    ax.set_ylabel("")
    ax.set_xlabel("")
    fig.savefig(f"{FIG}/fig14_sweep_heatmap_12state_down.svg",
                bbox_inches="tight")
    fig.savefig(f"{FIG}/fig14_sweep_heatmap_12state_down.png", dpi=150,
                bbox_inches="tight")
print("Done.")
