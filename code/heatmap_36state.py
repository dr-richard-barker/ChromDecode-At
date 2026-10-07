"""Cross-dataset heatmap at PCSD 36-state resolution (mark-level detail).

Re-runs state enrichment for each DE dataset against the PCSD 36-state TSS
assignments and generates up/down cross-dataset heatmaps with state
descriptions (epigenetic marks + preferential location) as labels.
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
DATA = f"{ROOT}/data"

tss = pd.read_csv(f"{OUT}/gene_state_assignments.csv")
st = tss.dropna(subset=["pcsd_state"]).rename(
    columns={"pcsd_state": "st"})[["agi", "st"]]
st["st"] = st["st"].astype(int)
bg = st["agi"].unique()

desc = pd.read_csv(f"{DATA}/pcsd/state_descriptions.tsv", sep="\t")
label = {r["state"]: f"S{r['state']}: {r['marks'][:38]}" for _, r in desc.iterrows()}

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
        e.to_csv(f"{SWEEP}/{acc}_enrichment36_{kind}.csv", index=False)
        for _, er in e.iterrows():
            rows.append({"accession": acc, "kind": kind, "state": er["st"],
                         "odds_ratio": er["odds_ratio"], "padj": er["padj"],
                         "n_set": er["n_in_set"]})
        sig = e[e["padj"] < 0.05].sort_values("odds_ratio", ascending=False)
        tops = ", ".join(f"S{int(r['st'])}({r['odds_ratio']:.1f})"
                         for _, r in sig.head(3).iterrows())
        print(f"{acc} {kind} ({len(genes)} genes): significant: {tops or 'none'}")
long = pd.DataFrame(rows)
long.to_csv(f"{OUT}/sweep_enrichment36_long.csv", index=False)

# order states by genomic-mark logic: active marks -> polycomb -> accessible
# -> quiet -> heterochromatin
state_order = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 22, 25, 26, 27, 28,
               11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24,
               29, 30, 31, 32, 33, 34, 35, 36]


def make_heatmap(kind, fname):
    d = long[long["kind"] == kind].pivot_table(
        index="accession", columns="state", values="odds_ratio")
    d_fdr = long[long["kind"] == kind].pivot_table(
        index="accession", columns="state", values="padj")
    cols = [s for s in state_order if s in d.columns]
    d = d[cols].reindex(sorted(d.index))
    d_fdr = d_fdr[cols].reindex(d.index)
    # display: clip odds to [0, 5] for color scale; annotate raw value
    disp = d.clip(0, 5)
    annot = d.copy().astype(object)
    for i in range(d.shape[0]):
        for j in range(d.shape[1]):
            o = d.iloc[i, j]
            q = d_fdr.iloc[i, j]
            if np.isnan(o):
                annot.iloc[i, j] = ""
                continue
            star = ""
            if not np.isnan(q):
                star = ("***" if q < 1e-3 else "**" if q < 1e-2 else
                        "*" if q < 0.05 else "")
            txt = "inf" if np.isinf(o) else f"{o:.1f}"
            annot.iloc[i, j] = f"{txt}{star}"
    fig, ax = plt.subplots(figsize=(13, max(3, 0.42 * len(d) + 1.5)))
    sns.heatmap(disp, annot=annot, fmt="", cmap="RdBu_r", center=1,
                vmin=0, vmax=5, ax=ax, linewidths=0.3,
                cbar_kws={"label": f"Odds ratio ({kind}-regulated DEGs, "
                                    "color clipped at 5)"})
    ax.set_title(f"PCSD 36-state enrichment of {kind}-regulated DEGs\n"
                 "across OSDR Arabidopsis datasets")
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_xticklabels([label.get(int(t.get_text().split(".")[0]),
                                  t.get_text())
                        for t in ax.get_xticklabels()],
                       rotation=45, ha="right", fontsize=7)
    fig.savefig(f"{FIG}/{fname}.svg", bbox_inches="tight")
    fig.savefig(f"{FIG}/{fname}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


make_heatmap("up", "fig15_sweep_heatmap_36state_up")
make_heatmap("down", "fig16_sweep_heatmap_36state_down")
print("Done.")
