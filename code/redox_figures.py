"""fig26: motif-logo heatmap of top enriched/depleted motifs.
fig27: redox-module network (GO + motif-sharing edges).
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import pickle
import sys

sys.path.insert(0, f"{ROOT}/code")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

MOTIF = f"{ROOT}/data/motif_go"
OUT = f"{ROOT}/results"
FIG = f"{ROOT}/figures"
BASE = "#0279EE"
COLORS = {"A": "#ECE9E2", "C": "#FD9BED", "G": "#75A025", "T": "#0279EE"}

# ---------------------------------------------------------------- inputs
with open(f"{MOTIF}/motif_targets.pkl", "rb") as f:
    saved = pickle.load(f)
targets_idx, names_raw = saved["targets"], saved["names"]

import chromdecode as cd
genes = cd.load_gene_coordinates(
    f"{ROOT}/data/araencode/genes_TAIR10.txt")
gene_list = sorted(set(genes["agi"]))
motif_targets = {m: {gene_list[i] for i in h} for m, h in targets_idx.items()}

m = pd.read_csv(f"{OUT}/motif_enrichment.csv")
m["name"] = m["name"].str.replace("\t", " ", regex=False)
e = m[m["comparison"] == "target_vs_universe"]
sig = e[(e["padj"] < 0.05) & (e["n_target_hits"] > 0)]
enriched = sig[sig["odds_ratio"] > 1].sort_values("odds_ratio",
                                                  ascending=False).head(12)
depleted = sig[sig["odds_ratio"] < 1].sort_values("odds_ratio").head(12)
top = pd.concat([enriched, depleted]).drop_duplicates()

def load_pfm(path):
    rows = {}
    for line in open(path):
        if line.startswith(">"):
            continue
        parts = line.replace("[", " ").replace("]", " ").split()
        if len(parts) >= 2 and parts[0] in "ACGT":
            rows[parts[0]] = np.array([float(x) for x in parts[1:]])
    if len(rows) != 4:
        return None
    return np.vstack([rows[b] for b in "ACGT"])

# ---------------------------------------------------------------- fig26
n = len(top)
LABEL_W = 11.0   # data units for label column
LOGO_W = 30.0    # data units for logo column (max motif length)
fig, ax = plt.subplots(figsize=(11, n * 0.42 + 1.0))
for i, (_, r) in enumerate(top.iterrows()):
    y = n - 1 - i
    star = "**" if r["padj"] < 0.01 else ("*" if r["padj"] < 0.05 else "")
    ax.text(-0.5, y + 0.5, f"{r['name']}  ({r['motif']})",
            fontsize=7, va="center", ha="right")
    ax.text(LABEL_W + LOGO_W + 0.5, y + 0.5,
            f"OR {r['odds_ratio']:.1f}{star}",
            fontsize=7, va="center", ha="left",
            color=BASE if r["odds_ratio"] > 1 else "#666666")
    pfm = load_pfm(f"{MOTIF}/pfms/{r['motif']}.jaspar")
    if pfm is None:
        continue
    p = pfm / pfm.sum(axis=0, keepdims=True)
    L = p.shape[1]
    for pos in range(L):
        h = 0.0
        for b in "ACGT":
            pb = p["ACGT".index(b), pos]
            if pb > 0.01:
                ax.add_patch(plt.Rectangle(
                    (LABEL_W + pos, y + h), 1, pb,
                    facecolor=COLORS[b], edgecolor="none", linewidth=0))
                h += pb
ax.set_xlim(-LABEL_W, LABEL_W + LOGO_W + 6)
ax.set_ylim(0, n)
ax.axis("off")
ax.set_title("Top motifs in Polycomb-down promoters "
             "(blue = enriched, gray = depleted; * FDR<0.05, "
             "** FDR<0.01)", fontsize=9, loc="left")
fig.savefig(f"{FIG}/fig26_motif_logos.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig26_motif_logos.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("fig26 done")

# ---------------------------------------------------------------- fig27
import networkx as nx

tss = pd.read_csv(f"{OUT}/gene_state_assignments.csv")
state_of = tss.set_index("agi")["pcsd_state"]
POLY = set(range(11, 16))
de37 = pd.read_csv(f"{OUT}/sweep/OSD-37_DE.csv", index_col=0)
down37 = set(de37[de37["direction"] == "down"].index)

go = pd.read_csv(f"{OUT}/go_enrichment.csv")
gu = go[(go["comparison"] == "target_vs_universe") & (go["padj"] < 0.05)]
REDOX_KEYS = ("peroxidase", "heme", "hydrogen peroxide", "oxidant",
              "monooxygenase", "oxidoreductase", "iron ion", "lactoperoxidase")
redox_terms = gu[gu["name"].str.lower().str.contains("|".join(REDOX_KEYS))]

from collections import defaultdict
import gzip
KEEP_EV = {"IDA", "IMP", "IGI", "IEP", "HDA", "HMP", "HGI", "HEP",
           "TAS", "IC", "EXP", "IEA"}
gene_terms = defaultdict(set)
with gzip.open(f"{MOTIF}/tair.gaf.gz", "rt") as gf:
    for line in gf:
        if line.startswith("!"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 15 or parts[6] not in KEEP_EV:
            continue
        gene_terms[cd.clean_agi(parts[1])].add(parts[4])

module = sorted({g for g in down37 if state_of.get(g) in POLY
                 and gene_terms[g] & set(redox_terms["GO"])})
print(f"redox module: {len(module)} genes")

top_motifs = set(top["motif"])
G = nx.Graph()
for g in module:
    G.add_node(g)
for i, g1 in enumerate(module):
    for g2 in module[i + 1:]:
        shared_go = len(gene_terms[g1] & gene_terms[g2] &
                        set(redox_terms["GO"])) > 0
        shared_motif = any(g2 in motif_targets[mt] for mt in top_motifs
                           if g1 in motif_targets[mt])
        if shared_go or shared_motif:
            G.add_edge(g1, g2)

lfc_col = [c for c in de37.columns if "log2" in c.lower()
           and "fc" in c.lower()][0]
lfc = de37[lfc_col]
deg = dict(G.degree())
pos = nx.spring_layout(G, seed=42, k=0.6)
fig, ax = plt.subplots(figsize=(10, 8))
nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.25, edge_color="#666666")
sc = nx.draw_networkx_nodes(
    G, pos, ax=ax, node_size=[60 + 40 * deg[g] for g in G.nodes],
    node_color=[lfc.get(g, 0) for g in G.nodes], cmap="coolwarm",
    vmin=-6, vmax=6, edgecolors="k", linewidths=0.3)
labels = {g: g for g in G.nodes if deg[g] >= np.quantile(list(deg.values()),
                                                         0.9)
          or abs(lfc.get(g, 0)) > 5}
nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=6)
plt.colorbar(sc, ax=ax, label="log2 fold change (spaceflight, down)")
ax.axis("off")
ax.set_title(f"Redox module: {len(module)} Polycomb-down genes "
             f"({G.number_of_edges()} edges; edge = shared redox GO term "
             "or motif target)", fontsize=9)
fig.savefig(f"{FIG}/fig27_redox_network.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig27_redox_network.png", dpi=150, bbox_inches="tight")
pd.Series(module).to_csv(f"{OUT}/redox_module_genes.csv", index=False,
                         header=["agi"])
print("fig27 done")
