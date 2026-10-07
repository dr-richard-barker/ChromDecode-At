"""Regenerate fig25 with legend outside plot area."""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys
import pickle
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

with open(f"{MOTIF}/motif_targets.pkl", "rb") as f:
    saved = pickle.load(f)
targets_idx, names = saved["targets"], saved["names"]

# rebuild gene_list to map indices -> AGIs (same order as main script)
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

tss = pd.read_csv(f"{OUT}/gene_state_assignments.csv")
state_of = tss.set_index("agi")["pcsd_state"]
POLY = set(range(11, 16))
de37 = pd.read_csv(f"{OUT}/sweep/OSD-37_DE.csv", index_col=0)
down37 = set(de37[de37["direction"] == "down"].index)
ns37 = set(de37[de37["direction"] == "ns"].index)
target = {g for g in down37 if state_of.get(g) in POLY}
control1 = {g for g in down37 if state_of.get(g) not in POLY}
control2 = {g for g in ns37 if state_of.get(g) in POLY}

ov = pd.DataFrame({
    r["name"]: [len(motif_targets[r["motif"]] & s)
                for s in (target, control1, control2)]
    for _, r in top.iterrows()},
    index=["Polycomb-down", "non-Polycomb-down", "Polycomb-non-DE"]).T

fig, ax = plt.subplots(figsize=(7.5, max(3, 0.3 * len(ov))))
left = np.zeros(len(ov))
for col, color in zip(ov.columns, ["#FD9BED", "#0279EE", "#9e9e9e"]):
    ax.barh(ov.index, ov[col], left=left, color=color, label=col)
    left += ov[col].values
ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
ax.set_xlabel("Genes with motif hit in promoter")
fig.savefig(f"{FIG}/fig25_motif_overlap.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig25_motif_overlap.png", dpi=150, bbox_inches="tight")
print("fig25 regenerated")
