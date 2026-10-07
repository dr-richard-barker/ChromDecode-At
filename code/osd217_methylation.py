"""OSD-217 spaceflight methylation of the redox module.

Downloads processed per-gene WGBS tables (FLT vs GC, root + leaf,
CG/CHG/CHH), tests whether redox-module genes (Polycomb-down + redox GO)
change methylation vs background (Mann-Whitney, BH-FDR). fig28-29.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import gzip
import io
import os
import sys
import tarfile
import urllib.request

sys.path.insert(0, f"{ROOT}/code")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
from scipy import stats as sps

import chromdecode as cd

matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

OUT = f"{ROOT}/results"
FIG = f"{ROOT}/figures"
WGBS = f"{ROOT}/data/osd217"
os.makedirs(WGBS, exist_ok=True)
MIN_SITES = 3

TISSUES = {"ROOT": "root", "LEAF": "leaf"}
CONTEXTS = ["CG", "CHG", "CHH"]

module = set(pd.read_csv(f"{OUT}/redox_module_genes.csv")["agi"])
print(f"redox module: {len(module)} genes")

# ---------------------------------------------------------------- download
frames = []
for tissue in TISSUES:
    for chunk in range(1, 6):
        fn = f"GLDS-217_wgbs_GSE95594_{tissue}FT-vs-{tissue}GC_genes0" \
             f"{chunk}_full_tar.gz"
        dest = f"{WGBS}/{fn}"
        if not os.path.exists(dest):
            url = ("https://osdr.nasa.gov/geode-py/ws/studies/OSD-217/"
                   f"download?source=datamanager&file={fn}")
            print(f"  downloading {fn} ...")
            urllib.request.urlretrieve(url, dest)
        with tarfile.open(dest, "r:gz") as tf:
            for member in tf.getmembers():
                if not member.name.endswith(".xlsx"):
                    continue
                ctx = "CG" if ".CG." in member.name else \
                    "CHG" if ".CHG." in member.name else "CHH"
                df = pd.read_excel(
                    io.BytesIO(tf.extractfile(member).read()))
                df.columns = ["agi", "mrna", "chrom", "start", "end",
                              "sites", "diff"]
                df["context"] = ctx
                df["tissue"] = TISSUES[tissue]
                frames.append(df[["agi", "sites", "diff", "context",
                                  "tissue"]])
meth = pd.concat(frames, ignore_index=True)
meth = meth[meth["sites"] >= MIN_SITES]
print(f"methylation table: {len(meth):,} gene-context-tissue rows")

# ---------------------------------------------------------------- tests
rows = []
for tissue in TISSUES.values():
    for ctx in CONTEXTS:
        sub = meth[(meth["tissue"] == tissue) & (meth["context"] == ctx)]
        mod = sub[sub["agi"].isin(module)]["diff"]
        bg = sub[~sub["agi"].isin(module)]["diff"]
        u, p = sps.mannwhitneyu(mod, bg, alternative="two-sided")
        rows.append({"tissue": tissue, "context": ctx,
                     "n_module": len(mod), "n_background": len(bg),
                     "median_module": mod.median(),
                     "median_background": bg.median(),
                     "p_value": p})
tests = pd.DataFrame(rows)
tests["padj"] = cd.bh_adjust(tests["p_value"].values)
tests.to_csv(f"{OUT}/osd217_methylation_tests.csv", index=False)
print(tests.to_string(index=False))

# ---------------------------------------------------------------- fig28
fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
for ax, ctx in zip(axes, CONTEXTS):
    data, labels = [], []
    for tissue in ["root", "leaf"]:
        sub = meth[(meth["tissue"] == tissue) & (meth["context"] == ctx)]
        data.append(sub[~sub["agi"].isin(module)]["diff"])
        labels.append(f"{tissue} bg (n={len(data[-1]):,})")
        data.append(sub[sub["agi"].isin(module)]["diff"])
        labels.append(f"{tissue} module (n={len(data[-1])})")
    bp = ax.boxplot(data, tick_labels=labels, showfliers=False)
    for i in (1, 3):
        bp["medians"][i].set_color("#FD9BED")
        bp["medians"][i].set_linewidth(2)
    t = tests[(tests["context"] == ctx)]
    for k, (_, r) in enumerate(t.iterrows()):
        star = "**" if r["padj"] < 0.01 else ("*" if r["padj"] < 0.05
                                              else "n.s.")
        ax.text(k * 2 + 1.5, ax.get_ylim()[1] * 0.95, star, ha="center")
    ax.set_title(ctx)
    ax.tick_params(axis="x", rotation=45, labelsize=7)
axes[0].set_ylabel("Methylation difference FLT - GC")
fig.suptitle("Redox module vs background, spaceflight methylation "
             "(OSD-217)", fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIG}/fig28_methylation_module.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig28_methylation_module.png", dpi=150,
            bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- fig29
mod_meth = meth[meth["agi"].isin(module)]
pivot = mod_meth.pivot_table(index="agi",
                             columns=["tissue", "context"], values="diff")
lfc_col = [c for c in pd.read_csv(f"{OUT}/sweep/OSD-37_DE.csv",
                                  index_col=0).columns
           if "log2" in c.lower() and "fc" in c.lower()][0]
de37 = pd.read_csv(f"{OUT}/sweep/OSD-37_DE.csv", index_col=0)
order = pivot.index.to_series().map(de37[lfc_col]).sort_values().index
pivot = pivot.loc[order]
fig, ax = plt.subplots(figsize=(6, max(4, 0.22 * len(pivot))))
im = ax.imshow(pivot.values, aspect="auto", cmap="coolwarm",
               vmin=-0.3, vmax=0.3)
ax.set_xticks(range(len(pivot.columns)),
              [f"{t}\n{c}" for t, c in pivot.columns], fontsize=7)
ax.set_yticks(range(len(pivot)),
              [f"{g}" for g in pivot.index], fontsize=5)
plt.colorbar(im, ax=ax, label="FLT - GC methylation")
ax.set_title("Per-gene methylation change, redox module "
             "(rows sorted by log2FC)", fontsize=9)
fig.savefig(f"{FIG}/fig29_methylation_heatmap.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig29_methylation_heatmap.png", dpi=150,
            bbox_inches="tight")
print("fig28-29 done")
