"""SOG1 ChIP-seq promoter-overlap analysis and figures (fig33-34)."""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys
from collections import defaultdict

sys.path.insert(0, f"{ROOT}/code")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

import chromdecode as cd

OUT = f"{ROOT}/results"
FIG = f"{ROOT}/figures"
W = f"{ROOT}/data/osd496"
UPSTREAM = 2000

# ---------------------------------------------------------------- inputs
genes = cd.load_gene_coordinates(
    f"{ROOT}/data/araencode/genes_TAIR10.txt")
genes["tss"] = np.where(genes["strand"] == "+", genes["start"], genes["end"])
tss = pd.read_csv(f"{OUT}/gene_state_assignments.csv")
state_of = tss.set_index("agi")["pcsd_state"]
POLY = set(range(11, 16))
de37 = pd.read_csv(f"{OUT}/sweep/OSD-37_DE.csv", index_col=0)
down37 = set(de37[de37["direction"] == "down"].index)
ns37 = set(de37[de37["direction"] == "ns"].index)
module = set(pd.read_csv(f"{OUT}/redox_module_genes.csv")["agi"])

sets = {
    "redox_module": module,
    "polycomb_down": {g for g in down37 if state_of.get(g) in POLY},
    "nonpolycomb_down": {g for g in down37 if state_of.get(g) not in POLY},
    "background": {g for g in ns37 if state_of.get(g) in POLY},
}
print({k: len(v) for k, v in sets.items()})

def load_peaks(path):
    pk = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end", "name", "score",
                            "strand", "fc", "pval", "qval", "summit"])
    # HISAT2 index built from Ensembl-style headers ("1".."5")
    pk["chrom"] = pk["chrom"].map(
        lambda c: f"chr{c}" if str(c).isdigit() else c)
    return pk[pk["chrom"].isin(["chr1", "chr2", "chr3", "chr4", "chr5"])]

peaks = {
    "20min": load_peaks(f"{W}/peaks/SOG1_20min_peaks.narrowPeak"),
    "1h": load_peaks(f"{W}/peaks/SOG1_1h_peaks.narrowPeak"),
    "vsWT_20min": load_peaks(f"{W}/peaks/SOG1vsWT_20min_peaks.narrowPeak"),
    "vsWT_1h": load_peaks(f"{W}/peaks/SOG1vsWT_1h_peaks.narrowPeak"),
}
union = pd.concat([peaks["20min"], peaks["1h"]]).drop_duplicates(
    ["chrom", "start", "end"])
print("union peaks:", len(union))

# vectorized overlap: for each gene, promoter interval [tss-2000, tss)
# (+ strand) or [tss, tss+2000) (- strand)
def promoter_peak_hits(pk):
    hits = set()
    pk_by_chrom = {c: (g["start"].to_numpy(), g["end"].to_numpy())
                   for c, g in pk.groupby("chrom")}
    for r in genes.itertuples():
        if r.strand == "+":
            ps, pe = r.tss - UPSTREAM, r.tss
        else:
            ps, pe = r.tss, r.tss + UPSTREAM
        se = pk_by_chrom.get(r.chr)
        if se is None:
            continue
        starts, ends = se
        lo = np.searchsorted(ends, ps, side="right")
        hi = np.searchsorted(starts, pe, side="left")
        if hi > lo:
            hits.add(r.agi)
    return hits

# ---------------------------------------------------------------- tests
rows = []
for label, pk in [("20min", peaks["20min"]), ("1h", peaks["1h"]),
                  ("union", union), ("vsWT_union",
                                    pd.concat([peaks["vsWT_20min"],
                                               peaks["vsWT_1h"]])
                                    .drop_duplicates(
                                        ["chrom", "start", "end"]))]:
    hits = promoter_peak_hits(pk)
    for set_name, gset in sets.items():
        bg = set(genes["agi"]) - gset
        a = len(gset & hits)
        b = len(gset) - a
        c = len(bg & hits)
        d = len(bg) - c
        o, p = sps.fisher_exact([[a, b], [c, d]])
        rows.append({"peaks": label, "gene_set": set_name,
                     "n_genes": len(gset), "n_with_peak": a,
                     "odds_ratio": o, "p_value": p})
tests = pd.DataFrame(rows)
ok = tests["p_value"].notna()
tests.loc[ok, "padj"] = cd.bh_adjust(tests.loc[ok, "p_value"].values)
tests.to_csv(f"{OUT}/sog1_promoter_tests.csv", index=False)
print(tests.to_string(index=False))

# ---------------------------------------------------------------- fig33
fig, ax = plt.subplots(figsize=(8, 5))
labels_order = ["redox_module", "polycomb_down", "nonpolycomb_down",
                "background"]
width = 0.2
for k, (label, color) in enumerate([("20min", "#FD9BED"),
                                    ("1h", "#0279EE"),
                                    ("union", "#75A025"),
                                    ("vsWT_union", "#9e9e9e")]):
    sub = tests[tests["peaks"] == label].set_index("gene_set")
    vals = [sub.loc[s, "odds_ratio"] if s in sub.index else np.nan
            for s in labels_order]
    xs = np.arange(len(labels_order)) + (k - 1.5) * width
    ax.bar(xs, vals, width, label=label, color=color)
    for x, s, v in zip(xs, labels_order, vals):
        q = sub.loc[s, "padj"] if s in sub.index else 1
        if not np.isnan(v) and q < 0.05:
            ax.text(x, v + 0.05, "**" if q < 0.01 else "*", ha="center",
                    fontsize=8)
ax.axhline(1, color="k", lw=0.8, ls="--")
ax.set_xticks(range(len(labels_order)),
              ["redox module", "Polycomb-down", "non-Polycomb-down",
               "Polycomb non-DE"])
ax.set_ylabel("Odds ratio of SOG1 peak in promoter")
ax.legend(frameon=False, fontsize=8)
fig.savefig(f"{FIG}/fig33_sog1_enrichment.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig33_sog1_enrichment.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("fig33 done")

# ---------------------------------------------------------------- fig34
# SOG1 signal per 100-bp promoter bin (IP reads, scaled by library size)
import subprocess
def bin_coverage(bam, chrom, starts, ends):
    cmd = f"samtools view -c {bam} {chrom}:{starts}-{ends}"
    return int(subprocess.check_output(cmd, shell=True).strip())

libsize = {n: int(subprocess.check_output(
    f"samtools view -c {W}/bam/{n}.bam", shell=True))
    for n in ["IP_SOG1_20min", "IP_SOG1_1h"]}

def occupancy(gene_set, sample):
    scale = 1e6 / libsize[sample]
    mat = []
    agis = []
    for r in genes[genes["agi"].isin(gene_set)].itertuples():
        row = []
        for b in range(20):
            if r.strand == "+":
                s = r.tss - UPSTREAM + b * 100
            else:
                s = r.tss + b * 100
            row.append(bin_coverage(f"{W}/bam/{sample}.bam",
                                    r.chr.replace("chr", ""),
                                    s, s + 99) * scale)
        mat.append(row)
        agis.append(r.agi)
    return np.array(mat), agis

mat_mod, agis_mod = occupancy(sets["redox_module"], "IP_SOG1_1h")
mat_pc, agis_pc = occupancy(sets["polycomb_down"] - module, "IP_SOG1_1h")
mat_bg, agis_bg = occupancy(list(sets["background"])[:300], "IP_SOG1_1h")
mat_all = np.vstack([mat_mod, mat_pc, mat_bg])
classes = (["redox module"] * len(mat_mod) +
           ["Polycomb-down"] * len(mat_pc) + ["background"] * len(mat_bg))
order = np.argsort([c for c in classes])
mat_all = mat_all[order]
classes = [classes[i] for i in order]
vmax = np.quantile(mat_all, 0.99)
fig, ax = plt.subplots(figsize=(7, max(4, 0.06 * len(mat_all))))
im = ax.imshow(mat_all, aspect="auto", cmap="viridis", vmin=0,
               vmax=max(vmax, 1))
bounds = [0, len(mat_mod), len(mat_mod) + len(mat_pc), len(mat_all)]
for b in bounds[1:-1]:
    ax.axhline(b, color="w", lw=1)
ax.set_yticks([])
ax.set_xticks([0, 10, 19], ["-2 kb", "-1 kb", "TSS"])
ax.set_xlabel("Promoter bin (SOG1 IP, 1 h, RPM)")
cls_colors = {"redox module": "#FD9BED", "Polycomb-down": "#0279EE",
              "background": "#9e9e9e"}
for cname, y0, y1 in [("redox module", bounds[0], bounds[1]),
                      ("Polycomb-down", bounds[1], bounds[2]),
                      ("background", bounds[2], bounds[3])]:
    ax.text(-0.02, (y0 + y1) / 2 / len(mat_all), cname,
            transform=ax.get_yaxis_transform(), ha="right", va="center",
            fontsize=7, color=cls_colors[cname])
plt.colorbar(im, ax=ax, label="RPM per 100-bp bin")
fig.savefig(f"{FIG}/fig34_sog1_occupancy.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig34_sog1_occupancy.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("fig34 done")
