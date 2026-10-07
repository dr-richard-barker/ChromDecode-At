"""OSD-217 promoter-window methylation + ecotype-matched comparison.

- Parse winavg tables (100-bp methylation bins, FLT/GC x root/leaf x
  CG/CHG/CHH), cache to CSV.
- Promoter = 2 kb upstream of TSS (Ws coords + TAIR10 strand) = 20 bins.
- Tests: promoter methylation change (FLT-GC) module vs background,
  Mann-Whitney, BH-FDR. Modules: Col-0 redox module (23 genes) and
  Ws-defined module (down-DEGs in OSD-217 Ws RNA-seq in PCSD Polycomb
  states).
- fig30 metaplots, fig31 ecotype-matched boxplots, fig32 Ws-vs-Col
  concordance scatter.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import os
import sys
import tarfile
import urllib.request

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
WGBS = f"{ROOT}/data/osd217"
CACHE = f"{WGBS}/cache"
os.makedirs(CACHE, exist_ok=True)
BASE_URL = ("https://osdr.nasa.gov/geode-py/ws/studies/OSD-217/"
            "download?source=datamanager&file=")

# ---------------------------------------------------------------- inputs
module_col0 = set(pd.read_csv(f"{OUT}/redox_module_genes.csv")["agi"])
tss = pd.read_csv(f"{OUT}/gene_state_assignments.csv")
state_of = tss.set_index("agi")["pcsd_state"]
POLY = set(range(11, 16))

# Ws gene coordinates from the gene-body tables (already downloaded)
gene_coords = pd.read_excel(
    f"{WGBS}/ROOTFT.vs.ROOTGC.genes01.full.CG.xlsx") if os.path.exists(
    f"{WGBS}/ROOTFT.vs.ROOTGC.genes01.full.CG.xlsx") else None
if gene_coords is None:
    import tarfile as tfmod
    with tfmod.open(f"{WGBS}/GLDS-217_wgbs_GSE95594_ROOTFT-vs-ROOTGC_"
                    "genes01_full_tar.gz") as t:
        t.extract("ROOTFT.vs.ROOTGC.genes01.full.CG.xlsx", WGBS)
    gene_coords = pd.read_excel(
        f"{WGBS}/ROOTFT.vs.ROOTGC.genes01.full.CG.xlsx")
gene_coords.columns = ["agi", "mrna", "chrom", "start", "end", "sites",
                       "diff"]
genes_ws = gene_coords[["agi", "chrom", "start", "end"]].drop_duplicates(
    "agi")
# TAIR10 strand
g10 = cd.load_gene_coordinates(
    f"{ROOT}/data/araencode/genes_TAIR10.txt")
strand_of = g10.set_index("agi")["strand"]
genes_ws["strand"] = genes_ws["agi"].map(strand_of)
genes_ws = genes_ws.dropna(subset=["strand"])
genes_ws["tss"] = np.where(genes_ws["strand"] == "+", genes_ws["start"],
                           genes_ws["end"])
print(f"genes with Ws coords + strand: {len(genes_ws):,}")

# Ws RNA-seq DE (root and leaf)
de_ws = {}
for tissue in ["ROOT", "LEAF"]:
    fn = f"GLDS-217_rna-seq_GSE95586_{tissue}FT.vs.{tissue}GC_diff_exp.xlsx"
    dest = f"{WGBS}/{fn}"
    if not os.path.exists(dest):
        urllib.request.urlretrieve(BASE_URL + fn, dest)
    d = pd.read_excel(dest)
    d = d.rename(columns={"GeneID": "agi", "Qval": "qval", "logFC": "logfc"})
    d["qval"] = pd.to_numeric(d["qval"], errors="coerce")
    d["logfc"] = pd.to_numeric(d["logfc"], errors="coerce")
    de_ws[tissue.lower()] = d.dropna(subset=["qval", "logfc"]).drop_duplicates(
        "agi")
    print(f"Ws DE {tissue}: {len(de_ws[tissue.lower()]):,} genes, "
          f"{(de_ws[tissue.lower()]['qval'] < 0.05).sum():,} sig")

# Ws-defined module: down-DEGs in Polycomb states (root primary)
def ws_module(tissue_key):
    d = de_ws[tissue_key]
    down = d[(d["qval"] < 0.05) & (d["logfc"] < -1)]["agi"]
    return {g for g in down if state_of.get(g) in POLY}

module_ws_root = ws_module("root")
module_ws_leaf = ws_module("leaf")
print(f"Ws modules: root {len(module_ws_root)}, leaf {len(module_ws_leaf)}")

# ---------------------------------------------------------------- winavg
def parse_winavg(sample):
    """sample: e.g. ROOTFT. Returns dict context -> DataFrame(window)."""
    out = {}
    for ctx in ["CG", "CHG", "CHH"]:
        cache = f"{CACHE}/{sample}_{ctx}.csv"
        if os.path.exists(cache):
            out[ctx] = pd.read_csv(cache)
            continue
        fn = f"GLDS-217_wgbs_GSE95594_{sample}-winavg_tar.gz"
        dest = f"{WGBS}/{fn}"
        if not os.path.exists(dest):
            urllib.request.urlretrieve(BASE_URL + fn, dest)
        with tarfile.open(dest, "r:gz") as t:
            member = [m for m in t.getmembers()
                      if m.name.endswith(f"winavg_{ctx}.xlsx")][0]
            df = pd.read_excel(t.extractfile(member), header=None,
                               skiprows=1,
                               names=["chrom", "start", "end", "winavg"])
        df.to_csv(cache, index=False)
        out[ctx] = df
        print(f"  parsed {sample} {ctx}: {len(df):,} windows")
    return out

winavg = {s: parse_winavg(s) for s in
          ["ROOTFT", "ROOTGC", "LEAFFT", "LEAFGC"]}
print("winavg parsed and cached")

# ------------------------------------------------------- promoter windows
def promoter_bins(genes_sub):
    """Yield (agi, bin_index 0..19, window_start) for 2kb promoters."""
    recs = []
    for r in genes_sub.itertuples():
        if r.strand == "+":
            p0 = r.tss - 2000
        else:
            p0 = r.tss  # upstream on minus strand = higher coords;
            # window start = tss, extend to tss+2000
        for b in range(20):
            if r.strand == "+":
                ws_ = p0 + b * 100
            else:
                ws_ = r.tss + b * 100
            recs.append((r.agi, b, r.chrom, ws_))
    return recs

def promoter_meth(sample, ctx, genes_sub):
    """Overlap-weighted mean winavg over promoter bins per gene."""
    win = winavg[sample][ctx]
    vals = {}
    for chrom in genes_sub["chrom"].unique():
        w = win[win["chrom"] == chrom]
        starts = w["start"].to_numpy()
        ends = w["end"].to_numpy()
        v = w["winavg"].to_numpy()
        order = np.argsort(starts)
        starts, ends, v = starts[order], ends[order], v[order]
        g = genes_sub[genes_sub["chrom"] == chrom]
        for r in g.itertuples():
            acc, wsum = 0.0, 0.0
            for b in range(20):
                bs = (r.tss - 2000 + b * 100 if r.strand == "+"
                      else r.tss + b * 100)
                be = bs + 100
                lo = np.searchsorted(ends, bs, side="right")
                hi = np.searchsorted(starts, be, side="left")
                if hi > lo:
                    ov = (np.minimum(ends[lo:hi], be) -
                          np.maximum(starts[lo:hi], bs))
                    m = ov > 0
                    acc += float(np.dot(v[lo:hi][m], ov[m]))
                    wsum += float(ov[m].sum())
            if wsum > 0:
                vals[r.agi] = acc / wsum
    return vals

# restrict to genes present in winavg chromosomes
chroms_avail = {c for s in winavg for ctx in winavg[s]
                for c in winavg[s][ctx]["chrom"].unique()}
genes_ok = genes_ws[genes_ws["chrom"].isin(chroms_avail)]
print(f"genes on winavg chromosomes: {len(genes_ok):,}")

prom = []
for tissue, samples in [("root", ("ROOTFT", "ROOTGC")),
                        ("leaf", ("LEAFFT", "LEAFGC"))]:
    for ctx in ["CG", "CHG", "CHH"]:
        ft = promoter_meth(samples[0], ctx, genes_ok)
        gc_ = promoter_meth(samples[1], ctx, genes_ok)
        prom.extend((agi, tissue, ctx, ft[agi], gc_[agi], ft[agi] - gc_[agi])
                    for agi in set(ft) & set(gc_))
prom = pd.DataFrame(prom, columns=["agi", "tissue", "context", "flt", "gc",
                                   "change"])
print(f"promoter methylation table: {len(prom):,} rows")
prom.to_csv(f"{OUT}/osd217_promoter_methylation.csv", index=False)

# ---------------------------------------------------------------- tests
de37 = pd.read_csv(f"{OUT}/sweep/OSD-37_DE.csv", index_col=0)
ns37 = set(de37[de37["direction"] == "ns"].index)
down37 = set(de37[de37["direction"] == "down"].index)
module_col0_poly = {g for g in down37 if state_of.get(g) in POLY}
bg_poly = {g for g in ns37 if state_of.get(g) in POLY}

rows = []
for tissue in ["root", "leaf"]:
    d_ws = de_ws[tissue]
    # the DE xlsx lists only significant genes; absent = non-DE
    ws_ns = (set(prom["agi"]) - set(d_ws["agi"])
             - (module_ws_root if tissue == "root" else module_ws_leaf))
    m_ws = module_ws_root if tissue == "root" else module_ws_leaf
    for ctx in ["CG", "CHG", "CHH"]:
        sub = prom[(prom["tissue"] == tissue) & (prom["context"] == ctx)]
        for label, mod, bg in [
                ("col0_redox_module", module_col0, None),
                ("col0_polycomb_down", module_col0_poly, bg_poly),
                ("ws_polycomb_down", m_ws, ws_ns)]:
            if bg is None:
                bg = set(sub["agi"]) - mod  # all other genes
            a = sub[sub["agi"].isin(mod)]["change"]
            b = sub[sub["agi"].isin(bg)]["change"]
            if len(a) < 5:
                continue
            u, p = sps.mannwhitneyu(a, b, alternative="two-sided")
            # baseline (GC) comparison
            a0 = sub[sub["agi"].isin(mod)]["gc"]
            b0 = sub[sub["agi"].isin(bg)]["gc"]
            _, p0 = sps.mannwhitneyu(a0, b0, alternative="two-sided")
            rows.append({"tissue": tissue, "context": ctx, "set": label,
                         "n_module": len(a), "median_change_module":
                         a.median(),
                         "median_change_bg": b.median(), "p_change": p,
                         "median_gc_module": a0.median(),
                         "median_gc_bg": b0.median(), "p_baseline": p0})
tests = pd.DataFrame(rows)
ok = tests["p_change"].notna()
tests.loc[ok, "padj_change"] = cd.bh_adjust(tests.loc[ok, "p_change"].values)
okb = tests["p_baseline"].notna()
tests.loc[okb, "padj_baseline"] = cd.bh_adjust(
    tests.loc[okb, "p_baseline"].values)
tests.to_csv(f"{OUT}/osd217_promoter_methylation_tests.csv", index=False)
print(tests.to_string(index=False))

# ---------------------------------------------------------------- fig30
# bin-level table: for each gene x bin, FLT and GC methylation
# (overlap-weighted; vectorized per chromosome)
bin_rows = []
for tissue, samples in [("root", ("ROOTFT", "ROOTGC")),
                        ("leaf", ("LEAFFT", "LEAFGC"))]:
    for ctx in ["CG", "CHG", "CHH"]:
        # build per-gene bin intervals once
        bin_iv = []  # (agi, bin, chrom, start)
        for r in genes_ok.itertuples():
            for b in range(20):
                wstart = (r.tss - 2000 + b * 100 if r.strand == "+"
                          else r.tss + b * 100)
                bin_iv.append((r.agi, b, r.chrom, wstart))
        biv = pd.DataFrame(bin_iv, columns=["agi", "bin", "chrom", "bs"])
        biv["be"] = biv["bs"] + 100
        for sample, col in [(samples[0], "flt"), (samples[1], "gc")]:
            win = winavg[sample][ctx]
            out = np.full(len(biv), np.nan)
            for chrom in biv["chrom"].unique():
                w = win[win["chrom"] == chrom]
                starts = w["start"].to_numpy()
                ends = w["end"].to_numpy()
                v = w["winavg"].to_numpy()
                order = np.argsort(starts)
                starts, ends, v = starts[order], ends[order], v[order]
                m = biv["chrom"] == chrom
                bs = biv.loc[m, "bs"].to_numpy()
                be = biv.loc[m, "be"].to_numpy()
                lo = np.searchsorted(ends, bs, side="right")
                hi = np.searchsorted(starts, be, side="left")
                for k in np.flatnonzero(hi > lo):
                    a, b_ = lo[k], hi[k]
                    ov = (np.minimum(ends[a:b_], be[k]) -
                          np.maximum(starts[a:b_], bs[k]))
                    ok = ov > 0
                    if ok.any():
                        out[np.flatnonzero(m)[k]] = float(
                            np.dot(v[a:b_][ok], ov[ok]) / ov[ok].sum())
            biv[col] = out
        sub = biv.dropna(subset=["flt", "gc"])
        bin_rows.extend((r.agi, tissue, ctx, r.bin, r.flt, r.gc)
                        for r in sub.itertuples())
        print(f"  bins {tissue} {ctx}: {len(sub):,}")
bins_df = pd.DataFrame(bin_rows, columns=["agi", "tissue", "context",
                                          "bin", "flt", "gc"])
bins_df.to_csv(f"{CACHE}/promoter_bins.csv", index=False)
print(f"bin-level table: {len(bins_df):,} rows")

bins = np.arange(20)
fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
for i, tissue in enumerate(["root", "leaf"]):
    for j, ctx in enumerate(["CG", "CHG", "CHH"]):
        ax = axes[i][j]
        sub = bins_df[(bins_df["tissue"] == tissue) &
                      (bins_df["context"] == ctx)]
        for label, mod, color, ls in [
                ("redox module", module_col0, "#FD9BED", "-"),
                ("background", None, "#9e9e9e", "-")]:
            sel = sub[sub["agi"].isin(mod)] if mod else \
                sub[~sub["agi"].isin(module_col0)]
            prof_ft = sel.groupby("bin")["flt"].mean().reindex(range(20))
            prof_gc = sel.groupby("bin")["gc"].mean().reindex(range(20))
            ax.plot(bins, prof_ft, color=color, ls=ls,
                    label=f"{label} FLT")
            ax.plot(bins, prof_gc, color=color, ls="--", alpha=0.6,
                    label=f"{label} GC")
        ax.axvline(19.5, color="k", lw=0.6, ls=":")
        ax.text(19.7, 0.5, "TSS", fontsize=6, rotation=90,
                transform=ax.get_xaxis_transform())
        ax.set_title(f"{tissue} {ctx}", fontsize=9)
        if i == 0 and j == 0:
            ax.legend(fontsize=5, frameon=False)
        if i == 1:
            ax.set_xlabel("Promoter bin (2 kb upstream -> TSS)")
        if j == 0:
            ax.set_ylabel("Mean methylation")
fig.suptitle("Promoter methylation, redox module vs background (OSD-217)",
             fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIG}/fig30_promoter_metaplot.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig30_promoter_metaplot.png", dpi=150,
            bbox_inches="tight")
plt.close(fig)
print("fig30 done")

# ---------------------------------------------------------------- fig31
fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
for ax, ctx in zip(axes, ["CG", "CHG", "CHH"]):
    data, labels = [], []
    for tissue in ["root", "leaf"]:
        sub = prom[(prom["tissue"] == tissue) & (prom["context"] == ctx)]
        m_ws = module_ws_root if tissue == "root" else module_ws_leaf
        d_ws = de_ws[tissue]
        ws_ns = set(d_ws[d_ws["qval"] >= 0.05]["agi"])
        for mod, bg, lab in [
                (m_ws, ws_ns, f"{tissue} Ws-PC-down"),
                (module_col0_poly, bg_poly, f"{tissue} Col0-PC-down")]:
            a = sub[sub["agi"].isin(mod)]["change"]
            b = sub[sub["agi"].isin(bg)]["change"]
            data.append(a)
            labels.append(f"{lab}\n(n={len(a)})")
            data.append(b)
            labels.append("bg")
    bp = ax.boxplot(data, tick_labels=labels, showfliers=False)
    ax.set_title(ctx)
    ax.tick_params(axis="x", labelsize=6, rotation=45)
axes[0].set_ylabel("Promoter methylation change (FLT - GC)")
fig.suptitle("Ecotype-matched promoter methylation test", fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIG}/fig31_ecotype_matched.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig31_ecotype_matched.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("fig31 done")

# ---------------------------------------------------------------- fig32
d_root = de_ws["root"].set_index("agi")
common = d_root.index.intersection(de37.index)
lfc_col = [c for c in de37.columns
           if "log2" in c.lower() and "fc" in c.lower()][0]
sc_df = pd.DataFrame({
    "ws_logfc": d_root.loc[common, "logfc"],
    "col_logfc": de37.loc[common, lfc_col]})
fig, ax = plt.subplots(figsize=(6, 6))
pc = sc_df[sc_df.index.map(lambda g: state_of.get(g) in POLY) &
           sc_df.index.isin(down37)]
ax.scatter(sc_df["col_logfc"], sc_df["ws_logfc"], s=3, alpha=0.2,
           color="#9e9e9e", label="all genes")
ax.scatter(pc["col_logfc"], pc["ws_logfc"], s=6, alpha=0.6,
           color="#FD9BED", label="Col-0 Polycomb-down")
r = sc_df[["ws_logfc", "col_logfc"]].corr().iloc[0, 1]
ax.set_xlabel("Col-0 log2FC (OSD-37)")
ax.set_ylabel("Ws log2FC (OSD-217 root)")
ax.set_title(f"DE concordance Ws vs Col-0 (r = {r:.2f})", fontsize=10)
ax.legend(fontsize=7, frameon=False)
fig.savefig(f"{FIG}/fig32_ws_col_concordance.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig32_ws_col_concordance.png", dpi=150,
            bbox_inches="tight")
sc_df.to_csv(f"{OUT}/ws_vs_col_concordance.csv")
print(f"fig32 done (r={r:.2f})")
