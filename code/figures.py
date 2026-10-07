"""ChromDecode-At module 6: figure factory (SVG + PNG, colorblind-safe)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

# Phylo palette-derived, colorblind-safe categorical colors
PALETTE = ["#0279EE", "#FF9400", "#75A025", "#FD9BED", "#E9ED4C", "#000000",
           "#7570b3", "#e7298a"]

GROUP_COLORS = {
    "Active_transcribed": "#0279EE",
    "Polycomb_repressed": "#FD9BED",
    "Accessible_promoter": "#75A025",
    "Intergenic_quiet": "#9e9e9e",
    "Heterochromatin_TE": "#FF9400",
}


def _save(fig, out_base):
    fig.savefig(f"{out_base}.svg", bbox_inches="tight")
    fig.savefig(f"{out_base}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_state_distribution(state_counts, out_base):
    """Grouped bar chart: fraction of up/down/all genes per state group."""
    fig, ax = plt.subplots(figsize=(7, 4))
    groups = list(state_counts.index)
    x = np.arange(len(groups))
    width = 0.38
    ax.bar(x - width / 2, state_counts["up_frac"], width,
           label="up-regulated", color="#0279EE")
    ax.bar(x + width / 2, state_counts["down_frac"], width,
           label="down-regulated", color="#FF9400")
    ax.set_xticks(x, groups, rotation=20, ha="right")
    ax.set_ylabel("State-group enrichment in DE set\n(fraction in set / fraction in genome)")
    ax.axhline(1, color="k", lw=0.8, ls="--")
    ax.legend(frameon=False)
    _save(fig, out_base)


def fig_enrichment(enr, state_col, label_map=None, out_base=None,
                   title="Chromatin state enrichment"):
    """Horizontal odds-ratio chart with FDR significance markers."""
    d = enr.copy()
    if label_map is not None:
        d["label"] = d[state_col].map(label_map)
    else:
        d["label"] = d[state_col].astype(str)
    d = d.sort_values("odds_ratio")
    fig, ax = plt.subplots(figsize=(7, max(3, 0.3 * len(d))))
    colors = ["#0279EE" if o > 1 else "#9e9e9e" for o in d["odds_ratio"]]
    ax.barh(d["label"], d["odds_ratio"].clip(lower=0.01), color=colors)
    for i, (orv, q) in enumerate(zip(d["odds_ratio"], d["padj"])):
        if q < 0.05:
            star = "***" if q < 0.001 else ("**" if q < 0.01 else "*")
            ax.text(max(orv, 0.02), i, f" {star}", va="center", fontsize=9)
    ax.axvline(1, color="k", lw=0.8, ls="--")
    ax.set_xlabel("Odds ratio (Fisher exact, * FDR<0.05, ** <0.01, *** <0.001)")
    ax.set_title(title)
    _save(fig, out_base)


def fig_volcano_by_state(de, state_assign, out_base,
                         padj_col="padj", lfc_col="log2fc"):
    """Volcano plot colored by predicted/assigned state group."""
    d = de.join(state_assign, how="inner")
    d["neg_log_q"] = -np.log10(d[padj_col].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for grp, color in GROUP_COLORS.items():
        sub = d[d["state_group"] == grp]
        if len(sub):
            ax.scatter(sub[lfc_col], sub["neg_log_q"], s=3, alpha=0.35,
                       color=color, label=grp, rasterized=True)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.legend(frameon=False, markerscale=4, fontsize=8)
    _save(fig, out_base)


def fig_heatmap_features(mat, row_colors, out_base, cmap="viridis",
                         col_label="Sample"):
    """Clustered heatmap of expression (rows=genes) with state color bar."""
    g = sns.clustermap(mat, row_colors=row_colors, cmap=cmap,
                       figsize=(8, 7), yticklabels=False, xticklabels=True,
                       dendrogram_ratio=0.15, cbar_pos=(0.02, 0.8, 0.03, 0.18))
    g.ax_heatmap.set_xlabel(col_label)
    _save(g.figure, out_base)


def fig_locus_view(loci, out_base, title="Locus view"):
    """Genome-browser-style stacked tracks for one locus.

    loci: list of dicts {label, segments: DataFrame(chr,start,end,value),
    value_type: 'state' (categorical) or 'signal' (continuous)}.
    """
    n = len(loci)
    fig, axes = plt.subplots(n, 1, figsize=(9, 1.1 * n + 1), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, loc in zip(axes, loci):
        seg = loc["segments"]
        if loc["value_type"] == "state":
            for _, r in seg.iterrows():
                c = loc["colors"].get(str(r["value"]), "#cccccc")
                ax.axvspan(r["start"], r["end"], ymin=0.15, ymax=0.95,
                           color=c, alpha=0.9)
            ax.set_yticks([])
            ax.set_ylabel(loc["label"], rotation=0, ha="right", va="center",
                          fontsize=8)
        else:
            ax.fill_between(seg["pos"], 0, seg["value"], step="mid",
                            color="#0279EE", alpha=0.7)
            ax.set_ylabel(loc["label"], rotation=0, ha="right", va="center",
                          fontsize=8)
        ax.set_yticks([])
    axes[-1].set_xlabel("Chromosomal position (bp)")
    axes[0].set_title(title, fontsize=10)
    _save(fig, out_base)


def fig_confusion(cm_df, out_base):
    """Confusion matrix heatmap for the learned layer."""
    fig, ax = plt.subplots(figsize=(6, 5))
    cm_norm = cm_df.div(cm_df.sum(axis=1), axis=0)
    sns.heatmap(cm_norm, annot=True, fmt=".2f", cmap="Blues", ax=ax,
                cbar_kws={"label": "Row-normalized fraction"})
    ax.set_ylabel("True state group")
    ax.set_xlabel("Predicted state group")
    _save(fig, out_base)


def fig_som_projection(proj, out_base, hue_col="state_group"):
    """2D projection (PCA/UMAP) of genes in state-feature space."""
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for grp, color in GROUP_COLORS.items():
        sub = proj[proj[hue_col] == grp]
        if len(sub):
            ax.scatter(sub["dim1"], sub["dim2"], s=3, alpha=0.3, color=color,
                       label=grp, rasterized=True)
    ax.set_xlabel("Dimension 1")
    ax.set_ylabel("Dimension 2")
    ax.legend(frameon=False, markerscale=4, fontsize=8)
    _save(fig, out_base)
