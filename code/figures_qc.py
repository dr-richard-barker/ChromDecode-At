"""Legibility-audit regeneration of figures (see manuscript/FIGURE_AUDIT.md).

Every figure here is redrawn ONLY from archived tables:
  supplementary_tables/*.csv, data/derived/*.csv, results/sweep/*_DE.csv
No value is typed in by hand, simulated or re-sampled.  Filenames are the
originals so manuscript/supplementary references keep working.

Conventions used throughout
  * Odds ratios are drawn on a log2 axis / log colour scale centred on OR = 1.
  * OR = 0 (no DEG in the state) and OR = inf (every DEG in the state) are
    never passed to the colour map or to a bar length: they are drawn as
    explicit markers and labelled "0" / "∞".  No pseudocount is used.
  * Stars: * FDR < 0.05, ** < 0.01, *** < 0.001 (BH-FDR as stored in tables).
  * Cross-dataset heatmaps show only datasets with >= 10 DEGs in the given
    direction; the dropped ones are listed in the figure footnote.

Usage:  python code/figures_qc.py [fig10 fig12 ...]   (no args = all)
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.size"] = 9
matplotlib.rcParams["axes.titlesize"] = 10
matplotlib.rcParams["axes.labelsize"] = 9
matplotlib.rcParams["xtick.labelsize"] = 8.5
matplotlib.rcParams["ytick.labelsize"] = 8.5
matplotlib.rcParams["legend.fontsize"] = 8
matplotlib.rcParams["hatch.linewidth"] = 0.6

FIG = f"{ROOT}/figures"
ST = f"{ROOT}/supplementary_tables"
DER = f"{ROOT}/data/derived"
SWEEP = f"{ROOT}/results/sweep"

# Same categorical colours as code/figures.py (GROUP_COLORS).
GROUP_COLORS = {
    "Active_transcribed": "#0279EE",
    "Polycomb_repressed": "#FD9BED",
    "Accessible_promoter": "#75A025",
    "Intergenic_quiet": "#9e9e9e",
    "Heterochromatin_TE": "#FF9400",
}
GROUP_ORDER5 = ["Accessible_promoter", "Active_transcribed",
                "Polycomb_repressed", "Intergenic_quiet", "Heterochromatin_TE"]
STATE_ORDER12 = ["Active", "Transcription1", "Transcription2", "Bivalent",
                 "Repression", "Flanking", "Quiescent", "Heterochromatin"]
# heatmap_36state.py state_order (biological grouping)
STATE_ORDER36 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 22, 25, 26, 27, 28,
                 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24,
                 29, 30, 31, 32, 33, 34, 35, 36]
ENRICH_C, DEPLETE_C = "#0279EE", "#8a8a8a"
MIN_DEG = 10
DEG_RULE = "padj < 0.05, |log2FC| > 1"
OR_LO, OR_HI = 0.2, 5.0          # heatmap colour clip (log scale, symmetric)


# ----------------------------------------------------------------- helpers
def save(fig, name):
    fig.savefig(f"{FIG}/{name}.svg", bbox_inches="tight")
    fig.savefig(f"{FIG}/{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def stars(q):
    if q is None or not np.isfinite(q):
        return ""
    return "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else ""


def fmt_or(o):
    if not np.isfinite(o):
        return "∞" if o > 0 else "n/a"
    if o == 0:
        return "0"
    if o >= 10:
        return f"{o:.0f}"
    if o >= 1:
        return f"{o:.1f}"
    return f"{o:.2f}"


def nice(s):
    return str(s).replace("_", " ")


def sweep_summary():
    return pd.read_csv(f"{ST}/sweep_summary.csv").set_index("accession")


def n_deg(summ, acc, kind):
    return int(summ.loc[acc, "n_de_up" if kind == "up" else "n_de_down"])


def or_heatmap(ax, OR, Q, row_labels, col_labels, cax=None, fontsize=8,
               lo=OR_LO, hi=OR_HI, cbar_label=None):
    """Odds-ratio heatmap on a symmetric log2 colour scale.

    OR, Q: 2-D arrays (rows x cols).  0 -> hatched grey "0"; inf -> top colour
    with "∞"; NaN -> hatched white "n/a".
    """
    L = np.log2(hi)
    nr, nc = OR.shape
    val = np.full(OR.shape, np.nan)
    fin = np.isfinite(OR) & (OR > 0)
    val[fin] = np.log2(np.clip(OR[fin], lo, hi))
    val[np.isinf(OR) & (OR > 0)] = L
    cmap = plt.get_cmap("RdBu_r")
    norm = Normalize(-L, L)
    im = ax.imshow(np.ma.masked_invalid(val), cmap=cmap, norm=norm,
                   aspect="auto", interpolation="nearest")
    for i in range(nr):
        for j in range(nc):
            o, q = OR[i, j], Q[i, j]
            if o == 0:
                ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, facecolor="#e6e6e6",
                                       edgecolor="#9a9a9a", hatch="////", lw=0))
                txt, col = "0", "#333333"
            elif np.isnan(o):
                ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, facecolor="white",
                                       edgecolor="#bbbbbb", hatch="....", lw=0))
                txt, col = "n/a", "#333333"
            else:
                txt = fmt_or(o)
                col = "white" if abs(val[i, j]) > 0.62 * L else "#111111"
            txt += stars(q)
            ax.text(j, i, txt, ha="center", va="center", fontsize=fontsize,
                    color=col)
    ax.set_xticks(range(nc), col_labels)
    ax.set_yticks(range(nr), row_labels)
    ax.set_xticks(np.arange(-.5, nc, 1), minor=True)
    ax.set_yticks(np.arange(-.5, nr, 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    if cax is not None:
        cb = plt.colorbar(im, cax=cax)
        ticks = [t for t in [0.125, 0.2, 0.25, 0.5, 1, 2, 4, 5, 8]
                 if lo <= t <= hi]
        cb.set_ticks(np.log2(ticks), labels=[f"{t:g}" for t in ticks])
        cb.set_label(cbar_label or f"Odds ratio (log scale; colour clipped "
                     f"at {lo:g} and {hi:g})")
        cb.outline.set_visible(False)
    return im


def heat_keys(fig_or_ax, loc="lower left", bbox=(0, -0.02), ncol=3,
              has_zero=True, has_inf=True, **kw):
    handles = [Patch(facecolor="#e6e6e6", edgecolor="#9a9a9a", hatch="////",
                     label="0 = no DEG in state (OR 0)"),
               Patch(facecolor=plt.get_cmap("RdBu_r")(0.999),
                     label="∞ = every DEG in state")]
    handles = [h for h, k in zip(handles, (has_zero, has_inf)) if k]
    if not handles:
        return None
    return fig_or_ax.legend(handles=handles, loc=loc, bbox_to_anchor=bbox,
                            ncol=ncol, frameon=False, **kw)


def or_bars(ax, labels, OR, Q, K, N, fig_title=None):
    """Horizontal log2-odds bars; 0 and inf drawn as explicit markers.

    K = DEGs in state, N = DEGs with a state assignment."""
    OR = np.asarray(OR, float)
    fin = np.isfinite(OR) & (OR > 0)
    lx = np.where(fin, np.log2(np.where(fin, OR, 1)), np.nan)
    lo = min(np.nanmin(lx) if fin.any() else 0, -1) - 0.9
    hi = max(np.nanmax(lx) if fin.any() else 0, 1) + 0.9
    y = np.arange(len(OR))
    for i in range(len(OR)):
        sig = np.isfinite(Q[i]) and Q[i] < 0.05
        if OR[i] == 0:
            ax.plot(lo + 0.15, i, marker="<", ms=8, color=DEPLETE_C,
                    alpha=1 if sig else 0.45, mec="k", mew=0.5)
            ax.text(lo + 0.4, i, "OR 0", va="center", fontsize=7.5,
                    color="#333333")
        elif np.isinf(OR[i]):
            ax.barh(i, hi - 0.25, color=ENRICH_C, alpha=1 if sig else 0.4,
                    height=0.65)
            ax.plot(hi - 0.15, i, marker=">", ms=8, color=ENRICH_C, mec="k",
                    mew=0.5)
            ax.text(hi - 0.45, i, "∞", va="center", ha="right",
                    fontsize=10, color="white", fontweight="bold")
        else:
            c = ENRICH_C if OR[i] > 1 else DEPLETE_C
            ax.barh(i, lx[i], color=c, alpha=1 if sig else 0.4, height=0.65)
    # right-hand annotation column (outside axes: never collides with bars)
    for i in range(len(OR)):
        ax.text(1.02, i, f"{fmt_or(OR[i])}{stars(Q[i])}", transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=8)
        ax.text(1.20, i, f"{int(K[i])}/{int(N)}", transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=8, color="#444444")
    ax.text(1.02, len(OR) - 0.3, "OR", transform=ax.get_yaxis_transform(),
            fontsize=8, fontweight="bold", va="bottom")
    ax.text(1.20, len(OR) - 0.3, "DEGs in state", transform=ax.get_yaxis_transform(),
            fontsize=8, fontweight="bold", va="bottom")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.6, len(OR) + 0.15)
    ax.set_yticks(y, labels)
    ticks = [t for t in range(int(np.floor(lo)), int(np.ceil(hi)) + 1)]
    ax.set_xticks(ticks, [f"{2 ** t:g}" if t >= 0 else f"1/{2 ** -t}" for t in ticks])
    ax.set_xlabel("Odds ratio (log2 axis; Fisher exact test)")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    keys = [Patch(color=ENRICH_C, label="OR > 1, FDR < 0.05"),
            Patch(color=ENRICH_C, alpha=0.4, label="OR > 1, n.s."),
            Patch(color=DEPLETE_C, label="OR < 1, FDR < 0.05"),
            Patch(color=DEPLETE_C, alpha=0.4, label="OR < 1, n.s.")]
    from matplotlib.transforms import offset_copy
    ax.legend(handles=keys, loc="upper center", bbox_to_anchor=(0.5, 0),
              bbox_transform=offset_copy(ax.transAxes, fig=ax.figure, y=-36,
                                         units="points"),
              ncol=4, frameon=False, fontsize=7.5, handlelength=1.2,
              columnspacing=1.0)


def text_panel(name, title, lines):
    fig, ax = plt.subplots(figsize=(6.5, 2.2))
    ax.axis("off")
    ax.text(0.5, 0.78, title, ha="center", va="center", fontsize=11,
            fontweight="bold", transform=ax.transAxes)
    ax.text(0.5, 0.38, "\n".join(lines), ha="center", va="center",
            fontsize=9, transform=ax.transAxes, linespacing=1.5)
    save(fig, name)


# ----------------------------------------------------------- fig2 / fig12
def _enrichment_bar_fig(name, labels, OR, Q, K, N, title, note=None):
    order = np.argsort(np.where(np.isinf(OR), 1e9, OR))
    labels = [labels[i] for i in order]
    OR, Q, K = np.asarray(OR)[order], np.asarray(Q)[order], np.asarray(K)[order]
    fig, ax = plt.subplots(figsize=(6.2, 0.42 * len(OR) + 1.6))
    or_bars(ax, labels, OR, Q, K, N)
    if note:
        title = title + "\n" + note.strip()
    ax.set_title(title, loc="left", fontsize=9.5)
    save(fig, name)


def fig2():
    de = pd.read_csv(f"{ST}/OSD-37_DE_FLT_vs_GC.csv")
    ndeg = de["direction"].value_counts()
    for kind in ["up", "down"]:
        e = pd.read_csv(f"{ST}/enrichment_{kind}.csv")
        N = int(e["n_in_set"].sum())
        _enrichment_bar_fig(
            f"fig2_enrichment_{kind}", [nice(s) for s in e["st"]],
            e["odds_ratio"].values, e["padj"].values, e["n_in_set"].values, N,
            f"OSD-37 (v1, GeneLab VST), {kind}-regulated DEGs: "
            f"n = {ndeg.get(kind, 0)} ({N} with a state)\n"
            f"FLT vs GC, {DEG_RULE}; tables S1, S{2 if kind == 'up' else 3}")


def fig12():
    t = pd.read_csv(f"{ST}/sweep_enrichment5_corrected.csv")
    summ = sweep_summary()
    for acc in ["OSD-37", "OSD-120", "OSD-251", "OSD-314", "OSD-346"]:
        contrast = str(summ.loc[acc, "contrast"]).replace("parsed_condition: ", "")
        for kind in ["up", "down"]:
            name = f"fig12_{acc}_enrichment_{kind}"
            n = n_deg(summ, acc, kind)
            e = t[(t["accession"] == acc) & (t["kind"] == kind)]
            if n <= 1:
                where = ""
                if len(e) and (e["n_set"] > 0).any():
                    where = (" (in " + ", ".join(nice(s) for s in
                             e.loc[e["n_set"] > 0, "state_group"]) + ")")
                text_panel(name, f"{acc}: {kind}-regulated DEGs",
                           [f"n = {n} DEG{'s' if n != 1 else ''}{where} at {DEG_RULE}",
                            "State enrichment is not estimable from ≤ 1 gene;",
                            "panel intentionally left without a chart "
                            + ("(table S51 rows kept for completeness)."
                               if len(e) else "(no rows in table S51).")])
                continue
            e = e.set_index("state_group").reindex(GROUP_ORDER5)
            N = int(e["n_set"].sum())
            note = None
            if n < MIN_DEG:
                note = (f"Only {n} DEGs (< {MIN_DEG}): low power"
                        + ("; no state reaches FDR < 0.05"
                           if not (e["padj"] < 0.05).any() else ""))
            _enrichment_bar_fig(
                name, [nice(s) for s in e.index], e["odds_ratio"].values,
                e["padj"].values, e["n_set"].values, N,
                f"{acc}, {kind}-regulated DEGs: n = {n} ({N} with a state)\n"
                f"contrast as parsed: {contrast}; {DEG_RULE}; table S51",
                note=("\n\n" + note) if note else None)


# ----------------------------------------------------------------- fig9
def fig9():
    t = pd.read_csv(f"{ST}/learned_v2_transfer.csv", index_col=0)
    names = {"AraENCODE_seedling(held-out study)": "AraENCODE seedling\n(held-out study)",
             "AraENCODE_leaf": "AraENCODE leaf",
             "AraENCODE_root": "AraENCODE root",
             "OSD-37_spaceflight_seedling": "OSD-37 spaceflight\nseedling"}
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    x = np.arange(len(t))
    for i, (k, r) in enumerate(t.iterrows()):
        ax.plot([i - 0.32, i + 0.32], [r["baseline"]] * 2, color="k", lw=1.6,
                ls="--")
        ax.plot(i, r["accuracy"], "o", ms=9, color=ENRICH_C, mec="k", mew=0.5)
        d = r["accuracy"] - r["baseline"]
        ax.annotate(f"{r['accuracy']:.3f}\n(Δ {d:+.3f})", (i, r["accuracy"]),
                    xytext=(0, -12), textcoords="offset points", ha="center",
                    va="top", fontsize=8)
    ax.set_xticks(x, [f"{names.get(k, k)}\nn = {int(r['n']):,}"
                      for k, r in t.iterrows()])
    ax.set_xlim(-0.6, len(t) - 0.4)
    lo = min(t["accuracy"].min(), t["baseline"].min())
    hi = max(t["accuracy"].max(), t["baseline"].max())
    ax.set_ylim(lo - 0.03, hi + 0.012)
    ax.set_ylabel("Accuracy (fraction of genes,\n5 state groups)")
    ax.set_title("v2 learned layer: transfer accuracy vs majority-class baseline\n"
                 "No target exceeds the baseline (table S9)", fontsize=9.5)
    h = [Line2D([], [], marker="o", ls="", color=ENRICH_C, mec="k", ms=8,
                label="transfer accuracy"),
         Line2D([], [], color="k", ls="--", lw=1.6,
                label=f"majority-class baseline (S9: {t['baseline'].iloc[0]:.3f})")]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.3),
              ncol=2, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save(fig, "fig9_v2_transfer")


# ------------------------------------------------------------- fig6 / fig8
def _confusion(name, csv, title):
    cm = pd.read_csv(f"{ST}/{csv}", index_col=0)
    cm.index = [nice(i.replace("true_", "")) for i in cm.index]
    cm.columns = [nice(c.replace("pred_", "")) for c in cm.columns]
    norm = cm.div(cm.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    im = ax.imshow(norm.values, cmap="Blues", vmin=0, vmax=1)
    for i in range(norm.shape[0]):
        for j in range(norm.shape[1]):
            v = norm.values[i, j]
            ax.text(j, i, f"{v:.2f}\n({cm.values[i, j]:,})", ha="center",
                    va="center", fontsize=7.5,
                    color="white" if v > 0.55 else "#111111")
    lab = [s.replace(" ", "\n") for s in cm.columns]
    ax.set_xticks(range(len(lab)), [s.replace("\n", " ") for s in lab], fontsize=8, rotation=30, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(cm.index)), [s.replace(" ", "\n") for s in cm.index],
                  fontsize=8)
    ax.set_xlabel("Predicted state group")
    ax.set_ylabel("True state group")
    ax.set_title(title, fontsize=9.5)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Row-normalised fraction (0–1; count in brackets)")
    save(fig, name)


def fig6():
    _confusion("fig6_confusion_matrix", "learned_confusion_matrix.csv",
               "v1 classifier: cross-validated confusion matrix (table S5)")


def fig8():
    _confusion("fig8_v2_confusion_matrix", "learned_v2_confusion.csv",
               "v2 classifier: cross-validated confusion matrix (table S8)")


# ------------------------------------------------- fig10 / fig13 / fig14
def _dataset_heatmap(name, long, statecol, states, kind, title, table_ref,
                     col_label_fn, figw):
    summ = sweep_summary()
    d = long[long["kind"] == kind]
    accs = sorted(d["accession"].unique(), key=lambda a: -n_deg(summ, a, kind))
    keep = [a for a in accs if n_deg(summ, a, kind) >= MIN_DEG]
    ok = summ[summ["status"] == "ok"].index
    dropped = [f"{a} (n={n_deg(summ, a, kind)})" for a in sorted(ok)
               if a not in keep]
    OR = np.full((len(keep), len(states)), np.nan)
    Q = OR.copy()
    for i, a in enumerate(keep):
        s = d[d["accession"] == a].set_index(statecol)
        for j, st in enumerate(states):
            if st in s.index:
                OR[i, j], Q[i, j] = s.loc[st, "odds_ratio"], s.loc[st, "padj"]
    from matplotlib.transforms import offset_copy
    fig, ax = plt.subplots(figsize=(figw, 0.62 * len(keep) + 0.9))
    cax = ax.inset_axes([1.03, 0.0, 0.025, 1.0])
    or_heatmap(ax, OR, Q, [f"{a}\n(n = {n_deg(summ, a, kind)})" for a in keep],
               [col_label_fn(s) for s in states], cax=cax, fontsize=8.5,
               cbar_label="Odds ratio (log scale;\ncolour clipped at 0.2 and 5)")
    ax.set_xlabel("Chromatin state" + (" group" if statecol == "state_group" else ""))
    ax.set_title(title, fontsize=9.5)
    foot = (f"Rows: datasets with \u2265 {MIN_DEG} {kind}-regulated DEGs ({DEG_RULE}). "
            f"Not shown: {', '.join(dropped) if dropped else 'none'}. "
            f"* FDR < 0.05, ** < 0.01, *** < 0.001 (BH). Source: {table_ref}.")
    tr = offset_copy(ax.transAxes, fig=fig, y=-62, units="points")
    heat_keys(ax, loc="upper left", bbox=(0, 0), ncol=2, fontsize=7.5,
              has_zero=bool((OR == 0).any()), has_inf=bool(np.isinf(OR).any()),
              bbox_transform=tr, borderaxespad=0, handlelength=1.6)
    tr2 = offset_copy(ax.transAxes, fig=fig, y=-80, units="points")
    ax.text(0, 0, textwrap.fill(foot, int(figw * 13)), transform=tr2,
            fontsize=7.5, va="top")
    save(fig, name)


def fig10():
    t = pd.read_csv(f"{ST}/sweep_enrichment5_corrected.csv")
    _dataset_heatmap("fig10_sweep_heatmap_up", t, "state_group", GROUP_ORDER5,
                     "up", "5-group state enrichment of up-regulated DEGs "
                     "(corrected sweep)", "sweep_enrichment5_corrected.csv (S51)",
                     lambda s: nice(s).replace(" ", "\n"), 7.0)


def fig13_14():
    t = pd.read_csv(f"{ST}/sweep_enrichment12_long.csv")
    for kind, name in [("up", "fig13_sweep_heatmap_12state_up"),
                       ("down", "fig14_sweep_heatmap_12state_down")]:
        _dataset_heatmap(name, t, "state", STATE_ORDER12, kind,
                         f"AraENCODE 12-state model: enrichment of {kind}-regulated DEGs",
                         "sweep_enrichment12_long.csv (S11)",
                         lambda s: s.replace("Transcription", "Transcr.\n")
                         .replace("Heterochromatin", "Hetero-\nchromatin"), 7.8)


# --------------------------------------------------------- fig15 / fig16
def pcsd_labels():
    """PCSD state mark strings.  state_descriptions.tsv is not in the archive,
    so the 'S<n>: marks' strings are extracted from the text of the ORIGINAL
    fig15 SVG as committed in git (HEAD), falling back to the working copy.
    The original code truncated each mark list at 38 characters; a truncated
    final item is dropped and replaced by an ellipsis."""
    import re
    import subprocess
    rel = "figures/fig15_sweep_heatmap_36state_up.svg"
    try:
        svg = subprocess.run(["git", "-C", ROOT, "show", f"HEAD:{rel}"],
                             capture_output=True, text=True, check=True).stdout
    except Exception:
        svg = open(f"{ROOT}/{rel}").read()
    raw = {int(m.group(1)): m.group(2).strip()
           for m in re.finditer(r"S(\d+): ([^<]*)", svg)}
    if len(raw) != 36:
        raise RuntimeError(f"expected 36 PCSD state labels, found {len(raw)}")
    lab = {}
    for s, marks in raw.items():
        items = [x.strip() for x in marks.split(",") if x.strip()]
        trunc = len(marks) >= 38
        if trunc and len(items) > 1:
            items = items[:-1]
        short = ", ".join(items[:3])
        if trunc or len(items) > 3:
            short += ", \u2026"
        lab[s] = f"S{s}  {short}"
    return lab


def fig15_16():
    lab = pcsd_labels()
    t = pd.read_csv(f"{ST}/sweep_enrichment36_long.csv")
    t["state"] = t["state"].astype(int)
    gsa = pd.read_csv(f"{DER}/gene_state_assignments.csv").dropna(subset=["pcsd_state"])
    grp = gsa.groupby(gsa["pcsd_state"].astype(int))["state_group"].agg(
        lambda x: x.mode()[0])
    summ = sweep_summary()
    for kind, name in [("up", "fig15_sweep_heatmap_36state_up"),
                       ("down", "fig16_sweep_heatmap_36state_down")]:
        d = t[t["kind"] == kind]
        accs = sorted(d["accession"].unique(), key=lambda a: -n_deg(summ, a, kind))
        keep = [a for a in accs if n_deg(summ, a, kind) >= MIN_DEG]
        ok = summ[summ["status"] == "ok"].index
        dropped = [f"{a} (n={n_deg(summ, a, kind)})" for a in sorted(ok)
                   if a not in keep]
        states = STATE_ORDER36
        OR = np.full((len(states), len(keep)), np.nan)
        Q = OR.copy()
        for j, a in enumerate(keep):
            s = d[d["accession"] == a].set_index("state")
            for i, st in enumerate(states):
                if st in s.index:
                    OR[i, j], Q[i, j] = s.loc[st, "odds_ratio"], s.loc[st, "padj"]
        w = 1.25 * len(keep)
        fig = plt.figure(figsize=(3.6 + w + 1.2, 10.6))
        H = 0.80
        x0 = 3.5 / (3.6 + w + 1.2)
        ww = w / (3.6 + w + 1.2)
        ax = fig.add_axes([x0, 0.13, ww, H])
        sax = fig.add_axes([x0 - 0.022, 0.13, 0.018, H], sharey=ax)
        cax = fig.add_axes([x0 + ww + 0.04, 0.13 + H * 0.3, 0.03, H * 0.4])
        or_heatmap(ax, OR, Q, [lab[s] for s in states],
                   [f"{a}\n(n = {n_deg(summ, a, kind)})" for a in keep],
                   cax=cax, fontsize=8.5,
                   cbar_label="Odds ratio (log scale;\ncolour clipped at 0.2 and 5)")
        ax.xaxis.tick_top()
        ax.tick_params(axis="y", labelleft=False)
        # group strip + separators
        gcol = [GROUP_COLORS[grp[s]] for s in states]
        for i, c in enumerate(gcol):
            sax.add_patch(Rectangle((0, i - .5), 1, 1, color=c, lw=0))
        sax.set_xlim(0, 1)
        sax.set_xticks([])
        sax.set_yticks(range(len(states)), [lab[s] for s in states], fontsize=8)
        sax.tick_params(which="both", length=0)
        ax.tick_params(which="both", length=0)
        for sp in sax.spines.values():
            sp.set_visible(False)
        for i in range(1, len(states)):
            if grp[states[i]] != grp[states[i - 1]]:
                ax.axhline(i - 0.5, color="k", lw=1.2)
                sax.axhline(i - 0.5, color="k", lw=1.2)
        ax.set_title(f"PCSD 36-state enrichment of {kind}-regulated DEGs",
                     fontsize=10, pad=34)
        present = [g for g in GROUP_ORDER5 if g in set(grp[states])]
        gh = [Patch(color=GROUP_COLORS[g], label=nice(g)) for g in present]
        zh = [Patch(facecolor="#e6e6e6", edgecolor="#9a9a9a", hatch="////",
                    label="0 = no DEG in state"),
              Patch(facecolor=plt.get_cmap("RdBu_r")(0.999),
                    label="∞ = every DEG in state")]
        fig.legend(handles=gh, title="State group (left strip)", loc="lower left",
                   bbox_to_anchor=(0.02, 0.045), ncol=3, frameon=False,
                   fontsize=7.5, title_fontsize=8)
        zh = [h for h, k in zip(zh, ((OR == 0).any(), np.isinf(OR).any())) if k]
        fig.legend(handles=zh, loc="lower left", bbox_to_anchor=(0.02, 0.012),
                   ncol=2, frameon=False, fontsize=7.5)
        foot = (f"Columns: datasets with ≥ {MIN_DEG} {kind}-regulated DEGs "
                f"({DEG_RULE}); not shown: {', '.join(dropped) or 'none'}. "
                "* FDR < 0.05, ** < 0.01, *** < 0.001 (BH). Mark lists shortened "
                "(…); full lists in the PCSD state table. Source: "
                "sweep_enrichment36_long.csv (S12).")
        fig.text(0.02, 0.0, textwrap.fill(foot, 95), fontsize=7.2, va="top")
        save(fig, name)


# ----------------------------------------------------------------- fig19
def fig19():
    g = pd.read_csv(f"{ST}/osd314_threshold_grid.csv")
    cols = ["S11", "S12", "S13", "S15"]
    V = g[cols].values.astype(float)
    hi = 2.0 ** np.ceil(np.log2(np.nanmax(V)))
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    lv = np.log2(np.clip(V, 1, hi))
    im = ax.imshow(lv, cmap="Reds", vmin=0, vmax=np.log2(hi), aspect="auto")
    for i in range(V.shape[0]):
        for j in range(V.shape[1]):
            ax.text(j, i, f"{V[i, j]:.1f}", ha="center", va="center", fontsize=8.5,
                    color="white" if lv[i, j] > 0.6 * np.log2(hi) else "#111111")
    ax.set_xticks(range(len(cols)), [f"{c}" for c in cols])
    ax.set_yticks(range(len(g)), [f"padj < {r.padj:g}, |log2FC| > {r.lfc:g}"
                                  f"   (n = {int(r.n_up)})" for r in g.itertuples()])
    ax.set_xlabel("PCSD Polycomb state")
    ax.set_ylabel("DE threshold (n up-regulated genes)")
    for i in (3, 6):
        ax.axhline(i - 0.5, color="white", lw=2)
    cb = plt.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    ticks = [t for t in [1, 2, 4, 8, 16] if t <= hi]
    cb.set_ticks(np.log2(ticks), labels=[f"{t:g}" for t in ticks])
    cb.set_label("Odds ratio, up-regulated genes\n(log scale from OR = 1)")
    ax.set_title("OSD-314 Polycomb-state odds across DE thresholds (table S15)\n"
                 "S15 stores odds only (no FDR), so no significance marks",
                 fontsize=9.5)
    save(fig, "fig19_threshold_grid")


# ----------------------------------------------------------------- fig20
def fig20():
    """Deterministic running sum, recomputed exactly as diagnose_osd314.py
    (test 4) from results/sweep/OSD-314_DE.csv; checked against S16."""
    gsea = pd.read_csv(f"{ST}/osd314_gsea.csv").iloc[0]
    gsa = pd.read_csv(f"{DER}/gene_state_assignments.csv").dropna(subset=["pcsd_state"])
    state_of = gsa.set_index("agi")["pcsd_state"].astype(int)
    poly = set(state_of[state_of.isin([11, 12, 13, 15])].index)
    de = pd.read_csv(f"{SWEEP}/OSD-314_DE.csv", index_col=0)
    ranked = de["log2fc"].sort_values(ascending=False)
    hit = ranked.index.isin(poly)
    nh, nm = hit.sum(), len(hit) - hit.sum()
    run = np.cumsum(np.where(hit, 1 / nh, -1 / nm))
    if nh != int(gsea["n_poly_genes"]) or not np.isclose(run.max(), gsea["max_running_sum"]):
        raise RuntimeError(f"fig20 recompute mismatch: n={nh} vs {gsea['n_poly_genes']}, "
                           f"max={run.max()} vs {gsea['max_running_sum']}")
    n_perm = 1000  # diagnose_osd314.py test 4
    p = gsea["perm_p"]
    ptxt = (f"permutation p < {1 / n_perm:g} (0 of {n_perm:,} permutations "
            "≥ observed)") if p == 0 else f"permutation p = {p:.3g}"
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(6.6, 4.4), sharex=True,
                                  gridspec_kw={"height_ratios": [4, 0.7],
                                               "hspace": 0.05})
    ax.plot(np.arange(len(run)), run, color=ENRICH_C, lw=1.5)
    ax.axhline(0, color="k", lw=0.5)
    im = int(np.argmax(run))
    ax.plot(im, run[im], "o", color="k", ms=4)
    ax.annotate(f"max running sum = {run[im]:.3f}", (im, run[im]),
                xytext=(10, 0), textcoords="offset points", va="center", fontsize=8)
    ax.set_ylabel("Running enrichment sum")
    ax.set_title(f"OSD-314: Polycomb-state genes (S11/12/13/15, n = {nh:,})\n{ptxt}",
                 fontsize=9.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    rax.vlines(np.where(hit)[0], 0, 1, color="#c2187a", lw=0.25, alpha=0.6)
    rax.set_yticks([0.5], ["Polycomb\ngenes"], fontsize=7.5)
    rax.set_ylim(0, 1)
    rax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        rax.spines[s].set_visible(False)
    rax.set_xlabel(f"Gene rank by log2 fold change (1 = most up-regulated; "
                   f"{len(run):,} genes)")
    save(fig, "fig20_running_enrichment")


# ----------------------------------------------------------------- fig21
def fig21():
    h = pd.read_csv(f"{ST}/osd314_up_polycomb_genes.csv").sort_values("log2fc")
    marks = {"S11": "o", "S12": "s", "S13": "^", "S15": "D"}
    fig, ax = plt.subplots(figsize=(5.4, 0.2 * len(h) + 1.6))
    y = np.arange(len(h))
    for st, mk in marks.items():
        m = (h["state"] == st).values
        if m.any():
            ax.scatter(h["log2fc"].values[m], y[m], marker=mk, s=34,
                       color="#444444", label=f"{st} (n = {m.sum()})", zorder=3)
    ax.hlines(y, 0, h["log2fc"].values, color="#bbbbbb", lw=0.8, zorder=1)
    ax.set_yticks(y, [f"{r.agi} ({r.state})" for r in h.itertuples()], fontsize=7.5)
    ax.set_ylim(-0.7, len(h) - 0.3)
    ax.set_xlim(0, h["log2fc"].max() * 1.08)
    ax.set_xlabel("log2 fold change (OSD-314 up-regulated)")
    ax.set_ylabel("Gene (PCSD state)")
    ax.set_title(f"OSD-314 up-regulated genes in Polycomb states (n = {len(h)}; S17)",
                 fontsize=9.5)
    ax.legend(title="PCSD state", loc="lower right", frameon=False, fontsize=7.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save(fig, "fig21_polycomb_up_genes")


# ----------------------------------------------------------------- fig22
def fig22():
    m = pd.read_csv(f"{ST}/motif_enrichment.csv")
    u = m[m["comparison"] == "target_vs_universe"].copy()
    u["tf"] = u["name"].str.split("\t").str[-1].str.strip()
    sig = u[(u["padj"] < 0.05) & (u["n_target_hits"] > 0)]
    enr = sig[sig["odds_ratio"] > 1].sort_values("odds_ratio", ascending=False).head(12)
    dep = sig[sig["odds_ratio"] < 1].sort_values("odds_ratio").head(12)
    top = pd.concat([dep.iloc[::-1], enr.iloc[::-1]])
    # target-set size: an infinite OR means every target promoter has the hit
    n_target = int(u.loc[np.isinf(u["odds_ratio"]), "n_target_hits"].max()) \
        if np.isinf(u["odds_ratio"]).any() else None
    OR, Q = top["odds_ratio"].values, top["padj"].values
    fig, ax = plt.subplots(figsize=(6.6, 0.27 * len(top) + 1.8))
    N = n_target if n_target else np.nan
    or_bars(ax, [f"{r.tf} ({r.motif})" for r in top.itertuples()], OR, Q,
            top["n_target_hits"].values, N)
    ax.texts[-1].set_text("target hits")
    ax.tick_params(axis="y", labelsize=7.5)
    ax.axhline(len(dep) - 0.5, color="#cccccc", lw=0.8, ls=":")
    ax.set_title("JASPAR 2024 motifs: Polycomb-down promoters vs expressed genome\n"
                 f"top {len(enr)} enriched and {len(dep)} depleted at FDR < 0.05 "
                 f"(n target = {n_target}; table S18)", fontsize=9.5, loc="left")
    save(fig, "fig22_motif_enrichment")
    return n_target


# ----------------------------------------------------------------- fig23
def fig23():
    g = pd.read_csv(f"{ST}/go_enrichment.csv")
    d = g[(g["comparison"] == "target_vs_polycomb_nonDE") & (g["padj"] < 0.1)].head(15)
    d = d.iloc[::-1]
    ns = {"biological_process": "BP", "molecular_function": "MF",
          "cellular_component": "CC"}
    fig, ax = plt.subplots(figsize=(6.6, 0.36 * len(d) + 1.9))
    lq = -np.log10(d["padj"].values)
    sc = ax.scatter(d["fold"], range(len(d)), s=8 * d["n_in_set"], c=lq,
                    cmap="viridis", vmin=1, vmax=np.ceil(lq.max()),
                    edgecolor="k", lw=0.4, zorder=3)
    ax.set_yticks(range(len(d)), [textwrap.fill(f"{r['name']} [{ns[r['namespace']]}]", 44)
                                  for _, r in d.iterrows()], fontsize=7.5)
    ax.axvline(1, color="k", lw=0.8, ls="--")
    ax.set_xlim(0, d["fold"].max() * 1.15)
    ax.set_ylim(-0.8, len(d) - 0.2)
    ax.set_xlabel("Fold enrichment (Polycomb-down vs Polycomb non-DE)")
    ax.grid(axis="y", color="#eeeeee", zorder=0)
    cb = plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("−log10 FDR")
    sizes = [5, 10, 25, 50]
    hh = [plt.scatter([], [], s=8 * s, color="#bbbbbb", edgecolor="k", lw=0.4)
          for s in sizes]
    ax.legend(hh, [str(s) for s in sizes], title="Genes in term", loc="upper center",
              bbox_to_anchor=(0.5, -0.12 - 0.3 / len(d)), ncol=4, frameon=False,
              fontsize=7.5, title_fontsize=8)
    ax.set_title("GO over-representation in Polycomb-down genes (FDR < 0.1; table S19)\n"
                 "BP biological process, MF molecular function, CC cellular component",
                 fontsize=9.5, loc="left")
    save(fig, "fig23_go_dotplot")


# ----------------------------------------------------------------- fig33
def fig33():
    t = pd.read_csv(f"{ST}/sog1_promoter_tests.csv")
    sets = [("redox_module", "Redox module"), ("polycomb_down", "Polycomb-down"),
            ("nonpolycomb_down", "non-Polycomb-down"), ("background", "Polycomb non-DE")]
    peaks = [("20min", "20 min", "#FD9BED"), ("1h", "1 h", "#0279EE"),
             ("union", "20 min ∪ 1 h", "#75A025"),
             ("vsWT_union", "SOG1 vs WT, 20 min ∪ 1 h", "#7a7a7a")]
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    w = 0.2
    ymax = 0
    for k, (pk, plab, col) in enumerate(peaks):
        s = t[t["peaks"] == pk].set_index("gene_set")
        for i, (gs, _) in enumerate(sets):
            r = s.loc[gs]
            pct = 100 * r["n_with_peak"] / r["n_genes"]
            ymax = max(ymax, pct)
            x = i + (k - 1.5) * w
            ax.bar(x, pct, w * 0.92, color=col, label=plab if i == 0 else None)
            ax.text(x, pct, f" {int(r['n_with_peak'])}/{int(r['n_genes'])}{stars(r['padj'])}",
                    rotation=90, ha="center", va="bottom", fontsize=7)
    ax.set_ylim(0, max(ymax * 1.9, 0.6))
    ax.set_xticks(range(len(sets)), [f"{lab}\n(n = {int(t[t.gene_set == gs].n_genes.iloc[0]):,})"
                                     for gs, lab in sets])
    ax.set_ylabel("Promoters with a SOG1 peak (%)\n(2 kb upstream of TSS)")
    ax.set_title("SOG1 ChIP-seq (OSD-496): peak overlap per gene set (table S26)\n"
                 "Labels: genes with peak / genes in set; no test reaches FDR < 0.05",
                 fontsize=9.5)
    ax.legend(title="Peak set", frameon=False, loc="upper left", fontsize=7.5,
              title_fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save(fig, "fig33_sog1_enrichment")


ALL = {"fig2": fig2, "fig6": fig6, "fig8": fig8, "fig9": fig9, "fig10": fig10,
       "fig12": fig12, "fig13_14": fig13_14, "fig15_16": fig15_16,
       "fig19": fig19, "fig20": fig20, "fig21": fig21, "fig22": fig22,
       "fig23": fig23, "fig33": fig33}

if __name__ == "__main__":
    want = sys.argv[1:] or list(ALL)
    for k in want:
        ALL[k]()
