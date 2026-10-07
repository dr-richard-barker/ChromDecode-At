"""OSD-314 Polycomb-enrichment diagnostics: real signal or small-N noise?

Five tests from cached data:
  1. Empirical null calibration (10,000 random 47-gene draws)
  2. OSD-37 downsampling control (529 -> 47 genes, 1,000 draws)
  3. Threshold-sensitivity grid (padj x |log2fc|)
  4. Threshold-free GSEA-style running enrichment
  5. Gene-identity check (which up genes sit in S11/S15)
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
rng = np.random.RandomState(42)

tss = pd.read_csv(f"{OUT}/gene_state_assignments.csv")
st = tss.dropna(subset=["pcsd_state"]).rename(columns={"pcsd_state": "st"})[
    ["agi", "st"]]
st["st"] = st["st"].astype(int)
state_of = st.set_index("agi")["st"]
bg = st["agi"].unique()

de314 = pd.read_csv(f"{SWEEP}/OSD-314_DE.csv", index_col=0)
de37 = pd.read_csv(f"{SWEEP}/OSD-37_DE.csv", index_col=0)
up314 = de314[de314["direction"] == "up"].index
POLY = [11, 12, 13, 15]

results = []


def odds_for(gene_set, state):
    genes_in = set(state_of[state_of == state].index)
    a = len(set(gene_set) & genes_in)
    b = len(gene_set) - a
    c = len(genes_in) - a
    d = len(bg) - len(gene_set) - c
    from scipy import stats as sps
    o, p = sps.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return o, p, a


# ---------------------------------------------------------------- test 1
print("[1] Empirical null calibration (10,000 draws of n=47) ...")
n_perm = 10_000
null = {s: np.zeros(n_perm) for s in POLY}
for i in range(n_perm):
    draw = rng.choice(bg, size=len(up314), replace=False)
    for s in POLY:
        o, _, _ = odds_for(draw, s)
        null[s][i] = o
emp = {}
for s in POLY:
    obs_o, obs_p, a = odds_for(up314, s)
    emp_p = (null[s] >= obs_o).mean()
    emp[f"S{s}"] = {"observed_odds": obs_o, "fisher_p": obs_p,
                    "empirical_p": emp_p, "n_genes_in_state": a,
                    "null_mean": null[s].mean(),
                    "null_p99": np.quantile(null[s], 0.99)}
    print(f"  S{s}: observed odds={obs_o:.1f} (in {a} up genes), "
          f"null mean={null[s].mean():.2f}, null p99="
          f"{np.quantile(null[s], 0.99):.2f}, empirical p={emp_p:.4f}")
pd.DataFrame(emp).T.to_csv(f"{OUT}/osd314_null_calibration.csv")

# ---------------------------------------------------------------- test 2
print("[2] OSD-37 downsampling control (529 -> 47, 1,000 draws) ...")
up37 = de37[de37["direction"] == "up"].index
ds = {s: np.zeros(1000) for s in POLY}
for i in range(1000):
    sub = rng.choice(up37, size=len(up314), replace=False)
    for s in POLY:
        o, _, _ = odds_for(sub, s)
        ds[s][i] = o
ds_res = {}
for s in POLY:
    obs_o = emp[f"S{s}"]["observed_odds"]
    ds_res[f"S{s}"] = {"osd37_downsampled_mean": ds[s].mean(),
                       "osd37_downsampled_p99": np.quantile(ds[s], 0.99),
                       "osd314_observed": obs_o,
                       "frac_draws_ge_observed": (ds[s] >= obs_o).mean()}
    print(f"  S{s}: OSD-37@n=47 mean={ds[s].mean():.2f}, "
          f"p99={np.quantile(ds[s], 0.99):.2f}, "
          f"P(draw >= OSD-314 observed)={(ds[s] >= obs_o).mean():.4f}")
pd.DataFrame(ds_res).T.to_csv(f"{OUT}/osd314_downsampling_control.csv")

# ---------------------------------------------------------------- test 3
print("[3] Threshold-sensitivity grid ...")
counts = pd.read_csv(f"{SWEEP}/OSD-314_counts.csv", index_col=0)
counts.index = [cd.clean_agi(i) or i for i in counts.index]
counts = counts[~counts.index.str.startswith(("GLDS", "SRR"))]
cond = np.array([1 if "_0G" in s.upper() or s.upper().endswith(
    ("_0G", "_0G_")) or "0G" in s.upper().split("_") else 0
    for s in counts.columns])
# robust: reuse the sweep's contrast logic
fac_cond = None
meta = pd.read_csv(f"{SWEEP}/OSD-314_metadata.csv", dtype=str)
import run_osdr_sweep as R
fac = R.sample_factors(meta, list(counts.columns))
if fac is not None and "parsed_condition" in fac.columns:
    lv = fac["parsed_condition"].reindex(counts.columns)
    ctrl = lv.mode()[0]
    fac_cond = (lv != ctrl).astype(float).values
cond = fac_cond if fac_cond is not None else cond
mat = np.log2(counts.clip(lower=1))
grid_rows = []
for padj_t in [0.01, 0.05, 0.1]:
    for lfc_t in [0.5, 1.0, 2.0]:
        de = cd.differential_expression_any(mat, cond,
                                            padj_threshold=padj_t,
                                            lfc_threshold=lfc_t)
        up = de[de["direction"] == "up"].index
        row = {"padj": padj_t, "lfc": lfc_t, "n_up": len(up)}
        for s in POLY:
            if len(up) == 0:
                row[f"S{s}"] = np.nan
            else:
                o, _, _ = odds_for(up, s)
                row[f"S{s}"] = o
        grid_rows.append(row)
grid = pd.DataFrame(grid_rows)
grid.to_csv(f"{OUT}/osd314_threshold_grid.csv", index=False)
print(grid.to_string(index=False))

# ---------------------------------------------------------------- test 4
print("[4] Threshold-free running enrichment (GSEA-style) ...")
poly_genes = set(state_of[state_of.isin(POLY)].index)
lfc = de314["log2fc"]
ranked = lfc.sort_values(ascending=False)
hit = ranked.index.isin(poly_genes)
n_hits = hit.sum()
n_miss = len(hit) - n_hits
inc = np.where(hit, 1 / n_hits, -1 / n_miss)
run = np.cumsum(inc)
nes_obs = run.max()

perm_nes = np.zeros(1000)
genes_all = ranked.index.values
for i in range(1000):
    phit = rng.permutation(hit)
    pinc = np.where(phit, 1 / n_hits, -1 / n_miss)
    prun = np.cumsum(pinc)
    perm_nes[i] = prun.max() if abs(prun).max() == prun.max() else prun.min()
# standardize
pos = perm_nes[perm_nes > 0]
nes_std = (nes_obs - pos.mean()) / pos.std() if len(pos) else np.nan
perm_p = (perm_nes >= nes_obs).mean()
print(f"  Polycomb set: {n_hits} genes; observed max running sum="
      f"{nes_obs:.3f}; permutation p={perm_p:.4f}")
pd.DataFrame({"n_poly_genes": [n_hits], "max_running_sum": [nes_obs],
              "perm_p": [perm_p]}).to_csv(
    f"{OUT}/osd314_gsea.csv", index=False)

# ---------------------------------------------------------------- test 5
print("[5] Gene-identity check ...")
genes_fn = cd.load_gene_coordinates(f"{DATA}/araencode/genes_TAIR10.txt")
desc_map = {}
try:
    fn = pd.read_csv(f"{DATA}/araencode/TAIR10_functional_descriptions_20130831.txt",
                     sep="\t", header=None, usecols=[0, 1],
                     names=["agi", "desc"], on_bad_lines="skip")
    desc_map = fn.set_index("agi")["desc"]
except Exception as e:
    print(f"  (no functional descriptions: {e})")
hits = []
for g in up314:
    if g in state_of.index and state_of[g] in POLY:
        hits.append({"agi": g, "state": f"S{state_of[g]}",
                     "log2fc": de314.loc[g, "log2fc"],
                     "padj": de314.loc[g, "padj"],
                     "desc": desc_map.get(g, "")})
hits_df = pd.DataFrame(hits).sort_values("log2fc", ascending=False)
hits_df.to_csv(f"{OUT}/osd314_up_polycomb_genes.csv", index=False)
print(f"  {len(hits)} up-regulated genes in Polycomb states:")
print(hits_df.to_string(index=False))

# ---------------------------------------------------------------- figures
print("[6] Figures ...")
# fig17: null histograms
fig, axes = plt.subplots(2, 2, figsize=(10, 7))
for ax, s in zip(axes.flat, POLY):
    ax.hist(null[s], bins=50, color="#9e9e9e", alpha=0.8)
    ax.axvline(emp[f"S{s}"]["observed_odds"], color="#FF9400", lw=2,
               label=f"observed = {emp[f'S{s}']['observed_odds']:.1f}")
    ax.axvline(np.quantile(null[s], 0.99), color="#0279EE", ls="--", lw=1,
               label="null 99th pct")
    ax.set_title(f"S{s} — empirical p = {emp[f'S{s}']['empirical_p']:.4f}")
    ax.set_xlabel("Odds ratio from random 47-gene draws")
    ax.legend(frameon=False, fontsize=8)
fig.suptitle("Null calibration: random 47-gene sets vs OSD-314 observed")
fig.tight_layout()
fig.savefig(f"{FIG}/fig17_null_calibration.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig17_null_calibration.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# fig18: downsampling distribution
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.hist(ds[11], bins=40, color="#9e9e9e", alpha=0.8,
        label="OSD-37 up genes downsampled to n=47")
ax.axvline(emp["S11"]["observed_odds"], color="#FF9400", lw=2,
           label=f"OSD-314 observed = {emp['S11']['observed_odds']:.1f}")
ax.axvline(np.quantile(ds[11], 0.99), color="#0279EE", ls="--", lw=1,
           label="downsampled 99th pct")
ax.set_xlabel("S11 odds ratio")
ax.set_ylabel("Draws")
ax.set_title("Small-N artifact control")
ax.legend(frameon=False, fontsize=8)
fig.savefig(f"{FIG}/fig18_downsampling_control.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig18_downsampling_control.png", dpi=150,
            bbox_inches="tight")
plt.close(fig)

# fig19: threshold grid heatmap
pv = grid.set_index(["padj", "lfc"])[[f"S{s}" for s in POLY]]
fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(pv, annot=True, fmt=".1f", cmap="RdBu_r", center=1, ax=ax,
            cbar_kws={"label": "Odds ratio (up-regulated genes)"})
ax.set_title("OSD-314 Polycomb-state enrichment across DE thresholds")
ax.set_ylabel("(padj, |log2fc| threshold)")
fig.savefig(f"{FIG}/fig19_threshold_grid.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig19_threshold_grid.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# fig20: running enrichment plot
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(range(len(run)), run, color="#0279EE", lw=1.5)
ax.axhline(0, color="k", lw=0.5)
hits_pos = np.where(hit)[0]
ax.vlines(hits_pos, 0, run[hits_pos], color="#FD9BED", alpha=0.3, lw=0.5)
ax.set_xlabel("Genes ranked by log2 fold change (up -> left)")
ax.set_ylabel("Running enrichment sum")
ax.set_title(f"Polycomb-state genes (S11/S12/S13/S15, n={n_hits}) in "
             f"OSD-314 — permutation p = {perm_p:.4f}")
fig.savefig(f"{FIG}/fig20_running_enrichment.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig20_running_enrichment.png", dpi=150,
            bbox_inches="tight")
plt.close(fig)

# fig21: dot plot of up genes in Polycomb states
if len(hits_df):
    fig, ax = plt.subplots(figsize=(8, 0.4 * len(hits_df) + 2))
    colors = {"S11": "#0279EE", "S12": "#FF9400", "S13": "#75A025",
              "S15": "#FD9BED"}
    ax.scatter(hits_df["log2fc"], range(len(hits_df)),
               c=[colors.get(s, "#9e9e9e") for s in hits_df["state"]], s=60)
    ax.set_yticks(range(len(hits_df)),
                  [f"{r.agi} ({r.state})" for r in hits_df.itertuples()],
                  fontsize=7)
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("log2 fold change (spaceflight vs control)")
    ax.set_title("OSD-314 up-regulated genes in Polycomb states")
    fig.savefig(f"{FIG}/fig21_polycomb_up_genes.svg", bbox_inches="tight")
    fig.savefig(f"{FIG}/fig21_polycomb_up_genes.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

print("Done.")
