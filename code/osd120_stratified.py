"""OSD-120 Col-0 WT: stratified flight-within-light analysis + interaction bound.

Design: balanced 2x2 factorial (FLT/GC x Alight/dark, 3 reps/cell, all Day 13).
1. Stratified DE: flight effect within each light regime (3v3).
2. Polycomb enrichment of down genes per stratum.
3. Gene-level interaction model: expr ~ cond + light + cond:light (12 samples).
4. Module bounds: redox module (23 genes) and OSD-218 Polycomb-down module.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys
sys.path.insert(0, f"{ROOT}/code")
import numpy as np
import pandas as pd
from scipy import stats
import chromdecode as cd

OUT = f"{ROOT}/results/root_replication"
SWEEP = f"{ROOT}/results/sweep"

tss_state = pd.read_csv(f"{ROOT}/results/gene_state_assignments.csv")
st = tss_state.dropna(subset=["state_group"])
BG = st["agi"].unique()

# ---- data ----
counts = pd.read_csv(f"{SWEEP}/OSD-120_counts.csv", index_col=0)
counts.index = [cd.clean_agi(i) or i for i in counts.index]
counts = counts[~counts.index.str.startswith(("GLDS", "SRR"))]
cols = [c for c in counts.columns
        if "Col-0" in c and "PhyD" not in c and "WS" not in c]
counts = counts[cols]
cond = np.array([1.0 if "_FLT_" in c else 0.0 for c in cols])
light = np.array([1.0 if "_Alight_" in c else 0.0 for c in cols])  # 1=Alight
print(f"n={len(cols)} | FLT {int(cond.sum())} GC {int((1-cond).sum())} | "
      f"Alight {int(light.sum())} dark {int((1-light).sum())}")
cells = pd.crosstab(cond, light)
print(cells)

def median_ratio_scale(counts):
    pos = counts[counts.sum(axis=1) > 0]
    ref = pos.replace(0, np.nan).median(axis=1)
    sf = pos.div(ref, axis=0).median(axis=0)
    return counts.div(sf / sf.mean(), axis=1)

expr = np.log2(median_ratio_scale(counts).clip(lower=1))

def ols_test(Y, X, coef_idx):
    """OLS per gene; returns beta, p for design column coef_idx."""
    betas, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ betas
    dof = X.shape[0] - np.linalg.matrix_rank(X)
    sigma2 = (resid ** 2).sum(axis=0) / dof
    XtX_inv = np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.maximum(sigma2 * XtX_inv[coef_idx, coef_idx], 1e-12))
    t = betas[coef_idx] / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    return betas[coef_idx], p

# ---- 1. stratified DE ----
strat_rows = []
strat_de = {}
for name, mask in [("Alight", light == 1), ("dark", light == 0)]:
    sub = expr.loc[:, mask]
    c = cond[mask]
    de = cd.differential_expression_any(sub, c)
    strat_de[name] = de
    n_up = int((de["direction"] == "up").sum())
    n_dn = int((de["direction"] == "down").sum())
    row = {"analysis": f"flight_within_{name}", "n": int(mask.sum()),
           "n_up": n_up, "n_down": n_dn}
    d = de[de["direction"] == "down"]
    if len(d):
        enr = cd.state_enrichment(st.rename(columns={"state_group": "st"}),
                                  d.index, background=BG, state_col="st")
        m = enr[enr["st"] == "Polycomb_repressed"]
        if len(m):
            row["pc_odds_down"] = m.iloc[0]["odds_ratio"]
            row["pc_padj_down"] = m.iloc[0]["padj"]
            row["pc_n_down_in_state"] = int(m.iloc[0]["n_in_set"])
    strat_rows.append(row)
    print(name, "->", row)

# ---- 2. interaction model (12 samples) ----
X = np.column_stack([np.ones(len(cols)), cond, light, cond * light])
Y = expr.values.T  # samples x genes
beta_i, p_i = ols_test(Y, X, 3)
inter = pd.DataFrame({"interaction_log2fc": beta_i, "p_value": p_i},
                     index=expr.index)
inter["padj"] = cd.bh_adjust(p_i)
n_inter_sig = int((inter["padj"] < 0.05).sum())
n_inter_big = int(((inter["padj"] < 0.05) &
                   (inter["interaction_log2fc"].abs() > 1)).sum())
q = inter["interaction_log2fc"].abs().quantile([0.5, 0.9, 0.95, 0.99])
print(f"interaction: sig {n_inter_sig}, sig&|lfc|>1 {n_inter_big}, "
      f"|lfc| quantiles {dict(q.round(2))}")

# ---- 3. module bounds ----
redox = pd.read_csv(f"{ROOT}/results/redox_module_genes.csv")
redox_genes = set(redox.iloc[:, 0])
# Polycomb-down module: OSD-218 down genes in Polycomb state
d218 = pd.read_csv(f"{OUT}/OSD-218_DE.csv", index_col=0)
pc_agi = set(st[st["state_group"] == "Polycomb_repressed"]["agi"])
pc_down_module = set(d218[(d218["direction"] == "down") &
                          d218.index.isin(pc_agi)].index)
print(f"modules: redox {len(redox_genes)}, pc_down {len(pc_down_module)}")

mod_rows = []
for mname, genes in [("redox_module", redox_genes),
                     ("pc_down_module", pc_down_module)]:
    g = [x for x in genes if x in inter.index]
    sub = inter.loc[g]
    # per-stratum log2fc for module genes
    la = strat_de["Alight"]["log2fc"].reindex(g)
    ld = strat_de["dark"]["log2fc"].reindex(g)
    mod_rows.append({
        "module": mname, "n_genes": len(g),
        "inter_padj_lt05": int((sub["padj"] < 0.05).sum()),
        "inter_max_abs_lfc": float(sub["interaction_log2fc"].abs().max()),
        "inter_median_abs_lfc": float(sub["interaction_log2fc"].abs().median()),
        "down_in_both_strata": int(((strat_de["Alight"]["direction"] == "down")
                                    & (strat_de["dark"]["direction"] == "down")
                                    ).reindex(g).fillna(False).sum()),
        "down_in_one_stratum": int((((strat_de["Alight"]["direction"] == "down")
                                     ^ (strat_de["dark"]["direction"] == "down"))
                                    ).reindex(g).fillna(False).sum()),
        "down_in_neither": int(((strat_de["Alight"]["direction"] != "down")
                                & (strat_de["dark"]["direction"] != "down")
                                ).reindex(g).fillna(False).sum()),
        "median_strat_lfc_Alight": float(la.median()),
        "median_strat_lfc_dark": float(ld.median()),
    })
    print(mod_rows[-1])

# ---- save ----
rows = []
for r in strat_rows:
    rows.append({"component": "stratified", **r})
rows.append({"component": "interaction", "analysis": "cond_x_light",
             "n_sig_padj05": n_inter_sig,
             "n_sig_lfc1": n_inter_big,
             "abs_lfc_q50": float(q.loc[0.5]), "abs_lfc_q90": float(q.loc[0.9]),
             "abs_lfc_q95": float(q.loc[0.95]), "abs_lfc_q99": float(q.loc[0.99])})
for r in mod_rows:
    rows.append({"component": "module_bound", **r})
pd.DataFrame(rows).to_csv(f"{OUT}/osd120_stratified.csv", index=False)
inter.to_csv(f"{OUT}/osd120_interaction_genelevel.csv")
for name in ["Alight", "dark"]:
    strat_de[name].to_csv(f"{OUT}/osd120_stratDE_{name}.csv")
print("saved")
