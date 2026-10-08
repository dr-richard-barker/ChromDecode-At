"""Re-run the power analysis and the OSD-120 module bound with OSD-218 Col-0 only.

power_analysis.py (Tables S42-S43) subsampled all 32 OSD-218 samples, but 16 of
them are WS, not Col-0. osd120_stratified.py (Table S44) took its "Polycomb-down
module" from that pooled OSD-218 run. This script:

  1. reproduces S42/S43 (analytic + pooled empirical; same code, seed 42) and
     the S44 module-bound rows (pooled OSD-218 module), as a check;
  2. repeats the empirical power curve on the 16 Col-0 samples only (8 v 8, so
     k <= 8; k = 8 is the full Col-0 set and must reproduce Table S49);
  3. recomputes the OSD-120 module bound with the Col-0-only OSD-218 module.

The OSD-120 stratified DE and interaction model (S44 stratified/interaction rows,
S45-S47) use OSD-120 Col-0 samples only and are unchanged; they are recomputed
here only because the module bound needs them.

Inputs: results/sweep/OSD-218_counts.csv and OSD-120_counts.csv,
results/root_replication/OSD-218_sample_meta.csv and qc_OSD-218_*_DE.csv
(from qc_genotype_rerun.py), data/derived/. Outputs: results/qc/.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import re
import sys

sys.path.insert(0, f"{ROOT}/code")
import numpy as np
import pandas as pd
from scipy import stats as sps
import chromdecode as cd

QC = f"{ROOT}/results/qc"
RR = f"{ROOT}/results/root_replication"
SWEEP = f"{ROOT}/results/sweep"
T = f"{ROOT}/supplementary_tables"

tss = pd.read_csv(f"{ROOT}/data/derived/gene_state_assignments.csv")
st = tss.dropna(subset=["state_group"])
n_bg = st["agi"].nunique()
n_pc = st[st["state_group"] == "Polycomb_repressed"]["agi"].nunique()
p0 = n_pc / n_bg
st_g = st.rename(columns={"state_group": "st"})
BG = st_g["agi"].unique()


def median_ratio_scale(c):
    pos = c[c.sum(axis=1) > 0]
    ref = pos.replace(0, np.nan).median(axis=1)
    sf = pos.div(ref, axis=0).median(axis=0)
    return c.div(sf / sf.mean(), axis=1)


# ================================================================ power
def fisher_power(rng, n_deg, odds, n_draws=10000, n_tests=5):
    p1 = odds * p0 / (1 + odds * p0 - p0)
    hits = rng.binomial(n_deg, p1, n_draws)
    a, b = hits, n_deg - hits
    c = rng.binomial(n_bg - n_deg, p0, n_draws)
    d = (n_bg - n_deg) - c
    ps = np.array([sps.fisher_exact([[a[i], b[i]], [c[i], d[i]]], alternative="greater")[1]
                   for i in range(n_draws)])
    return float((cd.bh_adjust(ps) < 0.05).mean())


counts218 = pd.read_csv(f"{SWEEP}/OSD-218_counts.csv", index_col=0)
counts218.index = [cd.clean_agi(i) or i for i in counts218.index]
counts218 = counts218[~counts218.index.str.startswith(("GLDS", "SRR"))]
meta218 = pd.read_csv(f"{RR}/OSD-218_sample_meta.csv").set_index("sample")
meta218["genotype_qc"] = np.where(meta218["title"].str.contains(r"\bCol-0\b"), "Col-0",
                                  np.where(meta218["title"].str.contains(r"\bWS\b"), "WS", None))


def empirical(rng, pool, ks):
    cnt = counts218[[c for c in counts218.columns if c in pool.index]]
    meta = pool.loc[cnt.columns]
    flt = meta[meta["condition"] == "Spaceflight"].index.tolist()
    gc = meta[meta["condition"] == "Ground Control"].index.tolist()
    age = meta["age_days"]
    rows = []
    for k in ks:
        for nuis in [0] + ([2] if k >= 3 else []):
            for d in range(20):
                cols = list(rng.choice(flt, k, replace=False)) + list(rng.choice(gc, k, replace=False))
                cond = np.array([1] * k + [0] * k, dtype=float)
                cov = {"age_days": age.reindex(cols).values}
                for j in range(nuis):
                    cov[f"nuis{j}"] = rng.integers(0, 2, len(cols)).astype(float)
                de = cd.differential_expression_any(
                    np.log2(median_ratio_scale(cnt[cols]).clip(lower=1)), cond, covariates=cov)
                dn = de[de["direction"] == "down"]
                if len(dn) < 3:
                    rows.append({"k": k, "nuisance": nuis, "draw": d, "n_down": len(dn), "pc_padj": np.nan})
                    continue
                enr = cd.state_enrichment(st_g, dn.index, background=BG, state_col="st")
                m = enr[enr["st"] == "Polycomb_repressed"]
                rows.append({"k": k, "nuisance": nuis, "draw": d, "n_down": len(dn),
                             "pc_padj": float(m.iloc[0]["padj"]) if len(m) else np.nan})
    emp = pd.DataFrame(rows)
    summ = emp.groupby(["k", "nuisance"]).agg(
        n_draws=("n_down", "size"), mean_n_down=("n_down", "mean"),
        power=("pc_padj", lambda s: float((s < 0.05).mean()))).reset_index()
    return emp, summ


# 1. reproduction of S42/S43: analytic first, then pooled empirical, one RNG stream
rng = np.random.default_rng(42)
analytic = pd.DataFrame([{"component": "analytic", "odds": o, "n_deg": n, "power": fisher_power(rng, n, o)}
                         for o in [3.1, 7.0] for n in [5, 10, 15, 20, 30, 40, 60, 80, 120, 200, 300, 500]])
_, summ_pooled = empirical(rng, meta218, [2, 3, 4, 6, 8, 12, 16])
s43 = pd.read_csv(f"{T}/osd120_power_summary.csv")
s42 = pd.read_csv(f"{T}/osd120_power_analysis.csv")
chk = {
    "S42_analytic_max_abs_diff": float(np.abs(
        analytic.sort_values(["odds", "n_deg"])["power"].values
        - s42[s42.component == "analytic"].sort_values(["odds", "n_deg"])["power"].values).max()),
    "S43_power_max_abs_diff": float(np.abs(summ_pooled["power"].values - s43["power"].values).max()),
    "S43_mean_n_down_max_abs_diff": float(np.abs(summ_pooled["mean_n_down"].values - s43["mean_n_down"].values).max()),
}
print("power reproduction:", chk)

# 2. Col-0 only (8 v 8), fresh seed
rng = np.random.default_rng(42)
emp_c, summ_c = empirical(rng, meta218[meta218["genotype_qc"] == "Col-0"], [2, 3, 4, 6, 8])
summ_c.insert(0, "pool", "OSD-218 Col-0 only (8 v 8)")
summ_pooled.insert(0, "pool", "OSD-218 pooled Col-0 + WS (16 v 16; as published)")
col0_full = pd.read_csv(f"{RR}/qc_OSD-218_col0_only_DE.csv", index_col=0)
k8 = summ_c[(summ_c.k == 8) & (summ_c.nuisance == 0)].iloc[0]
chk["col0_k8_mean_n_down_vs_S49"] = f"{k8.mean_n_down:.1f} vs {int((col0_full.direction == 'down').sum())}"
print(summ_c.to_string(index=False))
pd.concat([summ_pooled, summ_c], ignore_index=True).to_csv(f"{QC}/power_summary_col0.csv", index=False)
emp_c.to_csv(f"{QC}/power_draws_col0.csv", index=False)

# ================================================================ OSD-120 bound
counts = pd.read_csv(f"{SWEEP}/OSD-120_counts.csv", index_col=0)
counts.index = [cd.clean_agi(i) or i for i in counts.index]
counts = counts[~counts.index.str.startswith(("GLDS", "SRR"))]
cols = [c for c in counts.columns if "Col-0" in c and "PhyD" not in c and "WS" not in c]
counts = counts[cols]
cond = np.array([1.0 if "_FLT_" in c else 0.0 for c in cols])
light = np.array([1.0 if "_Alight_" in c else 0.0 for c in cols])
expr = np.log2(median_ratio_scale(counts).clip(lower=1))

strat_de = {name: cd.differential_expression_any(expr.loc[:, mask], cond[mask])
            for name, mask in [("Alight", light == 1), ("dark", light == 0)]}
X = np.column_stack([np.ones(len(cols)), cond, light, cond * light])
Y = expr.values.T
betas, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
resid = Y - X @ betas
dof = X.shape[0] - np.linalg.matrix_rank(X)
se = np.sqrt(np.maximum(((resid ** 2).sum(axis=0) / dof) * np.linalg.pinv(X.T @ X)[3, 3], 1e-12))
p_i = 2 * sps.t.sf(np.abs(betas[3] / se), dof)
inter = pd.DataFrame({"interaction_log2fc": betas[3], "p_value": p_i}, index=expr.index)
inter["padj"] = cd.bh_adjust(p_i)

redox = set(pd.read_csv(f"{ROOT}/data/derived/redox_module_genes.csv").iloc[:, 0])
pc_agi = set(st[st["state_group"] == "Polycomb_repressed"]["agi"])


def module_row(label, genes):
    g = [x for x in genes if x in inter.index]
    sub = inter.loc[g]
    a, d = strat_de["Alight"], strat_de["dark"]
    dA, dD = (a["direction"] == "down").reindex(g).fillna(False), (d["direction"] == "down").reindex(g).fillna(False)
    return {"module": label, "n_genes": len(g),
            "inter_padj_lt05": int((sub["padj"] < 0.05).sum()),
            "inter_max_abs_lfc": float(sub["interaction_log2fc"].abs().max()),
            "inter_median_abs_lfc": float(sub["interaction_log2fc"].abs().median()),
            "down_in_both_strata": int((dA & dD).sum()), "down_in_one_stratum": int((dA ^ dD).sum()),
            "down_in_neither": int((~dA & ~dD).sum()),
            "median_strat_lfc_Alight": float(a["log2fc"].reindex(g).median()),
            "median_strat_lfc_dark": float(d["log2fc"].reindex(g).median())}


rows = [module_row("redox_module", redox)]
for tag, fn in [("pc_down_module_pooled_as_published", "qc_OSD-218_v8_as_published_DE.csv"),
                ("pc_down_module_col0_only", "qc_OSD-218_col0_only_DE.csv"),
                ("pc_down_module_ws_only", "qc_OSD-218_ws_only_DE.csv")]:
    d218 = pd.read_csv(f"{RR}/{fn}", index_col=0)
    rows.append(module_row(tag, set(d218[(d218.direction == "down") & d218.index.isin(pc_agi)].index)))
bound = pd.DataFrame(rows)
print(bound.round(3).to_string(index=False))
bound.to_csv(f"{QC}/osd120_module_bound_col0.csv", index=False)

s44 = pd.read_csv(f"{T}/osd120_stratified.csv")
pub = s44[s44.component == "module_bound"].set_index("module")
for ours, theirs in [("redox_module", "redox_module"), ("pc_down_module_pooled_as_published", "pc_down_module")]:
    r = bound.set_index("module").loc[ours]
    chk[f"S44_{theirs}_n_genes_match"] = int(r.n_genes) == int(pub.loc[theirs, "n_genes"])
    chk[f"S44_{theirs}_median_abs_inter_diff"] = float(abs(r.inter_median_abs_lfc - pub.loc[theirs, "inter_median_abs_lfc"]))
    chk[f"S44_{theirs}_median_Alight_diff"] = float(abs(r.median_strat_lfc_Alight - pub.loc[theirs, "median_strat_lfc_Alight"]))
pd.Series(chk).to_csv(f"{QC}/col0_power_bound_reproduction_check.csv", header=["value"])
print("checks:", chk)
