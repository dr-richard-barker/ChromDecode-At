"""Power analysis: replicates needed for a stable OSD-120 Polycomb enrichment test.

Step 1: analytic Fisher-test power vs number of down DEGs (odds 3.1 / 7.0).
Step 2: empirical subsampling of OSD-218 (16v16 Col-0 roots) at k reps/group.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys
sys.path.insert(0, f"{ROOT}/code")
import numpy as np
import pandas as pd
from scipy import stats as sps
import chromdecode as cd

OUT = f"{ROOT}/results/root_replication"
rng = np.random.default_rng(42)

# background Polycomb fraction among expressed, state-assigned genes
tss = pd.read_csv(f"{ROOT}/results/gene_state_assignments.csv")
st = tss.dropna(subset=["state_group"])
n_bg = st["agi"].nunique()
n_pc = st[st["state_group"] == "Polycomb_repressed"]["agi"].nunique()
p0 = n_pc / n_bg
print(f"background: {n_pc}/{n_bg} = {p0:.3f} Polycomb")

# ---------------------------------------------------------------- step 1
def fisher_power(n_deg, odds, p0, n_draws=10000, n_tests=5):
    """Power of the enrichment test given n_deg down genes."""
    p1 = odds * p0 / (1 + odds * p0 - p0)  # convert odds -> fraction
    hits = rng.binomial(n_deg, p1, n_draws)
    a = hits
    b = n_deg - hits
    c = rng.binomial(n_bg - n_deg, p0, n_draws)  # approx: bg Polycomb
    d = (n_bg - n_deg) - c
    ps = np.empty(n_draws)
    for i in range(n_draws):
        ps[i] = sps.fisher_exact([[a[i], b[i]], [c[i], d[i]]],
                                 alternative="greater")[1]
    padj = cd.bh_adjust(ps) if n_tests > 1 else ps
    return float((padj < 0.05).mean())

grid = [5, 10, 15, 20, 30, 40, 60, 80, 120, 200, 300, 500]
rows = []
for odds in [3.1, 7.0]:
    for n in grid:
        pw = fisher_power(n, odds, p0)
        rows.append({"component": "analytic", "odds": odds, "k_reps": None,
                     "n_deg": n, "power": pw})
        print(f"analytic odds={odds} n_deg={n}: power={pw:.3f}")
analytic = pd.DataFrame(rows)
for odds in [3.1, 7.0]:
    sub = analytic[analytic["odds"] == odds]
    need = sub[sub["power"] >= 0.8]["n_deg"].min()
    print(f"--> odds {odds}: n_DEG needed for 80% power = {need}")

# ---------------------------------------------------------------- step 2
counts = pd.read_csv(f"{ROOT}/results/sweep/OSD-218_counts.csv",
                     index_col=0)
counts.index = [cd.clean_agi(i) or i for i in counts.index]
counts = counts[~counts.index.str.startswith(("GLDS", "SRR"))]
meta = pd.read_csv(f"{OUT}/OSD-218_sample_meta.csv").set_index("sample")
keep = [c for c in counts.columns if c in meta.index]
counts = counts[keep]
meta = meta.loc[keep]
flt = meta[meta["condition"] == "Spaceflight"].index.tolist()
gc = meta[meta["condition"] == "Ground Control"].index.tolist()
age = meta["age_days"]

def median_ratio_scale(c):
    pos = c[c.sum(axis=1) > 0]
    ref = pos.replace(0, np.nan).median(axis=1)
    sf = pos.div(ref, axis=0).median(axis=0)
    return c.div(sf / sf.mean(), axis=1)

def run_draw(k, nuisance=0, n_draws=20):
    """Subsample k reps/group; return n_down and Polycomb padj."""
    res = []
    st_g = st.rename(columns={"state_group": "st"})
    bg = st_g["agi"].unique()
    for d in range(n_draws):
        f_s = rng.choice(flt, k, replace=False)
        g_s = rng.choice(gc, k, replace=False)
        cols = list(f_s) + list(g_s)
        sub = counts[cols]
        cond = np.array([1] * k + [0] * k, dtype=float)  # FLT=1, GC=0
        cov = {"age_days": age.reindex(cols).values}
        # nuisance covariates: random split to mimic light/day design burden
        for j in range(nuisance):
            cov[f"nuis{j}"] = rng.integers(0, 2, len(cols)).astype(float)
        de = cd.differential_expression_any(
            np.log2(median_ratio_scale(sub).clip(lower=1)), cond,
            covariates=cov)
        dn = de[de["direction"] == "down"]
        if len(dn) < 3:
            res.append({"k": k, "nuisance": nuisance, "draw": d,
                        "n_down": len(dn), "pc_padj": np.nan})
            continue
        enr = cd.state_enrichment(st_g, dn.index, background=bg,
                                  state_col="st")
        m = enr[enr["st"] == "Polycomb_repressed"]
        res.append({"k": k, "nuisance": nuisance, "draw": d,
                    "n_down": len(dn),
                    "pc_padj": float(m.iloc[0]["padj"]) if len(m) else np.nan})
    return res

ks = [2, 3, 4, 6, 8, 12, 16]
emp_rows = []
for k in ks:
    for nuis, tag in [(0, "clean")] + (
            [(2, "osd120_burden")] if k >= 3 else []):
        rr = run_draw(k, nuisance=nuis)
        emp_rows += rr
        df = pd.DataFrame(rr)
        pw = float((df["pc_padj"] < 0.05).mean())
        print(f"k={k} nuis={nuis}: mean n_down="
              f"{df['n_down'].mean():.0f}, power={pw:.2f}")
emp = pd.DataFrame(emp_rows)

summary = emp.groupby(["k", "nuisance"]).agg(
    n_draws=("n_down", "size"),
    mean_n_down=("n_down", "mean"),
    power=("pc_padj", lambda s: float((s < 0.05).mean())),
).reset_index()
summary["component"] = "empirical"
summary["odds"] = np.nan
summary["n_deg"] = summary["mean_n_down"]
print(summary.to_string(index=False))

out = pd.concat([analytic, emp], ignore_index=True)
out.to_csv(f"{OUT}/osd120_power_analysis.csv", index=False)
summary.to_csv(f"{OUT}/osd120_power_summary.csv", index=False)
print("saved")
