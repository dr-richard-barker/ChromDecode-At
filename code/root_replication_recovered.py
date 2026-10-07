"""Run root replication tests on OSD-193, 218, 406 (metadata recovered)."""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys
sys.path.insert(0, f"{ROOT}/code")
import numpy as np
import pandas as pd
import chromdecode as cd

OUT = f"{ROOT}/results/root_replication"
SWEEP = f"{ROOT}/results/sweep"
tss_state = pd.read_csv(f"{ROOT}/results/gene_state_assignments.csv")
st = tss_state.dropna(subset=["state_group"])
BG = st["agi"].unique()

def median_ratio_scale(counts):
    pos = counts[counts.sum(axis=1) > 0]
    ref = pos.replace(0, np.nan).median(axis=1)
    sf = pos.div(ref, axis=0).median(axis=0)
    return counts.div(sf / sf.mean(), axis=1)

def run(acc, restrict=None, covariate_cols=("age_days",)):
    counts = pd.read_csv(f"{SWEEP}/{acc}_counts.csv", index_col=0)
    counts.index = [cd.clean_agi(i) or i for i in counts.index]
    counts = counts[~counts.index.str.startswith(("GLDS", "SRR"))]
    meta = pd.read_csv(f"{OUT}/{acc}_sample_meta.csv").set_index("sample")
    keep = [s for s in counts.columns if s in meta.index]
    counts = counts[keep]
    meta = meta.loc[keep]
    if restrict:
        mask = meta.apply(restrict, axis=1)
        keep = mask[mask].index.tolist()
        counts = counts[keep]
        meta = meta.loc[keep]
    cond = (meta["condition"] == "Spaceflight").astype(float)
    print(acc, "n=", len(keep), "| FLT:", int(cond.sum()),
          "GC:", int((1 - cond).sum()))
    cov = {}
    for c in covariate_cols:
        if c in meta.columns and meta[c].notna().any() \
                and meta[c].dropna().nunique() > 1:
            cov[c] = meta[c].values
    print("  covariates:", list(cov))
    de = cd.differential_expression_any(
        np.log2(median_ratio_scale(counts).clip(lower=1)), cond.values,
        covariates=cov or None)
    n_up = int((de["direction"] == "up").sum())
    n_dn = int((de["direction"] == "down").sum())
    print(f"  up {n_up} down {n_dn}")
    res = {"accession": acc, "n": len(keep), "n_de_up": n_up,
           "n_de_down": n_dn}
    de.to_csv(f"{OUT}/{acc}_DE.csv")
    d = de[de["direction"] == "down"]
    if len(d):
        enr = cd.state_enrichment(st.rename(columns={"state_group": "st"}),
                                  d.index, background=BG, state_col="st")
        enr.to_csv(f"{OUT}/{acc}_enrichment_down.csv", index=False)
        m = enr[enr["st"] == "Polycomb_repressed"]
        if len(m):
            res["pc_odds_down"] = m.iloc[0]["odds_ratio"]
            res["pc_padj_down"] = m.iloc[0]["padj"]
            print("  Polycomb_down odds:",
                  round(m.iloc[0]["odds_ratio"], 2), "padj:",
                  m.iloc[0]["padj"])
    return res

rows = []
# OSD-193: Col-0 WT only, age covariate
rows.append(run("OSD-193",
                restrict=lambda r: r["genotype"] == "Col-0",
                covariate_cols=("age_days",)))
# OSD-218: all Col-0, age covariate
rows.append(run("OSD-218", covariate_cols=("age_days",)))
# OSD-406: Col-0, suborbital (Virgin Galactic); no age covariate parsed
rows.append(run("OSD-406", covariate_cols=()))
pd.DataFrame(rows).to_csv(f"{OUT}/replication_recovered.csv", index=False)
print("saved replication_recovered.csv")
