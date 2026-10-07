"""QC re-run of the root replication with genotype-correct contrasts (2026-10 audit).

The recovered BioSample titles show that OSD-218 is 16 Col-0 + 16 WS (the v8
parser left WS as NaN, so all 32 samples entered the "Col-0" contrast), and
that OSD-406 is Col-0, WS and sku5 roots flown on two different suborbital
rockets (Virgin Galactic VP-03 and Blue Origin NS-12). This script:

  1. reproduces the original v8 calls (pipeline check: must match Table S29),
  2. re-runs each dataset restricted to Col-0, and with genotype as a covariate,
  3. writes results/root_replication/qc_genotype_rerun.csv.

DE and enrichment logic is the same as root_replication_recovered.run().
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import os
import re
import shutil
import sys

sys.path.insert(0, f"{ROOT}/code")
import numpy as np
import pandas as pd
import chromdecode as cd
import run_osdr_sweep as sw

OUT = f"{ROOT}/results/root_replication"
SWEEP = f"{ROOT}/results/sweep"
os.makedirs(OUT, exist_ok=True)
os.makedirs(SWEEP, exist_ok=True)
sw.SWEEP = SWEEP

# Inputs shipped in data/derived/ (Tables D1, D5, D7).
for fn, dst in [("gene_state_assignments.csv", f"{ROOT}/results"),
                ("OSD-218_sample_meta.csv", OUT),
                ("OSD-406_sample_meta.csv", OUT)]:
    if not os.path.exists(f"{dst}/{fn}"):
        shutil.copy(f"{ROOT}/data/derived/{fn}", dst)

tss_state = pd.read_csv(f"{ROOT}/results/gene_state_assignments.csv")
st = tss_state.dropna(subset=["state_group"])
BG = st["agi"].unique()


def genotype_from_title(t):
    t = t.lower()
    for pat, g in [(r"\bcol-?0\b", "Col-0"), (r"\bws\b", "WS"),
                   (r"\bsku5\b", "sku5"), (r"\bsku6\b", "sku6")]:
        if re.search(pat, t):
            return g
    return None


def rocket_from_title(t):
    t = t.lower()
    if "virgin galactic" in t:
        return "VirginGalactic"
    if "blue origin" in t:
        return "BlueOrigin"
    return None


def ensure_counts(acc):
    path = f"{SWEEP}/{acc}_counts.csv"
    if not os.path.exists(path):
        assays = sw.list_arabidopsis_rnaseq()
        url = sw.get_counts_url(acc, assays[acc])
        with open(path, "wb") as fh:
            fh.write(sw.fetch(url, 300))
    return path


def median_ratio_scale(counts):
    pos = counts[counts.sum(axis=1) > 0]
    ref = pos.replace(0, np.nan).median(axis=1)
    sf = pos.div(ref, axis=0).median(axis=0)
    return counts.div(sf / sf.mean(), axis=1)


def run(acc, label, restrict=None, covariate_cols=()):
    counts = pd.read_csv(ensure_counts(acc), index_col=0)
    counts.index = [cd.clean_agi(i) or i for i in counts.index]
    counts = counts[~counts.index.str.startswith(("GLDS", "SRR"))]
    meta = pd.read_csv(f"{OUT}/{acc}_sample_meta.csv").set_index("sample")
    meta["genotype_qc"] = meta["title"].map(genotype_from_title)
    meta["rocket"] = meta["title"].map(rocket_from_title)
    keep = [s for s in counts.columns if s in meta.index]
    counts, meta = counts[keep], meta.loc[keep]
    if restrict:
        mask = meta.apply(restrict, axis=1)
        keep = mask[mask].index.tolist()
        counts, meta = counts[keep], meta.loc[keep]
    cond = (meta["condition"] == "Spaceflight").astype(float)
    cov = {c: meta[c].values for c in covariate_cols
           if meta[c].notna().any() and meta[c].dropna().nunique() > 1}
    de = cd.differential_expression_any(
        np.log2(median_ratio_scale(counts).clip(lower=1)), cond.values,
        covariates=cov or None)
    res = {"accession": acc, "analysis": label, "n": len(keep),
           "n_flt": int(cond.sum()), "n_gc": int((1 - cond).sum()),
           "genotypes": ";".join(sorted(meta["genotype_qc"].dropna().unique())),
           "covariates": ";".join(cov), "n_de_up": int((de["direction"] == "up").sum()),
           "n_de_down": int((de["direction"] == "down").sum()),
           "pc_odds_down": np.nan, "pc_padj_down": np.nan}
    de.to_csv(f"{OUT}/qc_{acc}_{label}_DE.csv")
    d = de[de["direction"] == "down"]
    if len(d):
        enr = cd.state_enrichment(st.rename(columns={"state_group": "st"}),
                                  d.index, background=BG, state_col="st")
        m = enr[enr["st"] == "Polycomb_repressed"]
        if len(m):
            res["pc_odds_down"] = m.iloc[0]["odds_ratio"]
            res["pc_padj_down"] = m.iloc[0]["padj"]
    print(res)
    return res


rows = [
    # 1. reproduce v8 exactly (Table S29)
    run("OSD-218", "v8_as_published", covariate_cols=("age_days",)),
    run("OSD-406", "v8_as_published"),
    # 2. genotype-correct contrasts
    run("OSD-218", "col0_only", restrict=lambda r: r["genotype_qc"] == "Col-0",
        covariate_cols=("age_days",)),
    run("OSD-218", "ws_only", restrict=lambda r: r["genotype_qc"] == "WS",
        covariate_cols=("age_days",)),
    run("OSD-218", "all_genotype_covariate",
        covariate_cols=("age_days", "genotype_qc")),
    run("OSD-406", "col0_virgin_galactic",
        restrict=lambda r: r["genotype_qc"] == "Col-0"),
    run("OSD-406", "all_genotype_rocket_covariates",
        covariate_cols=("genotype_qc", "rocket")),
]
pd.DataFrame(rows).to_csv(f"{OUT}/qc_genotype_rerun.csv", index=False)
print("saved qc_genotype_rerun.csv")
