"""QC re-run of the OSDR sweep DE datasets and the OSD-314 design (2026-10 audit).

1. Re-runs run_osdr_sweep.process_dataset() for the five datasets with full
   DE (Table S10) and checks n_up/n_down and the 12-/36-state odds against
   Tables S10-S12 (pipeline check). This also regenerates the corrected
   5-group enrichment that was missing from the archive (the shipped
   sweep_enrichment_long.csv predates the log2 correction).
2. OSD-314 has three gravity arms (0g, 1g, 0.3g Mars) crossed with red/dark
   light. The sweep's token parser does not read "03g", so those samples have
   no parsed level; this re-run codes gravity and light explicitly from the
   sample names and fits 0g vs 1g and 0.3g vs 1g with light as a covariate.

Outputs (results/qc/): sweep_rerun_check.csv, sweep_enrichment5_corrected.csv,
osd314_design_rerun.csv, plus per-contrast DE tables.
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

QC = f"{ROOT}/results/qc"
SWEEP = f"{ROOT}/results/sweep"
os.makedirs(QC, exist_ok=True)
os.makedirs(SWEEP, exist_ok=True)
sw.SWEEP = SWEEP
if not os.path.exists(f"{ROOT}/results/gene_state_assignments.csv"):
    shutil.copy(f"{ROOT}/data/derived/gene_state_assignments.csv", f"{ROOT}/results")
tss = pd.read_csv(f"{ROOT}/results/gene_state_assignments.csv")

TABLES = f"{ROOT}/supplementary_tables"
s10 = pd.read_csv(f"{TABLES}/sweep_summary.csv").set_index("accession")
s11 = pd.read_csv(f"{TABLES}/sweep_enrichment12_long.csv")
s12 = pd.read_csv(f"{TABLES}/sweep_enrichment36_long.csv")


def enrich(de, kind, col):
    st = tss.dropna(subset=[col]).rename(columns={col: "st"})[["agi", "st"]]
    if col == "pcsd_state":
        st["st"] = st["st"].astype(int)
    genes = de[de["direction"] == kind].index
    if len(genes) == 0:
        return pd.DataFrame(columns=["st", "odds_ratio", "padj", "n_in_set"])
    return cd.state_enrichment(st, genes, background=st["agi"].unique(),
                               state_col="st")


def lookup(e, state):
    m = e[e["st"] == state]
    return (m.iloc[0]["odds_ratio"], m.iloc[0]["padj"]) if len(m) else (np.nan, np.nan)


# ---- 1. sweep re-run --------------------------------------------------------
assays = sw.list_arabidopsis_rnaseq()
DE_SETS = [a for a in s10.index if s10.loc[a, "status"] == "ok"]
check, long5 = [], []
for acc in DE_SETS:
    r = sw.process_dataset(acc, assays[acc], tss, None)
    de = pd.read_csv(f"{SWEEP}/{acc}_DE.csv", index_col=0)
    row = {"accession": acc, "status": r["status"],
           "n_up": r["n_de_up"], "n_down": r["n_de_down"],
           "S10_n_up": s10.loc[acc, "n_de_up"], "S10_n_down": s10.loc[acc, "n_de_down"]}
    for kind in ("up", "down"):
        e5 = enrich(de, kind, "state_group")
        for _, er in e5.iterrows():
            long5.append({"accession": acc, "kind": kind, "state_group": er["st"],
                          "odds_ratio": er["odds_ratio"], "padj": er["padj"],
                          "n_set": er["n_in_set"]})
        e12, e36 = enrich(de, kind, "araencode_state"), enrich(de, kind, "pcsd_state")
        ref12 = s11[(s11.accession == acc) & (s11.kind == kind)].drop_duplicates()
        ref36 = s12[(s12.accession == acc) & (s12.kind == kind)]
        d12 = [abs(lookup(e12, x.state)[0] - x.odds_ratio) for x in ref12.itertuples()
               if np.isfinite(x.odds_ratio)]
        d36 = [abs(lookup(e36, int(x.state))[0] - x.odds_ratio) for x in ref36.itertuples()
               if np.isfinite(x.odds_ratio)]
        row[f"max_abs_diff_12state_{kind}"] = max(d12) if d12 else np.nan
        row[f"max_abs_diff_36state_{kind}"] = max(d36) if d36 else np.nan
    print(row)
    check.append(row)
pd.DataFrame(check).to_csv(f"{QC}/sweep_rerun_check.csv", index=False)
pd.DataFrame(long5).to_csv(f"{QC}/sweep_enrichment5_corrected.csv", index=False)

# ---- 2. OSD-314 explicit design ---------------------------------------------
counts = pd.read_csv(f"{SWEEP}/OSD-314_counts.csv", index_col=0)
counts.index = [cd.clean_agi(i) or i for i in counts.index]
counts = counts[~counts.index.str.startswith(("GLDS", "SRR"))]


def gravity(s):
    m = re.search(r"_(0g|1g|03g)_", s, re.I)
    return {"0g": "0g", "1g": "1g", "03g": "0.3g"}[m.group(1).lower()] if m else None


def light(s):
    return "RED" if "_RED_" in s.upper() else "DARK" if "_DARK_" in s.upper() else None


design = pd.DataFrame({"gravity": [gravity(s) for s in counts.columns],
                       "light": [light(s) for s in counts.columns]},
                      index=counts.columns)
design.to_csv(f"{QC}/osd314_sample_design.csv")
print(design.value_counts().to_string())

# What the sweep's coding did: samples without a parsed 0G/1G token.
fac = sw.sample_factors(sw.get_metadata("OSD-314"), list(counts.columns))
contrast = cd.detect_contrast_from_factors(fac)
levels = contrast["levels"].reindex(counts.columns)
sweep_cond = (levels != contrast["control"]).astype(int)
design["sweep_parsed_level"] = levels
design["sweep_coded_as_treatment"] = sweep_cond
design.to_csv(f"{QC}/osd314_sample_design.csv")

rows = []
for treat in ("0g", "0.3g"):
    for use_light in (False, True):
        sub = design[design["gravity"].isin(["1g", treat])]
        cond = (sub["gravity"] == treat).astype(float).values
        cov = {"light": sub["light"].values} if use_light else None
        de = cd.differential_expression_any(np.log2(counts[sub.index].clip(lower=1)),
                                            cond, covariates=cov)
        tag = f"{treat.replace('.', '')}_vs_1g" + ("_light" if use_light else "")
        de.to_csv(f"{QC}/osd314_{tag}_DE.csv")
        row = {"contrast": tag, "n_treat": int(cond.sum()), "n_ctrl": int((1 - cond).sum()),
               "n_up": int((de.direction == "up").sum()),
               "n_down": int((de.direction == "down").sum())}
        e5d, e5u = enrich(de, "down", "state_group"), enrich(de, "up", "state_group")
        e12d, e36u = enrich(de, "down", "araencode_state"), enrich(de, "up", "pcsd_state")
        row["polycomb5_down_odds"], row["polycomb5_down_padj"] = lookup(e5d, "Polycomb_repressed")
        row["polycomb5_up_odds"], row["polycomb5_up_padj"] = lookup(e5u, "Polycomb_repressed")
        row["repression12_down_odds"], row["repression12_down_padj"] = lookup(e12d, "Repression")
        row["S11_up_odds"], row["S11_up_padj"] = lookup(e36u, 11)
        print(row)
        rows.append(row)
pd.DataFrame(rows).to_csv(f"{QC}/osd314_design_rerun.csv", index=False)
print("done")
