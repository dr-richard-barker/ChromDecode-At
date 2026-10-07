"""Patch runner for OSD-411 and OSD-624 (unnormalized counts only)."""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys, os, re, io, json, urllib.request
sys.path.insert(0, f"{ROOT}/code")
import numpy as np
import pandas as pd
import chromdecode as cd
import run_osdr_sweep as sw

OUT = f"{ROOT}/results/root_replication"
tss_state = pd.read_csv(f"{ROOT}/results/gene_state_assignments.csv")

def fetch(url, timeout=600):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read()

def median_ratio_scale(counts):
    """DESeq2-style median-of-ratios size-factor scaling (counts > 0)."""
    pos = counts[counts.sum(axis=1) > 0]
    geo_ref = pos.replace(0, np.nan).median(axis=1)
    ratios = pos.div(geo_ref, axis=0)
    sf = ratios.median(axis=0)
    sf = sf / sf.mean()
    return counts.div(sf, axis=1)

def run_unnorm(acc, assay, sample_rename=None):
    url = (f"https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/{acc}"
           f"/assay/{urllib.request.quote(assay, safe='()')}/files/?format=json")
    files = json.loads(fetch(url))[acc]["assays"][assay]["files"]
    cu = next((info["URL"] for k, info in files.items()
               if "RSEM_Unnormalized_Counts_GLbulkRNAseq" in k))
    counts = pd.read_csv(io.BytesIO(fetch(cu)), index_col=0)
    counts.index = [cd.clean_agi(i) or i for i in counts.index]
    counts = counts[~counts.index.str.startswith(("GLDS", "SRR"))]
    counts = median_ratio_scale(counts)
    if sample_rename:
        counts.columns = [sample_rename(c) for c in counts.columns]

    meta = sw.get_metadata(acc)
    fac = sw.sample_factors(meta, list(counts.columns))
    contrast = cd.detect_contrast_from_factors(fac) if fac is not None else None
    print(acc, "contrast:", None if contrast is None else
          (contrast["column"], contrast["control"], contrast["treatments"]))
    if contrast is None:
        return None
    col, ctrl, treats, levels = (contrast["column"], contrast["control"],
                                 contrast["treatments"], contrast["levels"])
    cond = (levels.reindex(counts.columns) != ctrl).astype(float).values
    cov = {c: fac[c].reindex(counts.columns).values
           for c in fac.columns if c != col}
    de = cd.differential_expression_any(np.log2(counts.clip(lower=1)),
                                        cond, covariates=cov)
    de.to_csv(f"{OUT}/{acc}_DE.csv")
    n_up = int((de["direction"] == "up").sum())
    n_dn = int((de["direction"] == "down").sum())
    print(acc, "up:", n_up, "down:", n_dn)
    res = {"accession": acc, "status": "ok", "n_de_up": n_up,
           "n_de_down": n_dn,
           "contrast": f"{col}: {ctrl} vs {treats[0]}"}
    st = tss_state.dropna(subset=["state_group"])
    bg = st["agi"].unique()
    for direction, fn in [("up", "enrichment_up"), ("down", "enrichment_down")]:
        d = de[de["direction"] == direction]
        if len(d) == 0:
            continue
        enr = cd.state_enrichment(
            st.rename(columns={"state_group": "st"}), d.index,
            background=bg, state_col="st")
        enr.to_csv(f"{OUT}/{acc}_{fn}.csv", index=False)
        if direction == "down":
            m = enr[enr["st"] == "Polycomb_repressed"]
            if len(m):
                res["pc_odds_down"] = m.iloc[0]["odds_ratio"]
                res["pc_padj_down"] = m.iloc[0]["padj"]
                print("  Polycomb_down odds:",
                      round(m.iloc[0]["odds_ratio"], 2),
                      "padj:", m.iloc[0]["padj"])
    return res

results = []
# OSD-411: Col-0 mature root / root tip
try:
    results.append(run_unnorm(
        "OSD-411",
        "OSD-411_transcription-profiling_rna-sequencing-(rna-seq)_Illumina"))
except Exception as e:
    print("OSD-411 ERR", e)
    results.append({"accession": "OSD-411", "status": f"error: {e}"})

# OSD-624: Col-0 roots, Flight vs Ground (sample names Flight-N-ROOTS-repN)
try:
    results.append(run_unnorm(
        "OSD-624",
        "OSD-624_transcription-profiling_rna-sequencing-(rna-seq)_Illumina",
        sample_rename=lambda s: re.sub(r"Flight-\d+", "FLT",
                                       re.sub(r"Ground-\d+", "GC", s, flags=re.I),
                                       flags=re.I)))
except Exception as e:
    print("OSD-624 ERR", e)
    results.append({"accession": "OSD-624", "status": f"error: {e}"})

pd.DataFrame(results).to_csv(f"{OUT}/replication_unnorm.csv", index=False)
print("done")
