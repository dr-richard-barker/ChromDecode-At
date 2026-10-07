"""OSDR-wide ChromDecode sweep: all Arabidopsis RNA-seq datasets.

For each accession: download GeneLab normalized counts, auto-detect the
primary contrast from metadata, run vectorized DE + state enrichment.
Fallback: expression-only decoding. Checkpoints per dataset.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import io
import json
import sys
import time
import urllib.request

sys.path.insert(0, f"{ROOT}/code")
import numpy as np
import pandas as pd
import chromdecode as cd

DATA = f"{ROOT}/data"
OUT = f"{ROOT}/results"
SWEEP = f"{OUT}/sweep"
API = "https://visualization.osdr.nasa.gov/biodata/api/v2"

import os
os.makedirs(SWEEP, exist_ok=True)


def fetch(url, timeout=120):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read()


def list_arabidopsis_rnaseq():
    """Return {accession: assay_name} for all Arabidopsis RNA-seq assays."""
    url = (f"{API}/query/metadata/?study.characteristics.organism="
           "Arabidopsis%20thaliana&format=csv")
    df = pd.read_csv(io.BytesIO(fetch(url, 180)), dtype=str)
    assays = {}
    for acc, assay in df[["id.accession", "id.assay name"]].dropna().values:
        if "rna" in assay.lower():
            assays.setdefault(acc, assay)
    return assays


def get_counts_url(acc, assay):
    url = f"{API}/dataset/{acc}/assay/{urllib.request.quote(assay, safe='()')}/files/?format=json"
    d = json.loads(fetch(url))
    files = d[acc]["assays"][assay]["files"]
    for f, info in files.items():
        if f.endswith("Normalized_Counts_GLbulkRNAseq.csv"):
            return info["URL"]
    return None


def get_metadata(acc):
    """Bulk per-sample metadata: characteristics + factor values."""
    url = (f"{API}/query/metadata/?id.accession={acc}"
           "&study.characteristics&study.factor.value&format=csv")
    return pd.read_csv(io.BytesIO(fetch(url, 180)), dtype=str)


def parse_condition_token(sample_name):
    """Extract a condition label from an OSDR sample name.

    Handles common plant-study conventions: FLT (spaceflight), GC (ground
    control), and gravity-level suffixes (0G microgravity, 1G control,
    0-08G partial gravity). Returns a label or None.
    """
    import re
    s = sample_name.upper()
    toks = re.split(r"[_\-]", s)
    for tok in reversed(toks):
        if tok in ("FLT", "SPACEFLIGHT", "SPACE"):
            return "Space Flight"
        if tok in ("GC", "GROUNDCONTROL", "GROUND"):
            return "Ground Control"
        if re.fullmatch(r"0G", tok):
            return "Microgravity (0G)"
        if re.fullmatch(r"1G", tok):
            return "1G Control"
        if re.fullmatch(r"0-\d+G", tok):
            return f"Partial Gravity ({tok})"
    return None


def sample_factors(meta, sample_names):
    """Build a samples x factors table from bulk metadata + sample names."""
    sample_col = next((c for c in meta.columns if "sample" in c.lower()
                       and "name" in c.lower()), None)
    if sample_col is None:
        return None
    meta = meta.set_index(sample_col)
    meta = meta[~meta.index.duplicated(keep="first")]
    factor_cols = [c for c in meta.columns
                   if "factor" in c.lower() or "characteristic" in c.lower()]
    factor_cols = [c for c in factor_cols
                   if meta[c].notna().any() and meta[c].nunique() > 1]
    rows = {}
    for s in sample_names:
        row = {}
        if s in meta.index:
            for c in factor_cols:
                row[c] = meta.loc[s, c]
        cond = parse_condition_token(s)
        if cond is not None:
            row["parsed_condition"] = cond
        rows[s] = row
    fac = pd.DataFrame.from_dict(rows, orient="index")
    keep = [c for c in fac.columns
            if fac[c].notna().any() and fac[c].nunique() > 1]
    fac = fac[keep]
    return fac if len(fac.columns) else None


def process_dataset(acc, assay, tss_state, log):
    t0 = time.time()
    res = {"accession": acc, "status": None, "n_de_up": 0, "n_de_down": 0,
           "contrast": None, "top_state_up": None, "top_fdr_up": None}
    try:
        # counts
        cpath = f"{SWEEP}/{acc}_counts.csv"
        if not os.path.exists(cpath):
            cu = get_counts_url(acc, assay)
            if cu is None:
                res["status"] = "no_processed_counts"
                return res
            data = fetch(cu, 600)
            open(cpath, "wb").write(data)
        counts = pd.read_csv(cpath, index_col=0)
        counts.index = [cd.clean_agi(i) or i for i in counts.index]
        counts = counts[~counts.index.str.startswith(("GLDS", "SRR"))]

        # metadata + contrast
        mpath = f"{SWEEP}/{acc}_metadata.csv"
        if not os.path.exists(mpath):
            meta = get_metadata(acc)
            meta.to_csv(mpath, index=False)
        else:
            meta = pd.read_csv(mpath, dtype=str)
        fac = sample_factors(meta, list(counts.columns))
        contrast = None
        if fac is not None:
            contrast = cd.detect_contrast_from_factors(fac)
        if contrast is None:
            # fallback: expression-only decoding
            enr = cd.expression_only_decode(np.log2(counts.clip(lower=1)),
                                            tss_state)
            enr.to_csv(f"{SWEEP}/{acc}_expr_only_enrichment.csv", index=False)
            res["status"] = "expression_only"
            res["contrast"] = "none-detected"
            res["top_state_up"] = enr.iloc[0]["st"]
            res["top_fdr_up"] = enr.iloc[0]["padj"]
            return res

        # DE with condition + covariates (log2 scale: linear models require
        # variance-stabilized values; raw normalized counts inflate DE)
        col, ctrl, treats, levels = (contrast["column"], contrast["control"],
                                     contrast["treatments"],
                                     contrast["levels"])
        cond = (levels.reindex(counts.columns) != ctrl).astype(float).values
        cov = {c: fac[c].reindex(counts.columns).values
               for c in fac.columns if c != col}
        # ecotype covariate from sample names (OSDR plant convention)
        import re as _re
        eco = [_re.search(r"(Col-0|Ler-0|Cvi-0|Ws-2)", s, _re.I).group(0).upper()
               if _re.search(r"(Col-0|Ler-0|Cvi-0|Ws-2)", s, _re.I) else None
               for s in counts.columns]
        if sum(e is not None for e in eco) > len(eco) * 0.8 and len(set(
                filter(None, eco))) > 1:
            cov["ecotype"] = eco
        de = cd.differential_expression_any(np.log2(counts.clip(lower=1)),
                                            cond, covariates=cov)
        de.to_csv(f"{SWEEP}/{acc}_DE.csv")
        n_up = int((de["direction"] == "up").sum())
        n_dn = int((de["direction"] == "down").sum())
        res.update({"status": "ok", "n_de_up": n_up, "n_de_down": n_dn,
                    "contrast": f"{col}: {ctrl} vs {treats[0]}"})

        if n_up + n_dn > 0:
            st = tss_state.dropna(subset=["state_group"])
            bg = st["agi"].unique()
            enr_up = cd.state_enrichment(st.rename(columns={"state_group": "st"}),
                                         de[de["direction"] == "up"].index,
                                         background=bg, state_col="st")
            enr_dn = cd.state_enrichment(st.rename(columns={"state_group": "st"}),
                                         de[de["direction"] == "down"].index,
                                         background=bg, state_col="st")
            enr_up.to_csv(f"{SWEEP}/{acc}_enrichment_up.csv", index=False)
            enr_dn.to_csv(f"{SWEEP}/{acc}_enrichment_down.csv", index=False)
            if n_up > 0:
                res["top_state_up"] = enr_up.iloc[0]["st"]
                res["top_fdr_up"] = enr_up.iloc[0]["padj"]
        else:
            res["status"] = "no_degs"
        return res
    except Exception as e:
        res["status"] = f"error: {type(e).__name__}: {e}"
        return res
    finally:
        res["seconds"] = round(time.time() - t0, 1)


def main():
    print("Enumerating Arabidopsis RNA-seq assays ...")
    assays = list_arabidopsis_rnaseq()
    print(f"  {len(assays)} accessions")
    print("Loading cached state assignments ...")
    tss = pd.read_csv(f"{OUT}/gene_state_assignments.csv")

    results = []
    for acc, assay in sorted(assays.items()):
        print(f"[{acc}] {assay}")
        r = process_dataset(acc, assay, tss, None)
        print(f"  -> {r['status']} ({r.get('seconds', '?')}s)")
        results.append(r)
        pd.DataFrame(results).to_csv(f"{SWEEP}/sweep_summary.csv", index=False)
    print("Sweep complete.")
    print(pd.DataFrame(results)[["accession", "status", "n_de_up",
                                 "n_de_down", "contrast"]].to_string(index=False))


if __name__ == "__main__":
    main()
