"""ChromDecode-At: Arabidopsis chromatin-state auto-decoder.

Takes an Arabidopsis gene-expression matrix or DEG table, infers associated
chromatin states from published epigenomic annotations (PCSD 36-state,
AraENCODE 12-state, PlantCADB ACRs), and produces human-readable figures.

Modules:
  1. ingest        - load expression/DEG tables, harmonize AGI IDs
  2. coordinates   - AGI -> TAIR10 coordinates
  3. states        - assign genes to chromatin states (interval intersection)
  4. enrichment    - state enrichment of gene sets vs background (Fisher + BH)
  5. learned       - elastic-net model: expression features -> state class
  6. figures       - SVG/PNG figure factory
"""

import numpy as np
import pandas as pd
from scipy import stats

CHROMS = [f"chr{i}" for i in range(1, 6)] + ["chrM", "chrC"]


# --------------------------------------------------------------------------
# Module 1: ingest
# --------------------------------------------------------------------------

def clean_agi(x):
    """Normalize an ID to primary AGI format (ATxGnnnnn)."""
    if not isinstance(x, str):
        return None
    x = x.strip().upper()
    for token in x.replace(";", ".").split("."):
        token = token.strip()
        if len(token) >= 9 and token.startswith("AT") \
                and token[2] in "12345MCG" and token[3] == "G" \
                and token[4:].isdigit():
            return token[:9]
    return None


def load_counts(path, kind="normalized"):
    """Load a GeneLab GLbulkRNAseq counts CSV -> DataFrame (genes x samples)."""
    df = pd.read_csv(path, index_col=0)
    df.index = [clean_agi(i) or i for i in df.index]
    return df


def differential_expression_vectorized(vst_counts, sample_names,
                                       padj_threshold=0.05,
                                       lfc_threshold=1.0,
                                       extra_covariates=None):
    """Vectorized linear-model DE on VST counts: condition + covariates.

    sample_names: list matching columns, parsed for FLT/GC condition and
    ecotype. extra_covariates: optional dict {name: array of level labels}.
    Returns DataFrame indexed by AGI with log2fc, p, padj, direction.
    """
    cond, eco = [], []
    for s in sample_names:
        s_up = s.upper()
        cond.append("FLT" if "_FLT_" in s_up or "FLT" in s_up.split("_")
                    else "GC")
        e = next((x for x in ["COL-0", "LER-0", "CVI-0", "WS-2"]
                  if x in s_up), "UNKNOWN")
        eco.append(e)
    cond = np.array(cond)
    eco = np.array(eco)

    parts = [pd.get_dummies(eco, drop_first=True).astype(float)]
    if extra_covariates:
        for name, levels in extra_covariates.items():
            if len(set(levels)) > 1 and len(set(levels)) <= 5:
                parts.append(pd.get_dummies(
                    pd.Series(levels, dtype=str), prefix=name,
                    drop_first=True).astype(float))
    design = pd.concat(parts, axis=1)
    design["condition"] = (cond == "FLT").astype(float)
    design = sm_add_constant(design)
    X = design.values.astype(float)
    n, p = X.shape
    dof = n - p
    if dof < 1:
        raise ValueError(f"Not enough samples: n={n}, params={p}")

    Y = vst_counts.values.T.astype(float)          # samples x genes
    betas, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)   # p x genes
    resid = Y - X @ betas
    sigma2 = (resid ** 2).sum(axis=0) / dof        # genes
    XtX_inv = np.linalg.pinv(X.T @ X)
    se = np.sqrt(sigma2 * XtX_inv[-1, -1])
    tvals = betas[-1] / np.maximum(se, 1e-12)
    pvals = 2 * stats.t.sf(np.abs(tvals), dof)

    de = pd.DataFrame({"log2fc": betas[-1], "p_value": pvals},
                      index=vst_counts.index)
    de["padj"] = bh_adjust(de["p_value"].values)
    de["direction"] = np.where(
        (de["padj"] < padj_threshold) & (de["log2fc"] > lfc_threshold), "up",
        np.where((de["padj"] < padj_threshold) & (de["log2fc"] < -lfc_threshold),
                 "down", "ns"))
    return de


def sm_add_constant(df):
    import statsmodels.api as sm
    return sm.add_constant(df, has_constant="add")


def differential_expression(vst_counts, sample_names, padj_threshold=0.05,
                            lfc_threshold=1.0):
    """Linear-model DE on VST counts: condition (FLT vs GC) + ecotype.

    sample_names: list matching columns, parsed for '_FLT_' / '_GC_' and
    ecotype token (Col-0, Ler-0, Cvi-0, Ws-2).
    Returns DataFrame indexed by AGI with log2fc, p, padj.
    """
    import statsmodels.api as sm

    cond, eco = [], []
    for s in sample_names:
        s_up = s.upper()
        cond.append("FLT" if "_FLT_" in s_up or "FLT" in s_up.split("_") else "GC")
        e = next((x for x in ["COL-0", "LER-0", "CVI-0", "WS-2"] if x in s_up),
                 "UNKNOWN")
        eco.append(e)
    cond = np.array(cond)
    eco = np.array(eco)
    assert set(cond) == {"FLT", "GC"}, f"condition parse failed: {set(cond)}"

    design = pd.get_dummies(eco, drop_first=True).astype(float)
    design["flt"] = (cond == "FLT").astype(float)
    design = sm.add_constant(design)
    design = design.astype(float)

    X = design.values
    n, p = X.shape
    dof = n - p
    res = []
    Y = vst_counts.values.T.astype(float)  # samples x genes
    for j, gene in enumerate(vst_counts.index):
        y = Y[:, j]
        beta, rss, _, _ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        sigma2 = resid @ resid / dof
        cov = sigma2 * np.linalg.pinv(X.T @ X)
        se = np.sqrt(np.diag(cov))
        tval = beta / se
        pval = 2 * stats.t.sf(np.abs(tval), dof)
        res.append((gene, beta[-1], pval[-1]))
    de = pd.DataFrame(res, columns=["gene", "log2fc", "p_value"]).set_index("gene")
    de["padj"] = bh_adjust(de["p_value"].values)
    de["direction"] = np.where(
        (de["padj"] < padj_threshold) & (de["log2fc"] > lfc_threshold), "up",
        np.where((de["padj"] < padj_threshold) & (de["log2fc"] < -lfc_threshold),
                 "down", "ns"))
    return de


def detect_contrast_from_factors(fac):
    """Auto-detect primary contrast from a factor DataFrame (samples x factors).

    Prefers a factor with a control-like level. Returns dict(column, control,
    treatments, levels Series indexed by sample) or None.
    """
    control_tokens = ["ground control", "mock", "control", "untreated",
                      "0 gy", "0 gray", "no treatment", "baseline"]
    best = None
    for c in fac.columns:
        vals = fac[c].dropna().astype(str)
        if vals.nunique() < 2 or vals.nunique() > 6:
            continue
        ctrl_level = next((v for v in vals.unique()
                           if any(t in v.lower() for t in control_tokens)),
                          None)
        if ctrl_level is None:
            continue
        treats = [v for v in vals.unique() if v != ctrl_level]
        score = (len(vals),)
        if best is None or score > best[0]:
            best = (score, c, ctrl_level, treats, vals)
    if best is None:
        return None
    _, col, ctrl, treats, vals = best
    return {"column": col, "control": ctrl, "treatments": treats,
            "levels": vals}


def differential_expression_any(counts, condition, covariates=None,
                                padj_threshold=0.05, lfc_threshold=1.0):
    """Vectorized DE: binary condition (1=treatment) + optional covariates.

    counts: genes x samples DataFrame. condition: array-like 0/1 per sample.
    covariates: dict {name: array of level labels per sample}; factors with
    2-5 levels enter as dummies.
    """
    parts = []
    if covariates:
        for name, levels in covariates.items():
            lv = pd.Series(levels, dtype=str)
            if lv.nunique() > 1 and lv.nunique() <= 5:
                parts.append(pd.get_dummies(lv, prefix=name,
                                            drop_first=True).astype(float))
    design = pd.concat(parts, axis=1) if parts else pd.DataFrame(
        index=counts.columns)
    design["condition"] = np.asarray(condition, dtype=float)
    design = sm_add_constant(design)
    X = design.values.astype(float)
    n, p = X.shape
    dof = n - p
    if dof < 1:
        raise ValueError(f"Not enough samples: n={n}, params={p}")

    Y = counts.values.T.astype(float)
    betas, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ betas
    sigma2 = (resid ** 2).sum(axis=0) / dof
    XtX_inv = np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.maximum(sigma2 * XtX_inv[-1, -1], 1e-12))
    tvals = betas[-1] / se
    pvals = 2 * stats.t.sf(np.abs(tvals), dof)

    de = pd.DataFrame({"log2fc": betas[-1], "p_value": pvals},
                      index=counts.index)
    de["padj"] = bh_adjust(de["p_value"].values)
    de["direction"] = np.where(
        (de["padj"] < padj_threshold) & (de["log2fc"] > lfc_threshold), "up",
        np.where((de["padj"] < padj_threshold) & (de["log2fc"] < -lfc_threshold),
                 "down", "ns"))
    return de


def detect_contrast(meta_csv_path):
    """Auto-detect the primary contrast from an OSDR metadata CSV.

    Looks for factor/characteristic columns; prefers a column with a
    control-like level (ground control, mock, 0 dose, control, untreated).
    Returns (column_name, control_level, treatment_level, levels_series)
    or None.
    """
    meta = pd.read_csv(meta_csv_path, dtype=str)
    meta.columns = [c.strip() for c in meta.columns]
    cand_cols = [c for c in meta.columns
                 if "factor" in c.lower() or "characteristic" in c.lower()
                 or "treatment" in c.lower() or "condition" in c.lower()]
    control_tokens = ["ground control", "mock", "control", "untreated",
                      "0 gy", "0 gray", "no treatment", "baseline"]
    best = None
    for c in cand_cols:
        vals = meta[c].dropna().astype(str)
        if vals.nunique() < 2 or vals.nunique() > 6:
            continue
        low = vals.str.lower()
        ctrl_level = next((v for v in vals.unique()
                           if any(t in v.lower() for t in control_tokens)),
                          None)
        if ctrl_level:
            treat_levels = [v for v in vals.unique() if v != ctrl_level]
            if len(treat_levels) >= 1:
                score = (0, len(vals))  # prefer control-bearing, larger n
                if best is None or score > best[0]:
                    best = (score, c, ctrl_level, treat_levels, vals)
    if best is None:
        return None
    _, col, ctrl, treats, vals = best
    return {"column": col, "control": ctrl, "treatments": treats,
            "levels": vals}


def expression_only_decode(counts, tss_state, top_frac=0.25):
    """Fallback decoding: state composition of top-expressed genes vs background."""
    mean_expr = counts.mean(axis=1)
    n_top = max(1, int(len(mean_expr) * top_frac))
    top_genes = mean_expr.sort_values(ascending=False).head(n_top).index
    st = tss_state.dropna(subset=["state_group"])
    enr = state_enrichment(st.rename(columns={"state_group": "st"}),
                           top_genes, state_col="st")
    return enr


def bh_adjust(pvals):
    """Benjamini-Hochberg FDR adjustment."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


# --------------------------------------------------------------------------
# Module 2: coordinates
# --------------------------------------------------------------------------

def load_gene_coordinates(path):
    """Load AraENCODE genes_TAIR10.txt -> DataFrame (primary AGI, one row/gene)."""
    g = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 3, 4],
                    names=["id", "chr", "start", "end", "strand"],
                    on_bad_lines="skip")
    g["agi"] = g["id"].map(clean_agi)
    g = g.dropna(subset=["agi"]).drop_duplicates("agi", keep="first")
    g = g[g["chr"].isin(CHROMS)]
    return g[["agi", "chr", "start", "end", "strand"]].reset_index(drop=True)


def gene_regions(genes, upstream=2000):
    """Gene body + upstream regulatory region per gene."""
    reg = genes.copy()
    plus = reg["strand"] == "+"
    reg["tss"] = np.where(plus, reg["start"], reg["end"])
    reg["rstart"] = np.where(plus, reg["start"] - upstream, reg["start"])
    reg["rend"] = np.where(plus, reg["end"], reg["end"] + upstream)
    reg["rstart"] = reg["rstart"].clip(lower=0)
    reg["gene_end"] = reg["end"]
    return reg[["agi", "chr", "rstart", "rend", "tss", "gene_end"]].rename(
        columns={"rstart": "start", "rend": "end"})


# --------------------------------------------------------------------------
# Module 3: state assignment
# --------------------------------------------------------------------------

def load_pcsd_states(data_dir, state_desc_path):
    """Load PCSD At_genes_S1..S36 -> gene->state DataFrame + descriptions."""
    rows = []
    for s in range(1, 37):
        g = pd.read_csv(f"{data_dir}/At_genes_S{s}", sep="\t", header=None,
                        usecols=[1], names=["agi"])
        g["agi"] = g["agi"].map(clean_agi)
        g = g.dropna()
        rows.append(pd.DataFrame({"agi": g["agi"].unique(), "pcsd_state": s}))
    gene_states = pd.concat(rows, ignore_index=True)
    # genes can appear in multiple states (segment overlap); keep all, dedupe later
    desc = pd.read_csv(state_desc_path, sep="\t")
    return gene_states, desc


def load_pcsd_segments(data_dir):
    """Load PCSD At_segments_S1..S36 -> DataFrame(chr,start,end,state)."""
    frames = []
    for s in range(1, 37):
        seg = pd.read_csv(f"{data_dir}/At_segments_S{s}", sep="\t", header=None,
                          names=["chr", "start", "end", "label"])
        seg["pcsd_state"] = s
        frames.append(seg)
    return pd.concat(frames, ignore_index=True)


def load_araencode_states(bedgz_path, tbl_path):
    """Load AraENCODE 12-state track -> DataFrame(chr,start,end,state_id,name)."""
    seg = pd.read_csv(bedgz_path, sep="\t", header=None, compression="gzip",
                      names=["chr", "start", "end", "state_id"])
    names = pd.read_csv(tbl_path, sep="\t", header=None,
                        names=["chr", "start", "end", "state_id", "state_name"])
    name_map = names.drop_duplicates("state_id").set_index("state_id")["state_name"]
    seg["state_name"] = seg["state_id"].map(name_map)
    return seg


def assign_state_by_tss(regions, segments, value_col="state_name",
                        out_col="state"):
    """Assign each gene the state of the segment covering its TSS."""
    seg_by_chr = {c: d.sort_values("start").reset_index(drop=True)
                  for c, d in segments.groupby("chr")}
    starts = {c: d["start"].values for c, d in seg_by_chr.items()}
    ends = {c: d["end"].values for c, d in seg_by_chr.items()}
    vals = {c: d[value_col].values for c, d in seg_by_chr.items()}

    out = []
    for chr_, tss in zip(regions["chr"], regions["tss"]):
        d = seg_by_chr.get(chr_)
        if d is None:
            out.append(None)
            continue
        i = np.searchsorted(starts[chr_], tss, side="right") - 1
        if i >= 0 and tss < ends[chr_][i]:
            out.append(vals[chr_][i])
        else:
            out.append(None)
    res = regions.copy()
    res[out_col] = out
    return res


def assign_acr_overlap(regions, acr_bed_path):
    """Flag genes whose regulatory region overlaps a PlantCADB ACR.

    Streams the whole-genome associated-genes file line by line (it mixes 37
    species and has ragged row lengths); Arabidopsis ACRs are identified by
    AGI IDs (ATxGnnnnn) occurring on chr1-5 rows.
    """
    import re
    import gzip as gz
    pat = re.compile(rb"AT[1-5MCG]G\d{5}")
    acr_agi = set()
    with gz.open(acr_bed_path, "rb") as fh:
        for line in fh:
            parts = line.split(b"\t", 1)[1]  # skip chr field for speed check
            if not line.startswith((b"chr1\t", b"chr2\t", b"chr3\t",
                                    b"chr4\t", b"chr5\t")):
                continue
            for m in pat.findall(line):
                acr_agi.add(m.decode())
    res = regions.copy()
    res["in_acr"] = res["agi"].isin(acr_agi)
    return res


# --------------------------------------------------------------------------
# Module 4: enrichment
# --------------------------------------------------------------------------

def state_enrichment(gene_states, gene_set, background=None, state_col="pcsd_state"):
    """Fisher exact enrichment of a gene set across states.

    gene_states: DataFrame(agi, state_col) - one row per gene (deduped).
    gene_set: iterable of AGI IDs (e.g. up-regulated DEGs).
    background: iterable; default all genes in gene_states.
    Returns DataFrame(state, n_set, n_bg, frac_set, frac_bg, odds_ratio, p, padj).
    """
    gs = gene_states.drop_duplicates("agi")
    if background is None:
        background = gs["agi"]
    bg = set(background)
    target = set(gene_set) & bg
    n_bg = len(bg)
    n_target = len(target)

    rows = []
    for st, grp in gs.groupby(state_col):
        st_genes = set(grp["agi"])
        a = len(st_genes & target)          # in set & in state
        b = n_target - a                    # in set & not state
        c = len((st_genes & bg) - target)   # not set & in state
        d = n_bg - n_target - c             # not set & not state
        table = [[a, b], [c, d]]
        odds, p = stats.fisher_exact(table, alternative="two-sided")
        rows.append((st, a, c, a / max(n_target, 1), c / max(n_bg, 1), odds, p))
    out = pd.DataFrame(rows, columns=[state_col, "n_in_set", "n_in_bg",
                                      "frac_in_set", "frac_in_bg",
                                      "odds_ratio", "p_value"])
    out["padj"] = bh_adjust(out["p_value"].values)
    return out.sort_values("p_value").reset_index(drop=True)
