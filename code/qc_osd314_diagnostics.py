"""Re-run the OSD-314 diagnostics and GO replication on the corrected contrasts.

The published diagnostics (diagnose_osd314.py -> Tables S13-S17) and GO
replication (motif_go_enrichment.py step 6b -> Table S21) used the sweep's
OSD-314 contrast, which pooled the 0.3g (Mars) arm with 0g. This script repeats
the same tests, with the same code paths and the same seed, for three contrasts:

  published   the sweep's pooled contrast (0g + 0.3g vs 1g; Table S57) --
              a reproduction check against S13-S17 and S21
  0g_light    0g vs 1g, light as covariate (Table S55)
  03g_light   0.3g vs 1g, light as covariate (Table S56)

The 0.3g contrast is the primary corrected one: its 1g controls come from the
same centrifuge run, whereas the 0g samples may be a separate batch.

Inputs: data/derived/gene_state_assignments.csv, results/sweep/OSD-314_counts.csv,
results/qc/osd314_sample_design.csv (both from qc_sweep_rerun.py),
results/sweep/OSD-37_DE.csv, and for GO data/motif_go/{TAIR10.fa.gz,
tair.gaf.gz, go.obo}. Outputs go to results/qc/.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import gzip
import sys
from collections import defaultdict

sys.path.insert(0, f"{ROOT}/code")
import numpy as np
import pandas as pd
from scipy import stats as sps
import chromdecode as cd

QC = f"{ROOT}/results/qc"
SWEEP = f"{ROOT}/results/sweep"
MOTIF = f"{ROOT}/data/motif_go"
TABLES = f"{ROOT}/supplementary_tables"

tss = pd.read_csv(f"{ROOT}/data/derived/gene_state_assignments.csv")
st = tss.dropna(subset=["pcsd_state"]).rename(columns={"pcsd_state": "st"})[["agi", "st"]]
st["st"] = st["st"].astype(int)
state_of = st.set_index("agi")["st"]
bg = st["agi"].unique()
POLY_DIAG = [11, 12, 13, 15]          # as in diagnose_osd314.py
POLY_GO = set(range(11, 16))          # as in motif_go_enrichment.py

de37 = pd.read_csv(f"{SWEEP}/OSD-37_DE.csv", index_col=0)
counts = pd.read_csv(f"{SWEEP}/OSD-314_counts.csv", index_col=0)
counts.index = [cd.clean_agi(i) or i for i in counts.index]
counts = counts[~counts.index.str.startswith(("GLDS", "SRR"))]
design = pd.read_csv(f"{QC}/osd314_sample_design.csv", index_col=0)
mat_all = np.log2(counts.clip(lower=1))

CONTRASTS = {
    "published": dict(de=f"{SWEEP}/OSD-314_DE.csv", samples=list(design.index),
                      cond=design["sweep_coded_as_treatment"].astype(float), light=False),
    "0g_light": dict(de=f"{QC}/osd314_0g_vs_1g_light_DE.csv",
                     samples=list(design.index[design.gravity.isin(["0g", "1g"])]),
                     cond=(design.gravity == "0g").astype(float), light=True),
    "03g_light": dict(de=f"{QC}/osd314_03g_vs_1g_light_DE.csv",
                      samples=list(design.index[design.gravity.isin(["0.3g", "1g"])]),
                      cond=(design.gravity == "0.3g").astype(float), light=True),
}


def odds_for(gene_set, state):
    genes_in = set(state_of[state_of == state].index)
    a = len(set(gene_set) & genes_in)
    b = len(gene_set) - a
    c = len(genes_in) - a
    d = len(bg) - len(gene_set) - c
    o, p = sps.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return o, p, a


def diagnostics(tag, spec):
    """Tests 1-5 of diagnose_osd314.py, same order and seed."""
    rng = np.random.RandomState(42)
    de = pd.read_csv(spec["de"], index_col=0)
    up = de[de["direction"] == "up"].index
    rows = {"contrast": tag, "n_up": len(up)}

    # 1. null calibration
    null = {s: np.zeros(10_000) for s in POLY_DIAG}
    for i in range(10_000):
        draw = rng.choice(bg, size=len(up), replace=False)
        for s in POLY_DIAG:
            null[s][i] = odds_for(draw, s)[0]
    t1 = []
    for s in POLY_DIAG:
        o, p, a = odds_for(up, s)
        t1.append({"contrast": tag, "state": f"S{s}", "observed_odds": o, "fisher_p": p,
                   "empirical_p": (null[s] >= o).mean(), "n_genes_in_state": a,
                   "null_mean": null[s].mean(), "null_p99": np.quantile(null[s], 0.99)})

    # 2. OSD-37 downsampling control: only defined when n_up <= OSD-37's 529 up genes
    up37 = de37[de37["direction"] == "up"].index
    t2 = []
    feasible = len(up) <= len(up37)
    ds = {s: np.zeros(1000) for s in POLY_DIAG}
    if feasible:
        for i in range(1000):
            sub = rng.choice(up37, size=len(up), replace=False)
            for s in POLY_DIAG:
                ds[s][i] = odds_for(sub, s)[0]
    for r in t1:
        s = int(r["state"][1:])
        t2.append({"contrast": tag, "state": r["state"], "osd314_observed": r["observed_odds"],
                   "osd37_downsampled_mean": ds[s].mean() if feasible else np.nan,
                   "osd37_downsampled_p99": np.quantile(ds[s], 0.99) if feasible else np.nan,
                   "frac_draws_ge_observed": (ds[s] >= r["observed_odds"]).mean() if feasible else np.nan,
                   "note": "" if feasible else f"not applicable: n_up={len(up)} > {len(up37)} OSD-37 up genes"})

    # 3. threshold grid on the same design as the contrast
    samples = spec["samples"]
    cond = spec["cond"].reindex(samples).values
    cov = {"light": design.loc[samples, "light"].values} if spec["light"] else None
    t3 = []
    for padj_t in [0.01, 0.05, 0.1]:
        for lfc_t in [0.5, 1.0, 2.0]:
            g = cd.differential_expression_any(mat_all[samples], cond, covariates=cov,
                                               padj_threshold=padj_t, lfc_threshold=lfc_t)
            u = g[g["direction"] == "up"].index
            row = {"contrast": tag, "padj": padj_t, "lfc": lfc_t, "n_up": len(u)}
            for s in POLY_DIAG:
                row[f"S{s}"] = odds_for(u, s)[0] if len(u) else np.nan
            t3.append(row)

    # 4. threshold-free running enrichment (both tails reported)
    poly_genes = set(state_of[state_of.isin(POLY_DIAG)].index)
    ranked = de["log2fc"].sort_values(ascending=False)
    hit = ranked.index.isin(poly_genes)
    n_hits, n_miss = hit.sum(), len(hit) - hit.sum()
    run = np.cumsum(np.where(hit, 1 / n_hits, -1 / n_miss))
    perm = np.zeros(1000)
    for i in range(1000):
        prun = np.cumsum(np.where(rng.permutation(hit), 1 / n_hits, -1 / n_miss))
        perm[i] = prun.max() if abs(prun).max() == prun.max() else prun.min()
    t4 = [{"contrast": tag, "n_poly_genes": int(n_hits), "max_running_sum": run.max(),
           "min_running_sum": run.min(), "perm_p": (perm >= run.max()).mean(),
           "n_perm": 1000}]

    # 5. gene identity
    t5 = [{"contrast": tag, "agi": g, "state": f"S{state_of[g]}",
           "log2fc": de.loc[g, "log2fc"], "padj": de.loc[g, "padj"]}
          for g in up if g in state_of.index and state_of[g] in POLY_DIAG]
    rows.update({"n_up_in_polycomb": len(t5)})
    return t1, t2, t3, t4, t5, rows


all_t = [[], [], [], [], []]
summary = []
for tag, spec in CONTRASTS.items():
    print(f"== {tag}")
    out = diagnostics(tag, spec)
    for k in range(5):
        all_t[k] += out[k]
    summary.append(out[5])
    print(pd.DataFrame(out[0])[["state", "observed_odds", "empirical_p", "n_genes_in_state"]].to_string(index=False))
names = ["null_calibration", "downsampling_control", "threshold_grid", "gsea", "up_polycomb_genes"]
for k, n in enumerate(names):
    pd.DataFrame(all_t[k]).to_csv(f"{QC}/osd314_{n}_corrected.csv", index=False)

# reproduction check against the published tables (S13-S16)
pub = pd.DataFrame(all_t[0]).query("contrast == 'published'").set_index("state")
s13 = pd.read_csv(f"{TABLES}/osd314_null_calibration.csv", index_col=0)
s16 = pd.read_csv(f"{TABLES}/osd314_gsea.csv")
s15 = pd.read_csv(f"{TABLES}/osd314_threshold_grid.csv")
grid_pub = pd.DataFrame(all_t[2]).query("contrast == 'published'").reset_index(drop=True)
chk = {
    "S13_max_abs_diff_observed_odds": float((pub["observed_odds"] - s13["observed_odds"]).abs().max()),
    "S13_max_abs_diff_empirical_p": float((pub["empirical_p"] - s13["empirical_p"]).abs().max()),
    "S15_max_abs_diff_odds": float(np.nanmax(np.abs(grid_pub[["S11", "S12", "S13", "S15"]].values
                                                    - s15[["S11", "S12", "S13", "S15"]].values))),
    "S15_n_up_match": bool((grid_pub["n_up"].values == s15["n_up"].values).all()),
    "S16_max_running_sum_diff": float(abs(pd.DataFrame(all_t[3]).iloc[0]["max_running_sum"]
                                          - s16["max_running_sum"].iloc[0])),
}
print("reproduction check:", chk)

# ---------------------------------------------------------------- GO replication
print("== GO replication")
seqs, cur, buf = {}, None, []
with gzip.open(f"{MOTIF}/TAIR10.fa.gz", "rt") as f:
    for line in f:
        if line.startswith(">"):
            if cur:
                seqs[cur] = "".join(buf)
            cur, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    seqs[cur] = "".join(buf)
seqs = {("chr" + k[3].lower() if k.startswith("Chr") else ("chr" + k if k.isdigit() else k.lower())): v
        for k, v in seqs.items()}
# gene coordinates rebuilt from D1 (region = gene body + 2 kb upstream)
g = tss.copy()
plus = g["gene_end"] == g["end"]
g["strand"] = np.where(plus, "+", "-")
g["gstart"] = np.where(plus, g["tss"], g["start"])
g["gend"] = np.where(plus, g["end"], g["tss"])
comp = str.maketrans("ACGTN", "TGCAN")
promoters = set()
for r in g.itertuples():
    s = seqs.get(r.chr)
    if s is None:
        continue
    if r.strand == "+":
        seq = s[max(0, r.gstart - 1 - 2000):r.gstart - 1]
    else:
        seq = s[r.gend:min(len(s), r.gend + 2000)]
    if len(seq) >= 200:
        promoters.add(r.agi)
print(f"  promoters: {len(promoters):,}")

KEEP_EV = {"IDA", "IMP", "IGI", "IEP", "HDA", "HMP", "HGI", "HEP", "TAS", "IC", "EXP", "IEA"}
gene_by_go = defaultdict(set)
with gzip.open(f"{MOTIF}/tair.gaf.gz", "rt") as gf:
    for line in gf:
        if line.startswith("!"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 15 or p[6] not in KEEP_EV:
            continue
        agi = cd.clean_agi(p[1])
        if agi:
            gene_by_go[p[4]].add(agi)
term_name, cur_term = {}, None
for line in open(f"{MOTIF}/go.obo"):
    line = line.strip()
    if line == "[Term]":
        cur_term = "pending"
    elif cur_term == "pending" and line.startswith("id: "):
        cur_term = line[4:]
    elif cur_term and cur_term != "pending" and line.startswith("name: "):
        term_name[cur_term] = line[6:]

state_any = tss.set_index("agi")["pcsd_state"]


def go_enrichment(target_set, background_set, min_genes=5):
    rows, bgs = [], set(background_set)
    for term, gset in gene_by_go.items():
        a = len(target_set & gset)
        if a < min_genes:
            continue
        c = len(bgs & gset)
        _, pv = sps.fisher_exact([[a, len(target_set) - a], [c, len(bgs) - c]], alternative="greater")
        rows.append({"GO": term, "name": term_name.get(term, ""), "n_in_set": a, "n_in_bg": c,
                     "fold": (a / len(target_set)) / (c / len(bgs)), "p_value": pv})
    e = pd.DataFrame(rows)
    e["padj"] = cd.bh_adjust(e["p_value"].values)
    return e.sort_values("p_value")


down37 = set(de37[de37.direction == "down"].index)
target37 = {x for x in down37 if x in state_any.index and state_any[x] in POLY_GO and x in promoters}
uni37 = promoters & set(de37.index)
go_uni = go_enrichment(target37, uni37)
s19 = pd.read_csv(f"{TABLES}/go_enrichment.csv").query("comparison == 'target_vs_universe'")
top37 = go_uni[go_uni.padj < 0.05].head(10)
chk["S19_target_n_now_vs_published"] = f"{len(target37)} vs 265"
chk["S19_n_sig_terms_now_vs_published"] = f"{int((go_uni.padj < 0.05).sum())} vs {int((s19.padj < 0.05).sum())}"
chk["S19_top10_terms_shared_with_published"] = len(set(top37.GO) & set(s19[s19.padj < 0.05].head(10).GO))
print("  GO check:", {k: v for k, v in chk.items() if k.startswith("S19")})

rep_rows = []
for tag, spec in CONTRASTS.items():
    de = pd.read_csv(spec["de"], index_col=0)
    t314 = {x for x in de[de.direction == "down"].index
            if x in state_any.index and state_any[x] in POLY_GO and x in promoters}
    u314 = promoters & set(de.index)
    for _, r in top37.iterrows():
        gset = gene_by_go[r["GO"]]
        a, c = len(t314 & gset), len(u314 & gset)
        _, pv = sps.fisher_exact([[a, len(t314) - a], [c, len(u314) - c]], alternative="greater")
        rep_rows.append({"contrast": tag, "GO": r["GO"], "name": r["name"],
                         "OSD37_fold": r["fold"], "OSD37_padj": r["padj"],
                         "OSD314_polycomb_down_n": len(t314), "OSD314_n": a,
                         "OSD314_fold": (a / max(len(t314), 1)) / (c / len(u314)) if c else np.nan,
                         "OSD314_p": pv})
rep = pd.DataFrame(rep_rows)
rep["OSD314_padj_within_contrast"] = rep.groupby("contrast")["OSD314_p"].transform(
    lambda x: cd.bh_adjust(x.values))
rep.to_csv(f"{QC}/osd314_go_replication_corrected.csv", index=False)
print(rep.pivot_table(index="name", columns="contrast", values="OSD314_fold").round(2).to_string())

s21 = pd.read_csv(f"{TABLES}/osd314_go_replication.csv")
pub_rep = rep.query("contrast == 'published'").set_index("GO")
shared = s21.set_index("GO").index.intersection(pub_rep.index)
chk["S21_terms_shared"] = len(shared)
chk["S21_max_abs_diff_OSD314_n"] = int((pub_rep.loc[shared, "OSD314_n"]
                                        - s21.set_index("GO").loc[shared, "OSD314_n"]).abs().max()) if len(shared) else None
pd.Series(chk).to_csv(f"{QC}/osd314_diagnostics_reproduction_check.csv", header=["value"])
print("final check:", chk)
