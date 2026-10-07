"""Motif + GO enrichment of Polycomb-repressed down-regulated genes.

Gene sets (OSD-37 primary):
  target   : down DEGs in PCSD Polycomb states S11-S15
  control1 : down DEGs in non-Polycomb states
  control2 : expressed, non-DE genes in S11-S15 (state-matched background)
Motifs: JASPAR plant PFMs -> log-odds PWMs, scan 2 kb promoters
(strand-aware), gene = target if any position p < 1e-4.
GO: hypergeometric ORA with expressed-gene universe, BH-FDR.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import gzip
import os
import sys
from collections import defaultdict

sys.path.insert(0, f"{ROOT}/code")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
matplotlib.rcParams["svg.fonttype"] = "none"

import chromdecode as cd

DATA = f"{ROOT}/data"
MOTIF = f"{DATA}/motif_go"
OUT = f"{ROOT}/results"
SWEEP = f"{OUT}/sweep"
FIG = f"{ROOT}/figures"
UPSTREAM = 2000
P_THRESHOLD = 1e-4
CODE = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3,
        ord("N"): 4, 10: 4, 13: 4}

# ---------------------------------------------------------------- promoters
print("[1] Extracting strand-aware 2 kb promoters ...")
seqs = {}
cur, buf = None, []
with gzip.open(f"{MOTIF}/TAIR10.fa.gz", "rt") as f:
    for line in f:
        if line.startswith(">"):
            if cur:
                seqs[cur] = "".join(buf)
            cur = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
# normalize chromosome names to chr1..chr5, chrM, chrC
seqs = {("chr" + k[3].lower() if k.startswith("Chr")
         else ("chr" + k if k.isdigit() else k.lower())): v
        for k, v in seqs.items()}
print(f"  genome: {len(seqs)} sequences, "
      f"{sum(len(v) for v in seqs.values()):,} bp")

genes = cd.load_gene_coordinates(f"{DATA}/araencode/genes_TAIR10.txt")
comp = str.maketrans("ACGTN", "TGCAN")
promoters = {}
for r in genes.itertuples():
    s = seqs.get(r.chr)
    if s is None:
        continue
    if r.strand == "+":
        start = max(0, r.start - 1 - UPSTREAM)
        seq = s[start:r.start - 1]
    else:
        end = min(len(s), r.end + UPSTREAM)
        seq = s[r.end:end][::-1].translate(comp)
    if len(seq) >= 200:
        promoters[r.agi] = seq.upper()
print(f"  promoters extracted: {len(promoters):,}")

# ---------------------------------------------------------------- gene sets
print("[2] Defining gene sets ...")
tss = pd.read_csv(f"{OUT}/gene_state_assignments.csv")
state_of = tss.set_index("agi")["pcsd_state"]
POLY = set(range(11, 16))

de37 = pd.read_csv(f"{SWEEP}/OSD-37_DE.csv", index_col=0)
down37 = set(de37[de37["direction"] == "down"].index)
ns37 = set(de37[de37["direction"] == "ns"].index)

target = {g for g in down37 if g in state_of.index
          and state_of[g] in POLY and g in promoters}
control1 = {g for g in down37 if g in state_of.index
            and state_of[g] not in POLY and g in promoters}
control2 = {g for g in ns37 if g in state_of.index
            and state_of[g] in POLY and g in promoters}
universe = set(promoters) & set(de37.index)
print(f"  target (Polycomb-down): {len(target)}")
print(f"  control1 (non-Polycomb-down): {len(control1)}")
print(f"  control2 (Polycomb-non-DE): {len(control2)}")
print(f"  universe (expressed, promoter available): {len(universe)}")

def gc_frac(seq):
    s = seq.replace("N", "")
    return (s.count("G") + s.count("C")) / len(s)

gc_stats = {k: float(np.mean([gc_frac(promoters[g]) for g in v]))
            for k, v in [("target", target), ("control1", control1),
                         ("control2", control2), ("universe", universe)]}
print("  mean GC:", {k: round(v, 3) for k, v in gc_stats.items()})
pd.Series(gc_stats).to_csv(f"{OUT}/motif_gc_content.csv")

# ------------------------------------------------------- concatenated codes
print("[3] Building scan arrays ...")
gene_list = sorted(promoters)
offsets = np.zeros(len(gene_list) + 1, dtype=np.int64)
chunks, rchunks = [], []
TABLE = np.full(256, 4, dtype=np.int8)
for b, v in [(65, 0), (67, 1), (71, 2), (84, 3)]:
    TABLE[b] = v  # A C G T; everything else (N, IUPAC, separators) -> 4
for i, g in enumerate(gene_list):
    s = promoters[g]
    offsets[i] = sum(len(c) for c in chunks)
    fwd = TABLE[np.frombuffer(s.encode(), dtype=np.uint8)]
    rev = fwd[::-1].copy()
    m = (rev >= 0) & (rev <= 3)
    rev[m] = 3 - rev[m]
    chunks.append(fwd)
    rchunks.append(rev)
offsets[len(gene_list)] = sum(len(c) for c in chunks)
codes = np.concatenate(chunks)
rcodes = np.concatenate(rchunks)
N = len(codes)
print(f"  concatenated scan length: {N:,} bp x 2 strands")

# ---------------------------------------------------------------- PWM scan
import pickle
CKPT = f"{MOTIF}/motif_targets.pkl"
SCAN_DONE = os.path.exists(CKPT)
def load_pfm(path):
    rows, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            parts = line.replace("[", " ").replace("]", " ").split()
            if len(parts) >= 2 and parts[0] in "ACGT":
                rows[parts[0]] = np.array([float(x) for x in parts[1:]])
    if len(rows) != 4 or len({len(v) for v in rows.values()}) != 1:
        return None, None
    return name, np.vstack([rows[b] for b in "ACGT"])

rng = np.random.RandomState(42)
all_seq = "".join(promoters[g][:500] for g in gene_list[:2000])
bg_freqs = np.array([all_seq.count(b) for b in "ACGT"], dtype=float)
bg_freqs /= bg_freqs.sum()
print("  promoter background: " +
      ", ".join(f"{b}={v:.2f}" for b, v in zip("ACGT", bg_freqs)))

def scan_motif(lo_pad, L, thr):
    """Return set of gene indices with a hit on either strand."""
    M = N - L + 1
    hits = np.zeros(M, dtype=bool)
    for arr in (codes, rcodes):
        score = np.zeros(M, dtype=np.float32)
        for j in range(L):
            score += lo_pad[arr[j:j + M], j]
        hits |= score >= thr
    idx = np.flatnonzero(hits)
    gene_idx = np.unique(np.searchsorted(offsets, idx, side="right") - 1)
    return set(int(i) for i in gene_idx)

if SCAN_DONE:
    with open(CKPT, "rb") as f:
        saved = pickle.load(f)
    motif_targets, motif_names = saved["targets"], saved["names"]
    print(f"[4] Loaded checkpoint: {len(motif_targets)} motifs with hits")
else:
    motif_targets, motif_names = {}, {}
pfm_files = sorted(os.listdir(f"{MOTIF}/pfms"))
pfms = {}
for fn in pfm_files:
    name, m = load_pfm(f"{MOTIF}/pfms/{fn}")
    if m is None or m.shape[1] < 5 or m.shape[1] > 30:
        continue
    p = (m + 0.5) / (m.sum(axis=0, keepdims=True) + 2)
    pfms[fn.replace(".jaspar", "")] = (name, np.log2(p / bg_freqs[:, None]))
print(f"  usable PFMs: {len(pfms)}")

# Blocked BLAS scan: pad motifs in a block to common length (zero log-odds
# padding contributes nothing), score = sum_j OH[:, j:j+M].T @ LO[j].
BLOCK_BP = 4_000_000
K_BLOCK = 32
items = [] if SCAN_DONE else sorted(pfms.items(),
                                    key=lambda kv: kv[1][1].shape[1])
n_done = 0
while n_done < len(items):
    block = items[n_done:n_done + K_BLOCK]
    L = max(lo.shape[1] for _, (_, lo) in block)
    K = len(block)
    LO = np.zeros((L, 5, K), dtype=np.float32)
    for k, (mid, (_, lo)) in enumerate(block):
        Lm = lo.shape[1]
        LO[:Lm, :4, k] = lo.T
        LO[:Lm, 4, k] = -50.0  # N / non-ACGT penalty within motif length
    thr = np.zeros(K)
    for k, (mid, (_, lo)) in enumerate(block):
        sim = rng.choice(4, size=(2000, lo.shape[1]), p=bg_freqs)
        thr[k] = np.quantile(
            np.pad(lo, ((0, 0), (0, 0)))[sim, np.arange(lo.shape[1])]
            .sum(axis=1), 1 - P_THRESHOLD)
    hits_block = [set() for _ in range(K)]
    pos = 0
    while pos < N:
        end = min(pos + BLOCK_BP, N)
        seg = codes[pos:end]
        rseg = rcodes[pos:end]
        M = len(seg) - L + 1
        if M > 0:
            OH = np.zeros((5, len(seg)), dtype=np.float32)
            OH[seg, np.arange(len(seg))] = 1.0
            ROH = np.zeros((5, len(rseg)), dtype=np.float32)
            ROH[rseg, np.arange(len(rseg))] = 1.0
            for strand_oh in (OH, ROH):
                S = np.zeros((M, K), dtype=np.float32)
                for j in range(L):
                    S += strand_oh[:, j:j + M].T @ LO[j]
                hit_mask = S >= thr[None, :]
                for k in np.flatnonzero(hit_mask.any(axis=0)):
                    gidx = np.unique(np.searchsorted(
                        offsets, pos + np.flatnonzero(
                            hit_mask[:, k]), side="right") - 1)
                    hits_block[k].update(int(i) for i in gidx)
        pos = end
    for k, (mid, (name, _)) in enumerate(block):
        if hits_block[k]:
            motif_targets[mid] = hits_block[k]
            motif_names[mid] = name
    n_done += K
    print(f"  scanned {n_done}/{len(items)} motifs "
          f"(L<={L}, {sum(1 for h in hits_block if h)} with hits)")
    with open(CKPT, "wb") as f:
        pickle.dump({"targets": motif_targets, "names": motif_names}, f)

# hits are gene indices into gene_list -> convert to AGIs
motif_targets = {mid: {gene_list[i] for i in hits}
                 for mid, hits in motif_targets.items()}

print(f"  motifs with >=1 target gene: {len(motif_targets)}")

# ---------------------------------------------------------------- enrichment
print("[5] Motif enrichment (three comparisons) ...")
def motif_enrichment(target_set, background_set, label):
    rows = []
    for mid, hits in motif_targets.items():
        a = len(target_set & hits)
        c = len(background_set & hits)
        o, p = sps.fisher_exact([[a, len(target_set) - a],
                                 [c, len(background_set) - c]])
        rows.append({"motif": mid, "name": motif_names[mid],
                     "n_target_hits": a, "n_bg_hits": c,
                     "odds_ratio": o, "p_value": p})
    e = pd.DataFrame(rows)
    e["padj"] = cd.bh_adjust(e["p_value"].values)
    e["comparison"] = label
    return e

e_vs_uni = motif_enrichment(target, universe - target, "target_vs_universe")
e_vs_c2 = motif_enrichment(target, control2, "target_vs_polycomb_nonDE")
e_c1 = motif_enrichment(control1, control2,
                        "nonPolycomb_down_vs_polycomb_nonDE")
motif_all = pd.concat([e_vs_uni, e_vs_c2, e_c1])
motif_all.to_csv(f"{OUT}/motif_enrichment.csv", index=False)
for lbl, e in [("vs universe", e_vs_uni), ("vs Polycomb-nonDE", e_vs_c2)]:
    sig = e[e["padj"] < 0.05].sort_values("odds_ratio", ascending=False)
    print(f"  {lbl}: {len(sig)} significant; top: " +
          "; ".join(f"{r['name']} ({r['odds_ratio']:.1f})"
                    for _, r in sig.head(5).iterrows()))

# ---------------------------------------------------------------- GO
print("[6] GO over-representation ...")
KEEP_EV = {"IDA", "IMP", "IGI", "IEP", "HDA", "HMP", "HGI", "HEP",
           "TAS", "IC", "EXP", "IEA"}
go_by_gene = defaultdict(set)
gene_by_go = defaultdict(set)
with gzip.open(f"{MOTIF}/tair.gaf.gz", "rt") as gf:
    for line in gf:
        if line.startswith("!"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 15 or parts[6] not in KEEP_EV:
            continue
        agi = cd.clean_agi(parts[1])
        if agi:
            go_by_gene[agi].add(parts[4])
            gene_by_go[parts[4]].add(agi)
print(f"  GO-annotated genes: {len(go_by_gene):,}, "
      f"terms: {len(gene_by_go):,}")

term_name, term_ns, cur_term = {}, {}, None
for line in open(f"{MOTIF}/go.obo"):
    line = line.strip()
    if line == "[Term]":
        cur_term = "pending"
    elif cur_term == "pending":
        if line.startswith("id: "):
            cur_term = line[4:]
            term_name[cur_term] = ""
            term_ns[cur_term] = ""
    elif cur_term:
        if line.startswith("name: "):
            term_name[cur_term] = line[6:]
        elif line.startswith("namespace: "):
            term_ns[cur_term] = line[11:]

def go_enrichment(target_set, background_set, label, min_genes=5):
    rows = []
    bg = set(background_set)
    n_bg = max(len(bg), 1)
    for term, gset in gene_by_go.items():
        a = len(target_set & gset)
        if a < min_genes:
            continue
        c = len(bg & gset)
        o, p = sps.fisher_exact([[a, len(target_set) - a],
                                 [c, n_bg - c]], alternative="greater")
        rows.append({"GO": term, "name": term_name.get(term, ""),
                     "namespace": term_ns.get(term, ""),
                     "n_in_set": a, "n_in_bg": c,
                     "fold": (a / max(len(target_set), 1)) / (c / n_bg),
                     "p_value": p})
    e = pd.DataFrame(rows)
    if not len(e):
        print(f"  [{label}] no terms with >= {min_genes} genes in set")
        return e
    e["padj"] = cd.bh_adjust(e["p_value"].values)
    e["comparison"] = label
    return e.sort_values("p_value")

go_uni = go_enrichment(target, universe, "target_vs_universe")
go_c2 = go_enrichment(target, control2, "target_vs_polycomb_nonDE")
pd.concat([go_uni, go_c2]).to_csv(f"{OUT}/go_enrichment.csv", index=False)
for lbl, e in [("vs universe", go_uni), ("vs Polycomb-nonDE", go_c2)]:
    if not len(e):
        continue
    sig = e[e["padj"] < 0.05]
    print(f"  {lbl}: {len(sig)} significant terms; top: " +
          "; ".join(f"{r['name'][:40]} ({r['fold']:.1f}x)"
                    for _, r in sig.head(5).iterrows()))

# ---------------------------------------------------------------- replication
print("[6b] OSD-314 replication (directional consistency) ...")
de314 = pd.read_csv(f"{SWEEP}/OSD-314_DE.csv", index_col=0)
down314 = set(de314[de314["direction"] == "down"].index)
ns314 = set(de314[de314["direction"] == "ns"].index)
target314 = {g for g in down314 if g in state_of.index
             and state_of[g] in POLY and g in promoters}
universe314 = set(promoters) & set(de314.index)
print(f"  OSD-314 Polycomb-down: {len(target314)}")
go314 = go_enrichment(target314, universe314, "OSD314_target_vs_universe")
# top OSD-37 terms: do they replicate directionally?
top37 = go_uni[go_uni["padj"] < 0.05].head(10)
rep_rows = []
n314 = max(len(target314), 1)
n_u314 = max(len(universe314), 1)
for _, r in top37.iterrows():
    gset = gene_by_go[r["GO"]]
    a314 = len(target314 & gset)
    c314 = len(universe314 & gset)
    _, p314 = sps.fisher_exact(
        [[a314, len(target314) - a314], [c314, len(universe314) - c314]],
        alternative="greater")
    rep_rows.append({"GO": r["GO"], "name": r["name"],
                     "OSD37_fold": r["fold"], "OSD37_padj": r["padj"],
                     "OSD314_n": a314,
                     "OSD314_fold": (a314 / n314) / (c314 / n_u314),
                     "OSD314_p": p314})
rep = pd.DataFrame(rep_rows)
rep.to_csv(f"{OUT}/osd314_go_replication.csv", index=False)
print(rep[["name", "OSD37_fold", "OSD314_n", "OSD314_fold",
           "OSD314_p"]].to_string(index=False))

# ---------------------------------------------------------------- figures
print("[7] Figures ...")
# clean JASPAR header tabs from motif names
motif_names = {k: v.replace("\t", " ") for k, v in motif_names.items()}
for e in (e_vs_uni, e_vs_c2, e_c1):
    e["name"] = e["name"].str.replace("\t", " ", regex=False)
sns_ok = True
try:
    import seaborn as sns  # noqa
except ImportError:
    sns_ok = False

# fig22: top motifs, target vs universe (the powered comparison; the
# state-matched comparison yields no BH-significant motifs)
sig = e_vs_uni[(e_vs_uni["padj"] < 0.05) & (e_vs_uni["n_target_hits"] > 0)]
enriched = sig[sig["odds_ratio"] > 1].sort_values("odds_ratio",
                                                  ascending=False).head(12)
depleted = sig[sig["odds_ratio"] < 1].sort_values("odds_ratio").head(12)
top = pd.concat([enriched, depleted]).drop_duplicates()
if len(top):
    fig, ax = plt.subplots(figsize=(8, max(3, 0.3 * len(top))))
    colors = ["#0279EE" if o > 1 else "#9e9e9e" for o in top["odds_ratio"]]
    labels = [f"{n} ({m})" for n, m in zip(top["name"], top["motif"])]
    ax.barh(labels, top["odds_ratio"], color=colors)
    for i, (o, q) in enumerate(zip(top["odds_ratio"], top["padj"])):
        if q < 0.05:
            ax.text(max(o, 0.02), i, " *" if q >= 0.01 else " **",
                    va="center")
    ax.axvline(1, color="k", lw=0.8, ls="--")
    ax.set_xlabel("Odds ratio (Polycomb-down vs expressed genome)")
    fig.savefig(f"{FIG}/fig22_motif_enrichment.svg", bbox_inches="tight")
    fig.savefig(f"{FIG}/fig22_motif_enrichment.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

# fig23: GO dot plot, target vs Polycomb-nonDE
sig_go = go_c2[go_c2["padj"] < 0.1].head(15) if len(go_c2) else []
if len(sig_go):
    fig, ax = plt.subplots(figsize=(9, max(3, 0.35 * len(sig_go))))
    sc = ax.scatter(sig_go["fold"], range(len(sig_go)),
                    s=30 * sig_go["n_in_set"],
                    c=-np.log10(sig_go["padj"].clip(lower=1e-300)),
                    cmap="viridis")
    ax.set_yticks(range(len(sig_go)),
                  [f"{r['name'][:45]} [{r['namespace'][:4]}]"
                   for _, r in sig_go.iterrows()], fontsize=8)
    ax.axvline(1, color="k", lw=0.8, ls="--")
    ax.set_xlabel("Fold enrichment")
    plt.colorbar(sc, ax=ax, label="-log10 FDR")
    fig.savefig(f"{FIG}/fig23_go_dotplot.svg", bbox_inches="tight")
    fig.savefig(f"{FIG}/fig23_go_dotplot.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

# fig24: GC content
fig, ax = plt.subplots(figsize=(6, 4))
ax.boxplot([[gc_frac(promoters[g]) for g in s] for s in
            [target, control1, control2]],
           tick_labels=["Polycomb-down", "non-Polycomb-down",
                        "Polycomb-non-DE"])
ax.set_ylabel("Promoter GC content")
fig.savefig(f"{FIG}/fig24_gc_content.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig24_gc_content.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# fig25: motif hit counts across the three gene sets
if len(top):
    ov = pd.DataFrame({
        r["name"]: [len(motif_targets[r["motif"]] & s)
                    for s in (target, control1, control2)]
        for _, r in top.iterrows()}, index=[
            "Polycomb-down", "non-Polycomb-down", "Polycomb-non-DE"]).T
    fig, ax = plt.subplots(figsize=(7, max(3, 0.3 * len(ov))))
    left = np.zeros(len(ov))
    for col, color in zip(ov.columns, ["#FD9BED", "#0279EE", "#9e9e9e"]):
        ax.barh(ov.index, ov[col], left=left, color=color, label=col)
        left += ov[col].values
    ax.legend(frameon=False)
    ax.set_xlabel("Genes with motif hit in promoter")
    fig.savefig(f"{FIG}/fig25_motif_overlap.svg", bbox_inches="tight")
    fig.savefig(f"{FIG}/fig25_motif_overlap.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

print("Done.")
