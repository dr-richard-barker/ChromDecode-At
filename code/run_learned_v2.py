"""Learned layer v2: tissue-matched training + honest evaluation.

Train on AraENCODE seedling expression features (tissue-matched to the
chromatin-state definitions), evaluate with nested CV and genomic-block CV,
then run transfer tests to OSD-37 (spaceflight seedlings), leaf, and root.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys
sys.path.insert(0, f"{ROOT}/code")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, f1_score
import chromdecode as cd
import learned as lrn
import figures as fg

DATA = f"{ROOT}/data"
OUT = f"{ROOT}/results"
FIG = f"{ROOT}/figures"

# ---------------------------------------------------------------- load
print("[1] Loading features and labels ...")
X_seed = pd.read_csv(f"{OUT}/features_v2_seedling.csv", index_col="agi")
tss = pd.read_csv(f"{OUT}/gene_state_assignments.csv")
labels = tss.set_index("agi")["state_group"]

common = X_seed.index.intersection(labels.dropna().index)
X = lrn.rank_transform(X_seed.loc[common])   # rank-transform ONCE, used everywhere
y = labels.loc[common].astype(str)
baseline = float(pd.Series(y).value_counts().max() / len(y))
print(f"    {len(X)} genes, {X.shape[1]} features, baseline={baseline:.3f}")
print(f"    class counts: {pd.Series(y).value_counts().to_dict()}")

# ---------------------------------------------------------------- nested CV
import os
metrics_path = f"{OUT}/learned_v2_metrics.csv"
cm_path = f"{OUT}/learned_v2_confusion.csv"
if os.path.exists(metrics_path) and os.path.exists(cm_path):
    print("[2] Loading cached CV metrics ...")
    all_m = pd.read_csv(metrics_path, index_col=0)["v2"].to_dict()
    metrics_rand = {k: all_m[k] for k in
                    ["outer_cv_accuracy", "outer_cv_macro_f1",
                     "majority_baseline"] if k in all_m}
    cm_rand = pd.read_csv(f"{OUT}/learned_v2_confusion.csv", index_col=0)
else:
    print("[2] Nested CV (random gene folds, inner C tuning) ...")
    metrics_rand, cm_rand, chosen = lrn.nested_cv(X, y)
    print(f"    outer accuracy={metrics_rand['outer_cv_accuracy']:.3f} "
          f"(baseline={baseline:.3f}) "
          f"macro-F1={metrics_rand['outer_cv_macro_f1']:.3f}")
    print(f"    chosen C: {metrics_rand['chosen_C']}")
    cm_rand.to_csv(f"{OUT}/learned_v2_confusion.csv")

# ---------------------------------------------------------------- block CV
if os.path.exists(cm_path):
    print("[3] Block CV metrics loaded from cache ...")
    metrics_block = {}
else:
    print("[3] Genomic-block CV (100-kb bins, whole blocks per fold) ...")
    genes = cd.load_gene_coordinates(f"{DATA}/araencode/genes_TAIR10.txt")
    folds = lrn.genomic_block_folds(X.index, genes=genes)
    mask = folds.values >= 0
    Xb, yb, fb = X.values[mask], y.values[mask], folds.values[mask]
    y_pred_b = np.empty_like(yb)
    for k in sorted(set(fb)):
        tr, te = np.where(fb != k)[0], np.where(fb == k)[0]
        if len(te) == 0 or len(set(yb[te])) < 2:
            continue
        m = lrn._make_model(1.0).fit(Xb[tr], yb[tr])
        y_pred_b[te] = m.predict(Xb[te])
    metrics_block = {
        "block_cv_accuracy": accuracy_score(yb, y_pred_b),
        "block_cv_macro_f1": f1_score(yb, y_pred_b, average="macro"),
        "block_majority_baseline": float(pd.Series(yb).value_counts().max() / len(yb)),
        "n_genes_block_cv": int(mask.sum()),
    }
    print(f"    block accuracy={metrics_block['block_cv_accuracy']:.3f} "
          f"(baseline={metrics_block['block_majority_baseline']:.3f})")

if not os.path.exists(cm_path):
    pd.DataFrame({**metrics_rand, **metrics_block}, index=["v2"]).T.to_csv(
        metrics_path)

# ---------------------------------------------------------------- final fit
print("[4] Final fit on all training genes ...")
model = lrn._make_model(1.0).fit(X, y)

# ---------------------------------------------------------------- transfer
print("[5] Transfer tests (no refitting, rank-transformed features) ...")
transfer = {}

def transfer_eval(feature_csv, name):
    Xt_raw = pd.read_csv(feature_csv, index_col=0)
    Xt_raw = Xt_raw.rename(columns={c: "seed_" + c for c in Xt_raw.columns
                                    if not c.startswith("seed_")
                                    and c not in
                                    ("tau", "tissue_max_min_fold",
                                     "seedling_vs_mean_fold")})
    common_t = Xt_raw.index.intersection(X.index)
    # align columns to training feature space (tissue-specificity features
    # are gene-intrinsic and shared; distribution features come from the
    # external tissue)
    Xte = Xt_raw.loc[common_t, X.columns.intersection(Xt_raw.columns)].copy()
    for c in X.columns.difference(Xte.columns):
        Xte[c] = X.loc[common_t, c]  # gene-intrinsic features from training table
    Xte = Xte[X.columns]
    Xte_r = lrn.rank_transform(Xte)
    yt = y.loc[common_t]
    acc = float((model.predict(Xte_r) == yt.values).mean())
    return {"accuracy": acc, "baseline": baseline, "n": len(common_t)}

transfer["AraENCODE_seedling(held-out study)"] = transfer_eval(
    f"{OUT}/features_v2_seedling.csv", "seedling")
transfer["AraENCODE_leaf"] = transfer_eval(
    f"{OUT}/features_v2_leaf.csv", "leaf")
transfer["AraENCODE_root"] = transfer_eval(
    f"{OUT}/features_v2_root.csv", "root")

# OSD-37 spaceflight seedlings (independent study, same tissue)
vst = cd.load_counts(f"{DATA}/osdr/GLDS-37_rna_seq_VST_Counts_GLbulkRNAseq.csv")
mat37 = np.log2(vst.clip(lower=1))
f37 = lrn.expression_features(mat37)
f37.to_csv(f"{OUT}/features_v2_osd37.csv")
transfer["OSD-37_spaceflight_seedling"] = transfer_eval(
    f"{OUT}/features_v2_osd37.csv", "osd37")

tr_df = pd.DataFrame(transfer).T
tr_df.to_csv(f"{OUT}/learned_v2_transfer.csv")
print(tr_df.round(3).to_string())

# ---------------------------------------------------------------- figures
print("[6] Figures ...")
fg.fig_confusion(cm_rand, f"{FIG}/fig8_v2_confusion_matrix")

fig, ax = plt.subplots(figsize=(7, 4))
names = list(tr_df.index)
accs = tr_df["accuracy"].astype(float).values
ax.bar(range(len(names)), accs, color="#0279EE", label="transfer accuracy")
ax.axhline(baseline, color="k", ls="--", lw=1, label="majority baseline")
ax.set_xticks(range(len(names)), names, rotation=15, ha="right")
ax.set_ylabel("Accuracy (state-group prediction)")
ax.legend(frameon=False)
fig.savefig(f"{FIG}/fig9_v2_transfer.svg", bbox_inches="tight")
fig.savefig(f"{FIG}/fig9_v2_transfer.png", dpi=150, bbox_inches="tight")
print("Done.")
