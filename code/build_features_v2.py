"""Build tissue-matched training features for learned layer v2.

Parses AraENCODE expression_TPM.txt (long format: agi, tissue, tpm) into:
  - seedling-block distribution features (tissue-matched to state definitions)
  - tissue-level medians for all 16 tissues (for specificity features)
Saves to $CHROMDECODE_ROOT/results/features_v2_seedling.csv (+ tissue medians).
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys
sys.path.insert(0, f"{ROOT}/code")
import numpy as np
import pandas as pd
import chromdecode as cd

DATA = f"{ROOT}/data"
OUT = f"{ROOT}/results"

print("Streaming TPM file ...")
seed_rows = []          # (agi, log2tpm) for Seedling
med_sums = {}           # tissue -> dict(agi -> list) too big; use two-pass:
tissue_vals = {}        # tissue -> {agi: [values]} only via aggregate

# Pass 1: collect seedling values per gene and per-tissue running sums for medians
# For exact medians per tissue we need all values; instead compute tissue MEDIANS
# from per-tissue value lists capped at memory-safe size using arrays per tissue.
from collections import defaultdict
tv = defaultdict(lambda: defaultdict(list))
n = 0
with open(f"{DATA}/araencode/expression_TPM.txt") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        agi, tissue, tpm = parts[0], parts[1], parts[2]
        agi = cd.clean_agi(agi)
        if agi is None:
            continue
        try:
            v = float(tpm)
        except ValueError:
            continue
        tv[tissue][agi].append(v)
        n += 1
        if n % 5_000_000 == 0:
            print(f"  {n:,} rows ...")

print(f"parsed {n:,} rows across {len(tv)} tissues")
for t in tv:
    print(f"  {t}: {len(tv[t]):,} genes, {sum(len(v) for v in tv[t].values()):,} values")

# Distribution features for selected tissues (seedling = training; leaf/root = transfer)
print("Computing distribution features for Seedling / Leaf / Root ...")
def dist_features(d):
    feat = {}
    for agi, vals in d.items():
        a = np.log2(np.clip(np.array(vals), 1e-3, None) + 1)
        feat[agi] = [a.mean(), a.std(), a.min(), a.max(),
                     np.quantile(a, 0.25), np.quantile(a, 0.75)]
    fs = pd.DataFrame.from_dict(feat, orient="index",
                                columns=["mean", "std", "min", "max",
                                         "q25", "q75"])
    fs["range"] = fs["max"] - fs["min"]
    fs["iqr"] = fs["q75"] - fs["q25"]
    fs.index.name = "agi"
    return fs

for tissue in ["Seedling", "Leaf", "Root"]:
    fs = dist_features(tv[tissue])
    fs.to_csv(f"{OUT}/features_v2_{tissue.lower()}.csv")
    print(f"  {tissue}: {len(fs):,} genes -> features_v2_{tissue.lower()}.csv")

fs = pd.read_csv(f"{OUT}/features_v2_seedling.csv", index_col="agi")

# Tissue medians -> specificity features
print("Computing tissue-specificity features ...")
tiss_med = pd.DataFrame({t: {agi: np.median(np.log2(np.clip(np.array(v), 1e-3, None) + 1))
                             for agi, v in d.items()} for t, d in tv.items()})
tiss_med.to_csv(f"{OUT}/features_v2_tissue_medians.csv")

eps = 1e-3
tmax = tiss_med.max(axis=1)
tmin = tiss_med.replace(0, np.nan).min(axis=1)
# tau (expression breadth, Yanai et al.): sum(1 - x_i/x_max) / (n_tissues - 1)
x = tiss_med.clip(lower=-10)
tau = (1 - x.div(tmax, axis=0)).sum(axis=1) / (tiss_med.shape[1] - 1)
fs["tau"] = tau.reindex(fs.index)
fs["tissue_max_min_fold"] = (tmax - tmin).reindex(fs.index)
fs["seedling_vs_mean_fold"] = (tiss_med["Seedling"] - tiss_med.mean(axis=1)).reindex(fs.index)
fs = fs.rename(columns={c: "seed_" + c for c in
                        ["mean", "std", "min", "max", "q25", "q75",
                         "range", "iqr"]})
fs.index.name = "agi"
fs.to_csv(f"{OUT}/features_v2_seedling.csv")
print(f"saved features for {len(fs):,} genes -> features_v2_seedling.csv")
print(fs.describe().round(2).to_string())
