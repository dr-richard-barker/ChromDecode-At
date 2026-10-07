"""ChromDecode-At end-to-end run on NASA OSDR OSD-37 (Arabidopsis spaceflight).

Input : GeneLab VST counts (48 samples, 4 ecotypes x FLT/GC)
Output: DE table, state assignments, enrichment tables, figure suite,
        learned-layer metrics -> $CHROMDECODE_ROOT/results + figures
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys
sys.path.insert(0, f"{ROOT}/code")
import numpy as np
import pandas as pd
import seaborn as sns
import chromdecode as cd
import learned as lrn
import figures as fg

DATA = f"{ROOT}/data"
OUT = f"{ROOT}/results"
FIG = f"{ROOT}/figures"

# ---------------------------------------------------------------- Module 1
print("[1] Loading OSD-37 VST counts ...")
vst = cd.load_counts(f"{DATA}/osdr/GLDS-37_rna_seq_VST_Counts_GLbulkRNAseq.csv")
print(f"    {vst.shape[0]} genes x {vst.shape[1]} samples")
de = cd.differential_expression(vst, list(vst.columns))
de.to_csv(f"{OUT}/OSD-37_DE_FLT_vs_GC.csv")
n_up = (de["direction"] == "up").sum()
n_down = (de["direction"] == "down").sum()
print(f"    DE: {n_up} up, {n_down} down (padj<0.05, |log2fc|>1)")

# ---------------------------------------------------------------- Module 2
print("[2] Mapping AGI -> TAIR10 coordinates ...")
genes = cd.load_gene_coordinates(f"{DATA}/araencode/genes_TAIR10.txt")
regions = cd.gene_regions(genes)

# ---------------------------------------------------------------- Module 3
print("[3] Assigning chromatin states ...")
pcsd_seg = cd.load_pcsd_segments(f"{DATA}/pcsd")
tss = cd.assign_state_by_tss(regions, pcsd_seg, value_col="pcsd_state",
                             out_col="pcsd_state")
tss["state_group"] = tss["pcsd_state"].map(lrn.state_group)

ae_seg = cd.load_araencode_states(
    f"{DATA}/araencode/seedling_12_dense.redefinecolor.bed.lst.gz",
    f"{DATA}/araencode/seedling_12_dense.redefinecolor.bed.lst.gz.tbl")
tss = cd.assign_state_by_tss(tss, ae_seg, value_col="state_name",
                             out_col="araencode_state")

try:
    acr = cd.assign_acr_overlap(tss, f"{DATA}/plantcadb/Whole_ACR_AssociatedGenes.bed.gz")
    tss["in_acr"] = acr["in_acr"]
    print(f"    genes overlapping a PlantCADB ACR: {tss['in_acr'].mean():.1%}")
except Exception as e:  # degrade gracefully; ACR flag is auxiliary
    print(f"    WARNING: PlantCADB ACR overlap unavailable ({e})")
    tss["in_acr"] = False
print(f"    PCSD coverage: {tss['pcsd_state'].notna().mean():.1%}, "
      f"AraENCODE coverage: {tss['araencode_state'].notna().mean():.1%}")
tss.to_csv(f"{OUT}/gene_state_assignments.csv", index=False)

# ---------------------------------------------------------------- Module 4
print("[4] State enrichment of DEGs ...")
bg = tss["agi"].unique()
enr_up = cd.state_enrichment(tss.rename(columns={"state_group": "st"}),
                             de[de["direction"] == "up"].index, background=bg,
                             state_col="st")
enr_down = cd.state_enrichment(tss.rename(columns={"state_group": "st"}),
                               de[de["direction"] == "down"].index,
                               background=bg, state_col="st")
enr_up.to_csv(f"{OUT}/enrichment_up.csv", index=False)
enr_down.to_csv(f"{OUT}/enrichment_down.csv", index=False)
print("    UP-regulated enrichment:")
print(enr_up.to_string(index=False))

# ---------------------------------------------------------------- Module 5
print("[5] Training learned layer (expression -> state group) ...")
mat = np.log2(vst.clip(lower=1))
common = mat.index.intersection(tss["agi"].dropna().unique())
X_raw = lrn.expression_features(mat.loc[common])
labels = tss.set_index("agi").loc[common, "state_group"]
keep = labels.notna()
X, y = X_raw[keep.values], labels[keep.values].astype(str)
# rank-transform features -> domain-invariant, comparable across datasets
Xr = lrn.rank_transform(X)
baseline = max(pd.Series(y).value_counts()) / len(y)
model, metrics, cm = lrn.train_state_model(Xr, y)
metrics["majority_baseline"] = baseline
print(f"    CV accuracy={metrics['cv_accuracy']:.3f} "
      f"(majority baseline={baseline:.3f}) "
      f"macro-F1={metrics['cv_macro_f1']:.3f} (n={metrics['n_samples']})")
print(f"    class counts: {metrics['class_counts']}")
cm.to_csv(f"{OUT}/learned_confusion_matrix.csv")
pd.Series(metrics).to_csv(f"{OUT}/learned_metrics.csv")

# External validation on AraENCODE TPM (independent tissues/samples)
print("    External validation on AraENCODE TPM matrix ...")
tpm_long = pd.read_csv(f"{DATA}/araencode/expression_TPM.txt", sep="\t",
                       header=None, names=["agi", "tissue", "tpm"])
tpm_long["agi"] = tpm_long["agi"].map(cd.clean_agi)
tpm_long = tpm_long.dropna(subset=["agi"])
g = tpm_long.groupby("agi")["tpm"]
tpm_feat = pd.DataFrame({
    "mean": np.log2(g.mean().clip(lower=1e-3) + 1),
    "std": np.log2(g.std().clip(lower=1e-3) + 1),
    "min": np.log2(g.min().clip(lower=1e-3) + 1),
    "max": np.log2(g.max().clip(lower=1e-3) + 1),
    "q25": np.log2(g.quantile(0.25).clip(lower=1e-3) + 1),
    "q75": np.log2(g.quantile(0.75).clip(lower=1e-3) + 1),
})
tpm_feat["range"] = tpm_feat["max"] - tpm_feat["min"]
tpm_feat["iqr"] = tpm_feat["q75"] - tpm_feat["q25"]
ext_common = tpm_feat.index.intersection(labels.dropna().index)
Xe = lrn.rank_transform(tpm_feat.loc[ext_common])
ye = labels.loc[ext_common].astype(str)
ext_pred = lrn.predict_states(model, Xe)
ext_acc = (ext_pred["predicted_state_group"] == ye).mean()
print(f"    External accuracy (AraENCODE TPM, {len(ext_common)} genes): "
      f"{ext_acc:.3f}")
pd.DataFrame({"external_accuracy": [ext_acc]}).to_csv(
    f"{OUT}/learned_external_validation.csv", index=False)

# Predict states for all expressed genes (the "decoded" output)
pred_all = lrn.predict_states(model, X)
pred_all.to_csv(f"{OUT}/predicted_state_groups_OSD-37.csv")

# ---------------------------------------------------------------- Module 6
print("[6] Generating figures ...")
# Fig 1: state distribution of up/down DEGs
tss_idx = tss.set_index("agi")
sets = {"up": de[de["direction"] == "up"].index, 
        "down": de[de["direction"] == "down"].index}
frac = {}
for name, ids in sets.items():
    st = tss_idx.loc[tss_idx.index.intersection(ids), "state_group"]
    frac[name] = st.value_counts(normalize=True)
frac_df = pd.DataFrame(frac).fillna(0)
all_frac = tss_idx["state_group"].value_counts(normalize=True)
frac_df["up_frac"] = frac_df.get("up", 0) / all_frac.reindex(frac_df.index).fillna(1)
frac_df["down_frac"] = frac_df.get("down", 0) / all_frac.reindex(frac_df.index).fillna(1)
fg.fig_state_distribution(frac_df[["up_frac", "down_frac"]],
                          f"{FIG}/fig1_state_distribution")

# Fig 2: enrichment bar charts
group_labels = {g: g for g in frac_df.index}
fg.fig_enrichment(enr_up, "st", group_labels, f"{FIG}/fig2_enrichment_up",
                  title="State enrichment: up-regulated DEGs (FLT vs GC)")
fg.fig_enrichment(enr_down, "st", group_labels, f"{FIG}/fig2_enrichment_down",
                  title="State enrichment: down-regulated DEGs (FLT vs GC)")

# Fig 3: volcano colored by state group
de_states = de.join(tss_idx["state_group"], how="inner")
fg.fig_volcano_by_state(de, tss_idx["state_group"], f"{FIG}/fig3_volcano_by_state")

# Fig 4: locus views for top DEGs
top_genes = pd.concat([de.sort_values("log2fc", ascending=False).head(3),
                       de.sort_values("log2fc").head(3)]).index.unique()
ae_colors = {s: c for s, c in zip(
    ae_seg["state_name"].unique(),
    fg.PALETTE * 3)}
for i, agi in enumerate(top_genes):
    row = regions[regions["agi"] == agi]
    if not len(row):
        continue
    r = row.iloc[0]
    pad = 10000
    win = (max(0, r["start"] - pad), r["end"] + pad)
    pc = pcsd_seg[(pcsd_seg["chr"] == r["chr"]) &
                  (pcsd_seg["end"] > win[0]) & (pcsd_seg["start"] < win[1])].copy()
    pc["value"] = pc["pcsd_state"].astype(str)
    ae = ae_seg[(ae_seg["chr"] == r["chr"]) &
                (ae_seg["end"] > win[0]) & (ae_seg["start"] < win[1])].copy()
    ae["value"] = ae["state_name"]
    desc = pd.read_csv(f"{DATA}/pcsd/state_descriptions.tsv", sep="\t")
    pc_colors = {str(s): c for s, c in zip(desc["state"],
                                           sns.color_palette("husl", 36))}
    fg.fig_locus_view(
        [{"label": "PCSD 36-state", "segments": pc, "value_type": "state",
          "colors": pc_colors},
         {"label": "AraENCODE 12-state", "segments": ae, "value_type": "state",
          "colors": ae_colors}],
        f"{FIG}/fig4_locus_{agi}",
        title=f"{agi} ({r['chr']}:{win[0]}-{win[1]})")

# Fig 5: heatmap of top-variable DEGs with state colors
top_var = mat.loc[de.index].var(axis=1).sort_values(ascending=False).head(50)
hm = mat.loc[top_var.index.intersection(tss_idx.index)]
row_colors = tss_idx.loc[hm.index, "state_group"].map(fg.GROUP_COLORS).fillna("#cccccc")
fg.fig_heatmap_features(hm, row_colors, f"{FIG}/fig5_heatmap_topDEGs")

# Fig 6: confusion matrix + 2D projection
fg.fig_confusion(cm, f"{FIG}/fig6_confusion_matrix")
from sklearn.decomposition import PCA
p = PCA(n_components=2, random_state=42)
emb = p.fit_transform(X)
proj = pd.DataFrame({"dim1": emb[:, 0], "dim2": emb[:, 1]}, index=X.index)
proj["state_group"] = y.values
fg.fig_som_projection(proj, f"{FIG}/fig7_state_space_projection")

print("Done. Results in results/, figures in figures/.")
