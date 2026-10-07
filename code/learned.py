"""ChromDecode-At module 5: learned layer.

Elastic-net multinomial model predicting a gene's chromatin-state group from
expression-distribution features. State groups collapse the 36 PCSD states
into 5 biologically interpretable classes.
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

STATE_GROUPS = {
    "Active_transcribed": list(range(1, 11)) + [22, 25, 26, 27, 28],
    "Polycomb_repressed": [11, 12, 13, 14, 15],
    "Accessible_promoter": [16, 17, 18, 19, 20, 21, 23, 24],
    "Intergenic_quiet": [29, 30],
    "Heterochromatin_TE": [31, 32, 33, 34, 35, 36],
}
STATE_TO_GROUP = {s: g for g, ss in STATE_GROUPS.items() for s in ss}


def state_group(state_id):
    if state_id is None or (isinstance(state_id, float) and np.isnan(state_id)):
        return None
    return STATE_TO_GROUP.get(int(state_id), None)


def expression_features(mat):
    """Distribution-summary features from a genes x samples log-expression matrix.

    mat: DataFrame (genes x samples), log-scale values (VST or log2 TPM).
    Returns DataFrame indexed like mat with columns
    [mean, std, min, max, q25, q75, range, iqr].
    """
    v = mat.values
    f = pd.DataFrame({
        "mean": v.mean(axis=1),
        "std": v.std(axis=1),
        "min": v.min(axis=1),
        "max": v.max(axis=1),
        "q25": np.quantile(v, 0.25, axis=1),
        "q75": np.quantile(v, 0.75, axis=1),
    }, index=mat.index)
    f["range"] = f["max"] - f["min"]
    f["iqr"] = f["q75"] - f["q25"]
    return f


def rank_transform(X):
    """Domain-invariant feature transform: per-column quantile ranks in [0,1].

    Makes feature marginal distributions identical across datasets
    (e.g. VST counts vs TPM), enabling cross-dataset transfer.
    """
    from scipy.stats import rankdata
    R = np.zeros_like(X.values, dtype=float)
    for j in range(X.shape[1]):
        R[:, j] = rankdata(X.values[:, j]) / len(X)
    return pd.DataFrame(R, index=X.index, columns=X.columns)


def genomic_block_folds(regions_index, n_blocks=5, block_size_bp=100_000,
                        genes=None, seed=42):
    """Assign genes to folds by genomic blocks (whole 100-kb bins per fold).

    genes: DataFrame with agi, chr, start (e.g. regions table). Genes missing
    from it get fold -1 (excluded from block-CV evaluation).
    Returns pd.Series: index=agi, value=fold id (0..n_blocks-1) or -1.
    """
    rng = np.random.RandomState(seed)
    g = genes.set_index("agi")
    blocks = (g["start"] // block_size_bp).astype(int)
    key = list(zip(g["chr"], blocks))
    uniq = sorted(set(key))
    rng.shuffle(uniq)
    block_fold = {b: i % n_blocks for i, b in enumerate(uniq)}
    folds = pd.Series([block_fold[k] for k in key], index=g.index)
    return folds.reindex(regions_index).fillna(-1).astype(int)


def nested_cv(X, y, n_outer=5, C_grid=(0.1, 1.0, 10.0), n_inner=3,
              random_state=42):
    """Nested CV: inner folds tune C, outer folds report honest metrics.

    Returns (outer_metrics dict, confusion matrix df, best C counter).
    """
    from sklearn.model_selection import StratifiedKFold
    outer = StratifiedKFold(n_splits=n_outer, shuffle=True,
                            random_state=random_state)
    y = np.asarray(y)
    y_pred = np.empty_like(y)
    chosen = []
    for tr, te in outer.split(X, y):
        inner = StratifiedKFold(n_splits=n_inner, shuffle=True,
                                random_state=random_state)
        best_c, best_score = None, -1
        for C in C_grid:
            scores = []
            for itr, ival in inner.split(X.iloc[tr], y[tr]):
                m = _make_model(C).fit(X.iloc[tr].iloc[itr], y[tr][itr])
                scores.append(accuracy_score(y[tr][ival],
                                             m.predict(X.iloc[tr].iloc[ival])))
            s = float(np.mean(scores))
            if s > best_score:
                best_score, best_c = s, C
        chosen.append(best_c)
        m = _make_model(best_c).fit(X.iloc[tr], y[tr])
        y_pred[te] = m.predict(X.iloc[te])
    classes = sorted(set(y))
    metrics = {
        "outer_cv_accuracy": accuracy_score(y, y_pred),
        "outer_cv_macro_f1": f1_score(y, y_pred, average="macro"),
        "majority_baseline": float(pd.Series(y).value_counts().max() / len(y)),
        "n_samples": len(y),
        "chosen_C": pd.Series(chosen).value_counts().to_dict(),
    }
    cm = confusion_matrix(y, y_pred, labels=classes)
    cm_df = pd.DataFrame(cm, index=[f"true_{c}" for c in classes],
                         columns=[f"pred_{c}" for c in classes])
    return metrics, cm_df, pd.Series(chosen).value_counts()


def _make_model(C):
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=0.5,
            max_iter=2000, class_weight="balanced", C=C,
            random_state=42)),
    ])


def train_state_model(X, y, n_splits=5, random_state=42):
    """Fit elastic-net multinomial logistic regression with stratified CV.

    Returns (fitted model on all data, cv_metrics dict, confusion matrix df).
    Feature scaling and the model are wrapped in a pipeline; no feature
    selection leakage (all features are per-gene summaries, no filtering).
    """
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=0.5,
            max_iter=2000, class_weight="balanced", C=1.0,
            random_state=random_state)),
    ])
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True,
                          random_state=random_state)
    y_pred = cross_val_predict(pipe, X, y, cv=skf, n_jobs=1)
    classes = sorted(pd.unique(y))
    metrics = {
        "cv_accuracy": accuracy_score(y, y_pred),
        "cv_macro_f1": f1_score(y, y_pred, average="macro"),
        "n_samples": len(y),
        "n_features": X.shape[1],
        "class_counts": pd.Series(y).value_counts().to_dict(),
    }
    cm = confusion_matrix(y, y_pred, labels=classes)
    cm_df = pd.DataFrame(cm, index=[f"true_{c}" for c in classes],
                         columns=[f"pred_{c}" for c in classes])
    pipe.fit(X, y)
    return pipe, metrics, cm_df


def predict_states(pipe, X):
    """Predict state group + class probabilities for new genes."""
    proba = pd.DataFrame(pipe.predict_proba(X), index=X.index,
                         columns=pipe.classes_)
    pred = proba.idxmax(axis=1).rename("predicted_state_group")
    conf = proba.max(axis=1).rename("confidence")
    return pd.concat([pred, conf, proba], axis=1)
