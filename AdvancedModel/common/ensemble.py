"""Independent H04 refits and a blend selected only from aligned training OOF."""

import json

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from threadpoolctl import threadpool_limits

from src.feature_engineering import FEATURE_COLUMNS
from .evaluation import ranking_metrics, validate_scores
from .training import load_mil_bundle, predict_bags


H04_PARAMETERS = {"learning_rate": .03, "max_iter": 350, "max_leaf_nodes": 31,
                  "min_samples_leaf": 20, "l2_regularization": 2.,
                  "class_weight": None, "early_stopping": False}


def make_h04(seed=42):
    return HistGradientBoostingClassifier(**H04_PARAMETERS, random_state=seed)


def run_h04_cv(store, features, run_dir, *, threads=4):
    """Fresh H04 fits: no full-training fitted model is used to generate OOF."""
    if not features.index.equals(store.sites.set_index(["dataset", "transcript_id", "transcript_position"]).index):
        raise ValueError("Tabular/read-bag identities are misaligned.")
    with (run_dir / "h04_cv_started.json").open("x") as handle:
        json.dump({"parameters": H04_PARAMETERS, "seed": store.manifest["seed"]}, handle)
    train_ids = store.ids("train")
    folds = store.sites.assessment_fold.to_numpy()
    scores = np.full(len(store.sites), np.nan)
    results = []
    for fold in range(1, store.manifest["cv_folds"] + 1):
        fit, assess = train_ids[folds[train_ids] != fold], train_ids[folds[train_ids] == fold]
        if set(store.sites.iloc[fit].group_id) & set(store.sites.iloc[assess].group_id):
            raise ValueError("H04 fitting/assessment groups overlap.")
        model = make_h04(store.manifest["seed"])
        with threadpool_limits(limits=threads):
            model.fit(features.iloc[fit][FEATURE_COLUMNS], store.sites.iloc[fit].label)
            scores[assess] = model.predict_proba(features.iloc[assess][FEATURE_COLUMNS])[:, 1]
        results.append({"fold": fold, **ranking_metrics(store.sites.iloc[assess].label, scores[assess])})
        print(f"H04 fold {fold}: AP={results[-1]['average_precision']:.4f}", flush=True)
        pd.DataFrame(results).to_csv(run_dir / "h04_cv_fold_metrics.csv", index=False)
    predictions = store.indexed_metadata(train_ids).assign(score=scores[train_ids])
    predictions.to_parquet(run_dir / "h04_oof_predictions.parquet")
    return predictions, pd.DataFrame(results)


def blend_scores(mil_scores, tree_scores, mil_weight):
    if not np.isfinite(mil_weight) or not 0 <= mil_weight <= 1:
        raise ValueError("MIL weight must lie in [0, 1].")
    _, mil_scores = validate_scores(np.zeros(len(mil_scores)), mil_scores)
    _, tree_scores = validate_scores(np.zeros(len(mil_scores)), tree_scores)
    return mil_weight * mil_scores + (1 - mil_weight) * tree_scores


def select_oof_blend(mil_oof, tree_oof, weights=None):
    """Tune convex weight by mean fold AP, ROC-AUC tie-break, then lower MIL weight.

    Metrics at the selected weight are development/selection scores, not an
    independent estimate of the complete weight-selection procedure.
    """
    if not mil_oof.index.is_unique or not mil_oof.index.equals(tree_oof.index):
        raise ValueError("OOF site identities/order must match exactly.")
    for column in ["label", "group_id", "assessment_fold", "split"]:
        if not mil_oof[column].equals(tree_oof[column]):
            raise ValueError(f"OOF metadata mismatch: {column}")
    if not mil_oof["split"].eq("train").all():
        raise ValueError("Blend selection only accepts outer-training OOF predictions.")
    if not mil_oof.groupby("group_id").assessment_fold.nunique().eq(1).all():
        raise ValueError("OOF genes/constructs cross assessment folds.")
    fold_ids = sorted(mil_oof.assessment_fold.unique())
    if len(fold_ids) < 2 or min(fold_ids) < 1:
        raise ValueError("Need at least two valid OOF folds.")
    weights = np.linspace(0., 1., 21) if weights is None else weights
    rows = []
    for weight in weights:
        scores = blend_scores(mil_oof.score, tree_oof.score, weight)
        fold_metrics = []
        for fold in fold_ids:
            mask = mil_oof.assessment_fold.eq(fold).to_numpy()
            if mil_oof.loc[mask, "label"].nunique() != 2:
                raise ValueError("Every OOF assessment fold must contain both classes.")
            fold_metrics.append(ranking_metrics(mil_oof.loc[mask, "label"], scores[mask]))
        table = pd.DataFrame(fold_metrics)
        rows.append({"mil_weight": float(weight), "h04_weight": float(1 - weight),
                     **{f"{column}_mean": table[column].mean() for column in table},
                     **{f"{column}_sd": table[column].std(ddof=1) for column in table}})
    candidates = pd.DataFrame(rows).sort_values(
        ["average_precision_mean", "roc_auc_mean", "mil_weight"], ascending=[False, False, True],
    ).reset_index(drop=True)
    return float(candidates.iloc[0].mil_weight), candidates


def fit_final_h04(store, features, run_dir, *, threads=4):
    ids = store.ids("train")
    model = make_h04(store.manifest["seed"])
    with threadpool_limits(limits=threads):
        model.fit(features.iloc[ids][FEATURE_COLUMNS], store.sites.iloc[ids].label)
    joblib.dump({"model": model, "model_name": "H04 full-training component", "feature_columns": FEATURE_COLUMNS,
                 "data_fingerprint": store.manifest["fingerprint"], "kmer_residuals": False},
                run_dir / "h04.joblib", compress=3)
    return model


def predict_h04(model, features, ids, threads=4):
    with threadpool_limits(limits=threads):
        return model.predict_proba(features.iloc[ids][FEATURE_COLUMNS])[:, 1]


def save_ensemble_bundle(run_dir, mil_weight, threshold, store):
    if not (run_dir / "model.pt").is_file() or not (run_dir / "h04.joblib").is_file():
        raise ValueError("Save both fitted component models first.")
    payload = {"format_version": 1, "mil_path": "model.pt", "h04_path": "h04.joblib",
               "mil_weight": float(mil_weight), "threshold": float(threshold),
               "data_fingerprint": store.manifest["fingerprint"],
               "weight_selection": "highest mean training OOF fold AP; no validation/test labels"}
    (run_dir / "ensemble.json").write_text(json.dumps(payload, indent=2))


def predict_saved_ensemble(run_dir, store, ids, features, device="cpu"):
    specification = json.loads((run_dir / "ensemble.json").read_text())
    if specification.get("format_version") != 1:
        raise ValueError("Unsupported ensemble version.")
    expected_index = store.indexed_metadata(ids).index
    if not features.iloc[ids].index.equals(expected_index):
        raise ValueError("Ensemble feature/site identities differ.")
    model, normalizer, config, _ = load_mil_bundle(run_dir / specification["mil_path"], device)
    tree = joblib.load(run_dir / specification["h04_path"])["model"]
    return blend_scores(predict_bags(model, normalizer, store, ids, config),
                        predict_h04(tree, features, ids, config.threads), specification["mil_weight"])
