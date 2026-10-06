"""Aligned pooled/per-dataset metrics, validation thresholds and saved reports."""

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score, auc, average_precision_score, balanced_accuracy_score,
    confusion_matrix, ConfusionMatrixDisplay, f1_score, precision_recall_curve,
    precision_score, recall_score, roc_auc_score, roc_curve,
)


def validate_scores(labels, scores):
    labels, scores = np.asarray(labels), np.asarray(scores, dtype=float)
    if labels.ndim != 1 or scores.shape != labels.shape or not len(labels):
        raise ValueError("Scores must have the same nonempty one-dimensional shape as labels.")
    if not np.isin(labels, [0, 1]).all() or not np.isfinite(scores).all() or ((scores < 0) | (scores > 1)).any():
        raise ValueError("Expected binary labels and finite scores between zero and one.")
    return labels, scores


def ranking_metrics(labels, scores):
    labels, scores = validate_scores(labels, scores)
    if len(np.unique(labels)) < 2:
        return {"roc_auc": np.nan, "average_precision": np.nan, "pr_auc_trapezoid": np.nan}
    precision, recall, _ = precision_recall_curve(labels, scores)
    return {"roc_auc": roc_auc_score(labels, scores), "average_precision": average_precision_score(labels, scores),
            "pr_auc_trapezoid": auc(recall, precision)}


def select_threshold(labels, scores):
    labels, scores = validate_scores(labels, scores)
    if len(np.unique(labels)) != 2:
        raise ValueError("Threshold selection requires both classes in validation.")
    precision, recall, thresholds = precision_recall_curve(labels, scores)
    f1 = np.divide(2 * precision[:-1] * recall[:-1], precision[:-1] + recall[:-1],
                   out=np.zeros_like(thresholds), where=(precision[:-1] + recall[:-1]) > 0)
    best = np.flatnonzero(f1 == f1.max())[-1]
    return float(thresholds[best]), pd.DataFrame({"threshold": thresholds, "precision": precision[:-1],
                                               "recall": recall[:-1], "f1": f1})


def threshold_metrics(labels, scores, threshold):
    labels, scores = validate_scores(labels, scores)
    if not np.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError("Threshold must lie in [0, 1].")
    prediction = scores >= threshold
    tn, fp, fn, tp = confusion_matrix(labels, prediction, labels=[0, 1]).ravel()
    return {"threshold": threshold, "precision": precision_score(labels, prediction, zero_division=0),
            "recall": recall_score(labels, prediction, zero_division=0),
            "f1": f1_score(labels, prediction, zero_division=0),
            "accuracy": accuracy_score(labels, prediction),
            "balanced_accuracy": balanced_accuracy_score(labels, prediction),
            "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
            "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)}


def metric_report(metadata, scores, threshold=None):
    """metadata must be in prediction order, with dataset as a column."""
    labels, scores = validate_scores(metadata.label, scores)
    rows = []
    for name in ["pooled", *sorted(metadata.dataset.unique())]:
        mask = np.ones(len(metadata), dtype=bool) if name == "pooled" else metadata.dataset.eq(name).to_numpy()
        y = labels[mask]
        record = {"dataset": name, "sites": int(mask.sum()), "groups": metadata.loc[mask, "group_id"].nunique(),
                  "positives": int(y.sum()), "positive_rate": y.mean(),
                  "metric_status": "ok" if len(np.unique(y)) == 2 else "single_class",
                  "prior_ap_reference": y.mean(), "always_positive_f1": 2 * y.mean() / (1 + y.mean()),
                  **ranking_metrics(y, scores[mask])}
        if "assessment_fold" in metadata:
            record["assessment_folds"] = metadata.loc[mask, "assessment_fold"].loc[lambda s: s > 0].nunique()
        if threshold is not None:
            record.update(threshold_metrics(y, scores[mask], threshold))
        rows.append(record)
    return pd.DataFrame(rows)


def bootstrap_intervals(metadata, scores, threshold, repeats=500, seed=42):
    """Resample groups within human/synthetic strata; a lone construct stays fixed."""
    labels, scores = validate_scores(metadata.label, scores)
    codes, names = pd.factorize(metadata.group_id, sort=True)
    groups = metadata[["group_id", "group_kind"]].drop_duplicates()
    if groups.group_id.duplicated().any():
        raise ValueError("Each group must have exactly one group kind.")
    kinds = groups.set_index("group_id").group_kind.reindex(names).to_numpy()
    pools = [np.flatnonzero(kinds == name) for name in sorted(set(kinds))]
    rng, samples = np.random.default_rng(seed), []
    for _ in range(repeats):
        chosen = np.concatenate([pool[rng.integers(len(pool), size=len(pool))] for pool in pools])
        weights = np.bincount(chosen, minlength=len(names))[codes]
        if len(np.unique(labels[weights > 0])) < 2:
            continue
        samples.append({"average_precision": average_precision_score(labels, scores, sample_weight=weights),
                        "roc_auc": roc_auc_score(labels, scores, sample_weight=weights),
                        "f1": f1_score(labels, scores >= threshold, sample_weight=weights, zero_division=0)})
    if repeats < 1 or len(samples) < repeats * .9:
        raise ValueError("Insufficient valid grouped bootstrap draws.")
    draws = pd.DataFrame(samples)
    point = {**ranking_metrics(labels, scores), **threshold_metrics(labels, scores, threshold)}
    return pd.DataFrame({"estimate": pd.Series({name: point[name] for name in draws}),
                         "lower_95": draws.quantile(.025), "upper_95": draws.quantile(.975)})


def save_evaluation(store, val_scores, test_scores, threshold, run_dir, *, bootstrap_repeats=500):
    results = {}
    for split, scores in [("val", val_scores), ("test", test_scores)]:
        ids = store.ids(split)
        report = metric_report(store.sites.iloc[ids], scores, threshold)
        report.to_csv(run_dir / f"{split}_metrics_by_dataset.csv", index=False)
        predictions = store.indexed_metadata(ids).assign(score=scores, prediction=(scores >= threshold).astype("int8"))
        predictions.to_parquet(run_dir / f"{split}_predictions.parquet")
        results[split] = report
    intervals = bootstrap_intervals(store.sites.iloc[store.ids("test")], test_scores, threshold,
                                    repeats=bootstrap_repeats, seed=store.manifest["seed"])
    intervals.to_csv(run_dir / "test_group_bootstrap_intervals.csv")
    results["intervals"] = intervals
    summary = {"threshold": threshold, "validation": results["val"].to_dict(orient="records"),
               "test": results["test"].to_dict(orient="records"),
               "limits": ["Internal evaluation of previously explored datasets; test is not used for tuning.",
                          "Only one synthetic construct in test; pooled intervals cannot estimate between-construct uncertainty.",
                          "data2 binary labels mean positive mixture, not every read modified.",
                          "Attention weights and weakly supervised read scores are not validated read-level modification truth."]}
    (run_dir / "evaluation.json").write_text(json.dumps(summary, indent=2))
    return results


def plot_evaluation(metadata, scores, threshold, output_path=None):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    for name, table in metadata.groupby("dataset", sort=True):
        mask = metadata.dataset.eq(name).to_numpy()
        y, values = table.label.to_numpy(), np.asarray(scores)[mask]
        if len(np.unique(y)) < 2:
            continue
        p, r, _ = precision_recall_curve(y, values)
        fpr, tpr, _ = roc_curve(y, values)
        axes[0].step(r, p, where="post", label=f"{name}: AP={average_precision_score(y, values):.3f}")
        axes[1].plot(fpr, tpr, label=f"{name}: AUC={roc_auc_score(y, values):.3f}")
    axes[0].set(title="Test precision–recall by dataset", xlabel="Recall", ylabel="Precision", xlim=(0, 1), ylim=(0, 1.02))
    axes[1].set(title="Test ROC by dataset", xlabel="False positive rate", ylabel="True positive rate", xlim=(0, 1), ylim=(0, 1.02))
    axes[1].plot([0, 1], [0, 1], "--", color="gray")
    for ax in axes[:2]:
        ax.legend()
    ConfusionMatrixDisplay.from_predictions(metadata.label, np.asarray(scores) >= threshold,
                                            labels=[0, 1], colorbar=False, ax=axes[2])
    axes[2].set_title(f"Pooled test at threshold {threshold:.3f}")
    fig.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    return fig
