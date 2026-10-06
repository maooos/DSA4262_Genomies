"""Exercise notebook metric helpers without running the full data pipeline."""

import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import (
    accuracy_score,
    auc,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
)


@pytest.fixture(scope="module")
def helpers():
    notebook = json.loads(
        (Path(__file__).resolve().parents[1]
         / "notebooks/data0_modelling_pipeline.ipynb").read_text()
    )
    namespace = {
        **globals(), "BOOTSTRAP_REPEATS": 20, "SEED": 42,
    }
    helper_names = {
        "ranking_metrics", "threshold_metrics", "select_f1_threshold",
        "gene_bootstrap_metrics",
        "dataset_metrics",
    }
    for cell in notebook["cells"]:
        if cell["cell_type"] != "code":
            continue
        functions = [
            node for node in ast.parse("".join(cell["source"])).body
            if isinstance(node, ast.FunctionDef) and node.name in helper_names
        ]
        if functions:
            exec(compile(ast.Module(body=functions, type_ignores=[]),
                         "<notebook-metrics>", "exec"), namespace)
    assert helper_names <= namespace.keys()
    return namespace


def test_constant_baseline_ap_equals_prevalence(helpers):
    labels = pd.Series([0, 0, 0, 1])
    result = helpers["ranking_metrics"](labels, np.full(4, 0.1))
    assert result["average_precision"] == pytest.approx(0.25)
    assert result["roc_auc"] == pytest.approx(0.5)
    # Linear interpolation is different from AP for a constant predictor.
    assert result["pr_auc_trapezoid"] == pytest.approx(0.625)


def test_ranking_rejects_nonfinite_or_misaligned_scores(helpers):
    for scores in [[0.1], [0.1, np.nan]]:
        with pytest.raises(ValueError, match="aligned"):
            helpers["ranking_metrics"]([0, 1], scores)


def test_threshold_matches_exhaustive_search_and_largest_tie(helpers):
    labels = np.array([1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1])
    scores = np.linspace(0.99, 0.09, len(labels))
    candidates = np.unique(scores)
    brute_f1 = np.array([f1_score(labels, scores >= t) for t in candidates])
    expected = candidates[brute_f1 == brute_f1.max()].max()
    threshold, table = helpers["select_f1_threshold"](labels, scores)
    assert threshold == expected == scores[0]
    assert len(table) == len(candidates)
    assert f1_score(labels, scores >= threshold) == pytest.approx(table.f1.max())


def test_threshold_equality_predicts_positive_and_counts_are_aligned(helpers):
    metrics = helpers["threshold_metrics"]([0, 1, 1, 0], [0.2, 0.5, 0.3, 0.6], 0.5)
    assert [metrics[name] for name in ["tn", "fp", "fn", "tp"]] == [1, 1, 1, 1]
    assert metrics["precision"] == metrics["recall"] == metrics["f1"] == 0.5


def test_gene_bootstrap_matches_explicit_whole_gene_resampling(helpers):
    labels = pd.Series([0, 1, 0, 1, 0, 1], index=list("abcdef"))
    genes = pd.Series(["A", "A", "B", "B", "C", "C"], index=labels.index)
    scores = np.array([0.05, 0.9, 0.8, 0.7, 0.4, 0.3])
    actual = helpers["gene_bootstrap_metrics"](
        labels, scores, genes, threshold=0.5, n_bootstrap=20, seed=7,
    )
    rng = np.random.default_rng(7)
    gene_names = sorted(genes.unique())
    expected = []
    for _ in range(20):
        selected = rng.integers(0, len(gene_names), size=len(gene_names))
        rows = np.concatenate([
            np.flatnonzero(genes.to_numpy() == gene_names[i]) for i in selected
        ])
        expected.append({
            "average_precision": average_precision_score(labels.iloc[rows], scores[rows]),
            "roc_auc": roc_auc_score(labels.iloc[rows], scores[rows]),
            "f1": f1_score(labels.iloc[rows], scores[rows] >= 0.5, zero_division=0),
        })
    pd.testing.assert_frame_equal(actual, pd.DataFrame(expected))
    with pytest.raises(ValueError, match="align"):
        helpers["gene_bootstrap_metrics"](labels, scores, genes.iloc[::-1], 0.5)


def test_pooled_reports_preserve_dataset_alignment_and_shared_threshold(helpers):
    index = pd.MultiIndex.from_tuples(
        [(name, "same_tx", pos) for name in ["data0", "data1", "data2"] for pos in [10, 11]],
        names=["dataset", "transcript_id", "transcript_position"],
    )
    labels = pd.Series([0, 1, 0, 1, 0, 1], index=index)
    metadata = pd.DataFrame({"group_id": ["human:A"] * 4 + ["synthetic:B"] * 2}, index=index)
    result = helpers["dataset_metrics"](labels, [.1, .9, .6, .4, .2, .7], metadata, .5).set_index("dataset")
    assert result.loc["data0", "average_precision"] == 1
    assert result.loc["data1", "average_precision"] == .5
    assert result.loc["data2", "f1"] == 1
    assert result.threshold.eq(.5).all() and result.sites.eq(2).all()
    with pytest.raises(ValueError, match="ordered index"):
        helpers["dataset_metrics"](labels, np.zeros(6), metadata.iloc[::-1])


def test_stratified_bootstrap_keeps_single_construct_present(helpers):
    labels = pd.Series([0, 0, 0, 0, 1, 1])
    groups = pd.Series(["A", "A", "B", "B", "C", "C"])
    strata = pd.Series(["human"] * 4 + ["synthetic"] * 2)
    # Without stratification, omitting the only synthetic group loses a class.
    samples = helpers["gene_bootstrap_metrics"](
        labels, [.1, .2, .3, .4, .8, .9], groups, .5,
        n_bootstrap=40, strata=strata,
    )
    assert len(samples) == 40
    assert samples.f1.eq(1).all()


@pytest.fixture(scope="module")
def cv_helpers():
    import warnings
    from datetime import datetime
    from time import perf_counter
    from zoneinfo import ZoneInfo
    from sklearn.compose import ColumnTransformer
    from sklearn.dummy import DummyClassifier
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from threadpoolctl import threadpool_limits
    from xgboost import XGBClassifier

    namespace = {**globals(), **locals(), "CV_FOLDS": 3, "SEED": 42,
                 "MODEL_THREADS": 1, "RUN_ID": "unit-test",
                 "LOCAL_TZ": ZoneInfo("Asia/Kuala_Lumpur"),
                 "FEATURE_SETS": {"synthetic": ["signal_a", "signal_b"]}}
    needed = {"make_group_folds", "build_experiment_model", "run_cv_experiment",
              "rank_experiments", "append_csv_records", "ranking_metrics", "positive_scores"}
    notebook = json.loads((Path(__file__).resolve().parents[1]
                          / "notebooks/data0_modelling_pipeline.ipynb").read_text())
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            definitions = [node for node in ast.parse("".join(cell["source"])).body
                           if isinstance(node, ast.FunctionDef) and node.name in needed]
            if definitions:
                exec(compile(ast.Module(body=definitions, type_ignores=[]),
                             "<notebook-cv>", "exec"), namespace)
    assert needed <= namespace.keys()
    return namespace


@pytest.fixture
def grouped_samples():
    rng = np.random.default_rng(11)
    genes = pd.Series(np.repeat([f"gene_{i}" for i in range(18)], 4), name="gene_id")
    labels = pd.Series(np.tile([0, 0, 0, 1], 18), name="label", dtype="int8")
    X = pd.DataFrame({"signal_a": rng.normal(size=len(genes)) + labels,
                      "signal_b": np.repeat(np.arange(18) * 10, 4) + rng.normal(size=len(genes)),
                      "excluded_column": np.arange(len(genes))})
    return X, labels, genes


def test_group_cv_is_disjoint_complete_and_reproducible(cv_helpers, grouped_samples):
    X, y, genes = grouped_samples
    folds = cv_helpers["make_group_folds"](X, y, genes)
    again = cv_helpers["make_group_folds"](X, y, genes)
    seen = np.zeros(len(y), dtype=int)
    for (fit_idx, assess_idx), (fit_again, assess_again) in zip(folds, again):
        assert set(genes.iloc[fit_idx]).isdisjoint(genes.iloc[assess_idx])
        np.testing.assert_array_equal(fit_idx, fit_again)
        np.testing.assert_array_equal(assess_idx, assess_again)
        seen[assess_idx] += 1
    np.testing.assert_array_equal(seen, np.ones(len(y), dtype=int))
    with pytest.raises(ValueError, match="ordered index"):
        cv_helpers["make_group_folds"](X, y.iloc[::-1], genes)


def test_xgboost_weight_uses_only_supplied_fit_labels(cv_helpers):
    spec = {"family": "xgb", "feature_set": "synthetic", "weighting": "ratio",
            "params": {"n_estimators": 2, "max_depth": 1}}
    model, details = cv_helpers["build_experiment_model"](spec, [0] * 18 + [1] * 2)
    assert model.named_steps["classifier"].get_params()["scale_pos_weight"] == 9
    assert details["fit_negatives"] == 18 and details["fit_positives"] == 2
    model2, _ = cv_helpers["build_experiment_model"](spec, [0] * 6 + [1] * 2)
    assert model2.named_steps["classifier"].get_params()["scale_pos_weight"] == 3
    sqrt_model, _ = cv_helpers["build_experiment_model"](
        {**spec, "weighting": "sqrt_ratio"}, [0] * 18 + [1] * 2)
    assert sqrt_model.named_steps["classifier"].get_params()["scale_pos_weight"] == 3
    with pytest.raises(ValueError, match="exactly classes"):
        cv_helpers["build_experiment_model"](spec, [0, 0])


def test_cv_scaler_is_fit_per_fold_and_reports_sample_sd(cv_helpers, grouped_samples, tmp_path):
    X, y, genes = grouped_samples
    folds = cv_helpers["make_group_folds"](X, y, genes)
    spec = {"experiment_id": "unit_lr", "family": "logistic", "feature_set": "synthetic",
            "weighting": "balanced", "params": {"C": 1.0},
            "parent": None, "stage": "test", "changed": "unit test"}
    original_factory = cv_helpers["build_experiment_model"]
    fitted_models = []

    def record_model(specification, fit_labels):
        model, details = original_factory(specification, fit_labels)
        fitted_models.append(model)
        return model, details

    cv_helpers["build_experiment_model"] = record_model
    try:
        summary, rows, oof = cv_helpers["run_cv_experiment"](spec, X, y, genes, folds, tmp_path)
    finally:
        cv_helpers["build_experiment_model"] = original_factory
    assert np.isfinite(oof).all()
    for (fit_idx, assess_idx), model, row in zip(folds, fitted_models, rows.to_dict("records")):
        np.testing.assert_allclose(model.named_steps["scale"].mean_,
                                   X.iloc[fit_idx][["signal_a", "signal_b"]].mean())
        assert model.named_steps["scale"].n_samples_seen_ == len(fit_idx)
        assert model.named_steps["classifier"].n_features_in_ == 2
        assert row["average_precision"] == pytest.approx(average_precision_score(y.iloc[assess_idx], oof[assess_idx]))
    for metric in ["roc_auc", "pr_auc_trapezoid", "average_precision"]:
        assert summary[metric + "_mean"] == pytest.approx(rows[metric].mean())
        assert summary[metric + "_sd"] == pytest.approx(rows[metric].std(ddof=1))
    assert len(pd.read_csv(tmp_path / "fold_results.csv")) == len(folds)
    restored = pd.read_parquet(tmp_path / "oof_unit_lr.parquet")
    np.testing.assert_allclose(restored.oof_score, oof)


def test_cv_selection_uses_ap_then_roc_and_stable_id(cv_helpers):
    table = pd.DataFrame({"experiment_id": ["D", "B", "A", "C"],
                          "average_precision_mean": [.5, .6, .6, .6],
                          "roc_auc_mean": [.99, .8, .8, .7]})
    assert cv_helpers["rank_experiments"](table).experiment_id.tolist() == ["A", "B", "C", "D"]


def test_experiment_logging_appends_history(cv_helpers, tmp_path):
    path = tmp_path / "history.csv"
    cv_helpers["append_csv_records"]([{"run_id": "first", "roc_auc": .8, "pr_auc": .4}], path)
    cv_helpers["append_csv_records"]([{"run_id": "second", "roc_auc": .9, "pr_auc": .5}], path)
    records = pd.read_csv(path)
    assert records.run_id.tolist() == ["first", "second"]
    assert records.pr_auc.tolist() == [.4, .5]
