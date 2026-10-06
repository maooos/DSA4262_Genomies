"""Check score alignment, reusable inputs, and failures before publishing CSVs."""

import gzip
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pytest
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline

from src.data_parser import parse_to_parquet
from src.feature_engineering import RAW_FEATURE_NAMES, build_features
from src.inference import KEYS, features_from_reads, predict_to_csv


@pytest.fixture
def inputs(tmp_path):
    records, pieces = [], []
    for i, n_reads in enumerate((1, 7, 2, 5, 3, 9)):
        values = np.random.default_rng(i).uniform(0.1, 2, (n_reads, 9))
        records.append({f"tx{i}": {str(i + 100): {"AGAACAT": values.tolist()}}})
        reads = pd.DataFrame(values, columns=RAW_FEATURE_NAMES)
        reads[KEYS[0]], reads[KEYS[1]] = f"tx{i}", i + 100
        reads["sequence"], reads["n_reads"] = "AGAACAT", n_reads
        pieces.append(reads)
    raw = tmp_path / "sites.json.gz"
    with gzip.open(raw, "wt") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")
    reads = pd.concat(pieces, ignore_index=True)
    # Independent, full-table construction matching the training notebook.
    features = build_features(reads).set_index(KEYS).join(
        reads.groupby(KEYS)[RAW_FEATURE_NAMES].mean().add_prefix("raw_mean__")
    )
    table = tmp_path / "features.parquet"
    features.to_parquet(table)
    parsed = tmp_path / "reads.parquet"
    parse_to_parquet(raw, parsed, batch_size=3)
    columns = list(features.columns)
    model = make_pipeline(
        ColumnTransformer([("features", "passthrough", ["n_reads", "central_log_dwell_mean"])]),
        LogisticRegression(),
    ).fit(features, [0, 0, 1, 1, 0, 1])
    artifact = tmp_path / "model.joblib"
    joblib.dump({"model": model, "feature_columns": columns, "model_name": "test model"}, artifact)
    return raw, parsed, table, artifact, features, model


@pytest.mark.parametrize("kind,index", [("reads", 0), ("reads", 1), ("features", 2)])
def test_formats_and_batch_boundaries_match_full_table(inputs, tmp_path, kind, index):
    *_, artifact, features, model = inputs
    output = tmp_path / "scores.csv"
    record = predict_to_csv(artifact, inputs[index], output, input_kind=kind, batch_size=3, threads=1)
    result = pd.read_csv(output).set_index(KEYS)
    pd.testing.assert_index_equal(result.index, features.index)
    np.testing.assert_allclose(result.score, model.predict_proba(features)[:, 1], rtol=1e-12, atol=1e-12)
    assert record["site_count"] == 6
    assert list(pd.read_csv(output)) == [*KEYS, "score"]
    assert json.loads(output.with_suffix(".metadata.json").read_text())["model_path"] == str(artifact)


def test_switch_model_and_csv_dataset(inputs, tmp_path):
    _, _, _, artifact, features, model = inputs
    # A second fitted model and a reordered subset must use the chosen model/site IDs.
    model.fit(features, [1, 1, 0, 0, 1, 0])
    second = tmp_path / "second.joblib"
    joblib.dump({"model": model, "feature_columns": list(features.columns)}, second)
    subset = features.iloc[[4, 0, 2]]
    path = tmp_path / "subset.csv"
    subset.to_csv(path)
    output = tmp_path / "second_scores.csv"
    predict_to_csv(second, path, output, batch_size=2, threads=1)
    result = pd.read_csv(output).set_index(KEYS)
    pd.testing.assert_index_equal(result.index, subset.index)
    np.testing.assert_allclose(result.score, model.predict_proba(subset)[:, 1])
    assert not np.allclose(result.score, joblib.load(artifact)["model"].predict_proba(subset)[:, 1])


@pytest.mark.parametrize("problem", ["duplicate", "missing", "nonfinite", "empty"])
def test_bad_input_does_not_leave_completed_output(inputs, tmp_path, problem):
    _, _, _, artifact, features, _ = inputs
    if problem == "duplicate":
        features = pd.concat([features, features.iloc[:1]])
    elif problem == "missing":
        features = features.drop(columns="n_reads")
    elif problem == "nonfinite":
        features.iloc[-1, 0] = np.nan
    else:
        features = features.iloc[:0]
    bad = tmp_path / "bad.parquet"
    features.to_parquet(bad)
    output = tmp_path / "bad.csv"
    with pytest.raises(ValueError):
        predict_to_csv(artifact, bad, output, batch_size=2, threads=1)
    assert not output.exists()
    assert not output.with_suffix(".metadata.json").exists()
    assert not list(tmp_path.glob(".*.partial"))


def test_existing_output_is_preserved(inputs, tmp_path):
    output = tmp_path / "existing.csv"
    output.write_text("keep me")
    with pytest.raises(FileExistsError):
        predict_to_csv(inputs[3], inputs[2], output)
    assert output.read_text() == "keep me"


def test_pooled_features_keep_overlapping_sites_separate(inputs, tmp_path):
    _, _, _, artifact, features, model = inputs
    pooled = pd.concat({"data0": features, "data1": features}, names=["dataset"])
    path = tmp_path / "pooled.parquet"
    pooled.to_parquet(path)
    output = tmp_path / "pooled.csv"
    metadata = predict_to_csv(artifact, path, output, batch_size=3, threads=1)
    keys = ["dataset", *KEYS]
    actual = pd.read_csv(output).set_index(keys)
    pd.testing.assert_index_equal(actual.index, pooled.index)
    np.testing.assert_allclose(actual.score, model.predict_proba(pooled)[:, 1])
    assert metadata["identifier_columns"] == keys
    assert len(actual) == 2 * len(features)


def test_incomplete_parquet_site_is_rejected(inputs):
    reads = pd.read_parquet(inputs[1]).iloc[:-1]
    with pytest.raises(ValueError, match="Incomplete"):
        features_from_reads(reads)


def test_actual_h04_matches_historical_scores_when_available(tmp_path):
    run = Path(__file__).resolve().parents[1] / (
        "data/processed/modelling/data0_gene_split_seed42/person_c_runs/20261005T103445_98fad904"
    )
    if not (run / "selected_model.joblib").exists():
        pytest.skip("Local model/data are excluded from Git")
    output = tmp_path / "H04.csv"
    predict_to_csv(run / "selected_model.joblib", run / "X_test.parquet", output, threads=1)
    expected = pd.read_parquet(run / "test_predictions.parquet")["score"]
    actual = pd.read_csv(output).set_index(KEYS)["score"]
    pd.testing.assert_index_equal(actual.index, expected.index)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
