"""Submission coverage must span all partitions and keep datasets separate."""

import json

import joblib
import numpy as np
import pandas as pd
import pytest
from sklearn.linear_model import LogisticRegression

from scripts.generate_submissions import generate_submissions
from src.modelling_data import DATASETS, OBSERVATION_KEYS, SPLIT_NAMES


@pytest.fixture
def submission_inputs(tmp_path):
    run, raw = tmp_path / "run", tmp_path / "raw"
    run.mkdir()
    # Identical site IDs across datasets, with different feature values.
    index = pd.MultiIndex.from_product(
        [DATASETS, ["shared_tx"], [100, 200, 300]], names=OBSERVATION_KEYS,
    )
    features = pd.DataFrame({"signal": np.arange(9, dtype=float)}, index=index)
    model = LogisticRegression().fit(features, [0, 0, 0, 0, 1, 1, 1, 1, 1])
    joblib.dump({
        "model": model, "feature_columns": ["signal"],
        "experiment": {"experiment_id": "H04"}, "datasets": list(DATASETS),
    }, run / "selected_model.joblib")
    for split, position in zip(SPLIT_NAMES, [100, 200, 300]):
        features.loc[index.get_level_values("transcript_position") == position].to_parquet(run / f"X_{split}.parquet")
    for dataset in DATASETS:
        directory = raw / dataset
        directory.mkdir(parents=True)
        name = "data.info.labelled" if dataset == "data0" else "data.info"
        # Metadata order differs from partition order; label values are unused.
        features.xs(dataset).iloc[[2, 0, 1]][[]].assign(label=1).to_csv(directory / name)
    return run, raw, tmp_path / "submission", features, model


def test_all_partitions_one_model_and_metadata_order(submission_inputs):
    run, raw, output, features, model = submission_inputs
    manifest = generate_submissions(run, raw_dir=raw, output_dir=output, threads=1)
    assert sorted(path.name for path in output.glob("*.csv")) == [
        f"Genomies_dataset{i}.csv" for i in range(3)
    ]
    for i, dataset in enumerate(DATASETS):
        saved = pd.read_csv(output / f"Genomies_dataset{i}.csv")
        assert saved.columns.tolist() == ["transcript_id", "transcript_position", "score"]
        assert saved.transcript_position.tolist() == [300, 100, 200]
        expected = model.predict_proba(features.xs(dataset).iloc[[2, 0, 1]])[:, 1]
        np.testing.assert_allclose(saved.score, expected, rtol=1e-14)
    assert [entry["site_count"] for entry in manifest["files"]] == [3, 3, 3]
    assert json.loads((output / "Genomies_submission_metadata.json").read_text()) == manifest


@pytest.mark.parametrize("problem", ["missing", "duplicate", "wrong_model"])
def test_invalid_run_does_not_publish_partial_submission(submission_inputs, problem):
    run, raw, output, features, _ = submission_inputs
    if problem == "missing":
        # A failure on the final dataset must leave no earlier published CSVs.
        frame = pd.read_parquet(run / "X_test.parquet").drop(index=("data2", "shared_tx", 300))
        frame.to_parquet(run / "X_test.parquet")
    elif problem == "duplicate":
        pd.concat([pd.read_parquet(run / "X_test.parquet"), features.iloc[:1]]).to_parquet(run / "X_test.parquet")
    else:
        bundle = joblib.load(run / "selected_model.joblib")
        bundle["experiment"]["experiment_id"] = "H03"
        joblib.dump(bundle, run / "selected_model.joblib")
    with pytest.raises(ValueError):
        generate_submissions(run, raw_dir=raw, output_dir=output, threads=1)
    assert not output.exists()


def test_existing_submission_is_preserved(submission_inputs):
    run, raw, output, *_ = submission_inputs
    output.mkdir()
    existing = output / "Genomies_dataset1.csv"
    existing.write_text("existing submission\n")
    with pytest.raises(FileExistsError):
        generate_submissions(run, raw_dir=raw, output_dir=output, threads=1)
    assert existing.read_text() == "existing submission\n"
    assert len(list(output.iterdir())) == 1
