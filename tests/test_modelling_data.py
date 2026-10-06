"""Small synthetic inputs only: pooled identity, leakage and batch boundaries."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scripts.run_modelling import run_notebook
from src.feature_engineering import RAW_FEATURE_NAMES
from src.modelling_data import (
    OBSERVATION_KEYS, assign_group_splits, build_pooled_features,
    load_site_metadata, validate_splits,
)


@pytest.fixture
def pooled_paths(tmp_path):
    paths = {}
    (tmp_path / "annotated").mkdir()
    for dataset in ("data0", "data1", "data2"):
        rows = []
        for group in range(4 if dataset == "data2" else 16):
            for label in (0, 1):
                synthetic = dataset == "data2"
                for read_idx in range(3):
                    row = {"transcript_id": f"mixture_{label}" if synthetic else f"TX{group}",
                           "transcript_position": group + 10 if synthetic else label + 10,
                           "sequence": "AAGACCA", "n_reads": 3, "read_idx": read_idx,
                           "gene_id": f"construct {group}" if synthetic else f"ENSG{group:06}",
                           "label": label}
                    row.update({name: 0.1 + group + read_idx + int(dataset[-1]) for name in RAW_FEATURE_NAMES})
                    if synthetic:
                        row.update(geo_reference_id=f"cc{group}", geo_reference_center_0based=group + 100,
                                   reference_mapping_status="inferred_unique_ordered_alignment",
                                   label_original=label * 0.25)
                    rows.append(row)
        paths[dataset] = tmp_path / "annotated" / f"{dataset}_reads.parquet"
        pd.DataFrame(rows).to_parquet(paths[dataset], index=False)
    return paths


def test_shared_genes_and_all_construct_mixtures_stay_together(pooled_paths):
    sites = load_site_metadata(pooled_paths)
    sites["split"] = assign_group_splits(sites)
    report = validate_splits(sites)
    assert sites.groupby("group_id").split.nunique().eq(1).all()
    human = sites.query("dataset != 'data2'")
    assert human.groupby("gene_id").dataset.nunique().eq(2).all()
    assert human.groupby("gene_id").split.nunique().eq(1).all()
    # Each anonymized mixture spans constructs and therefore spans partitions.
    synthetic = sites.query("dataset == 'data2'")
    assert synthetic.groupby("transcript_id").split.nunique().eq(3).all()
    assert report.loc["data2", "groups"].tolist() == [2, 1, 1]
    assert report[["positives", "negatives"]].gt(0).all().all()
    shuffled = sites.sample(frac=1, random_state=8)
    pd.testing.assert_series_equal(assign_group_splits(shuffled).sort_index(), sites["split"])
    # Distinct observations at an overlapping site are not dropped or merged.
    assert len(sites) == 16 * 2 * 2 + 4 * 2
    assert sites.set_index(OBSERVATION_KEYS).index.is_unique


def test_batched_features_keep_dataset_identity_and_match_whole_input(pooled_paths):
    sites = load_site_metadata(pooled_paths)
    tiny_batches = build_pooled_features(pooled_paths, sites, batch_size=2)
    whole_batches = build_pooled_features(pooled_paths, sites, batch_size=10000)
    pd.testing.assert_frame_equal(tiny_batches, whole_batches)
    assert tiny_batches.index.names == OBSERVATION_KEYS
    assert tiny_batches.shape == (len(sites), 162)
    assert np.isfinite(tiny_batches.to_numpy()).all()
    signal = "raw_mean__" + RAW_FEATURE_NAMES[0]
    first = tiny_batches.loc[("data0", "TX0", 10), signal]
    second = tiny_batches.loc[("data1", "TX0", 10), signal]
    assert second - first == pytest.approx(1)
    assert set(tiny_batches.columns).isdisjoint({"gene_id", "group_id", "label", "label_original", "dataset", "geo_reference_id"})


def test_cross_dataset_gene_conflict_is_rejected(pooled_paths):
    frame = pd.read_parquet(pooled_paths["data1"])
    frame.loc[frame.transcript_id.eq("TX0"), "gene_id"] = "wrong_gene"
    frame.to_parquet(pooled_paths["data1"], index=False)
    with pytest.raises(ValueError, match="inconsistent genes across"):
        load_site_metadata(pooled_paths)


@pytest.mark.parametrize("problem", ["fractional_label", "wrong_binary_label", "missing_construct", "partial"])
def test_incomplete_or_stale_annotations_are_rejected(pooled_paths, problem):
    path = pooled_paths["data2"]
    frame = pd.read_parquet(path)
    if problem == "fractional_label":
        frame["label"] = frame.label_original
    elif problem == "wrong_binary_label":
        frame["label"] = 0
    elif problem == "missing_construct":
        frame = frame.loc[frame.geo_reference_id.ne("cc0")]
    else:
        path.with_suffix(".audit.json").write_text(json.dumps({"max_sites_requested": 2}))
    frame.to_parquet(path, index=False)
    with pytest.raises(ValueError):
        load_site_metadata(pooled_paths)


def test_leaked_human_observation_is_rejected(pooled_paths):
    sites = load_site_metadata(pooled_paths)
    sites["split"] = assign_group_splits(sites)
    i = sites.index[sites.dataset.eq("data1")][0]
    sites.loc[i, "split"] = "test" if sites.loc[i, "split"] != "test" else "train"
    with pytest.raises(ValueError, match="Group leakage"):
        validate_splits(sites)


def test_inconsistent_repeated_read_annotation_is_rejected(pooled_paths):
    sites = load_site_metadata(pooled_paths)
    frame = pd.read_parquet(pooled_paths["data2"])
    frame.loc[1, "geo_reference_id"] = "wrong_construct"
    frame.to_parquet(pooled_paths["data2"], index=False)
    with pytest.raises(ValueError, match="inconsistent annotations"):
        build_pooled_features(pooled_paths, sites, batch_size=10)


def test_cli_preparation_boundary_cannot_execute_training(tmp_path):
    notebook = tmp_path / "notebook.ipynb"
    notebook.write_text(json.dumps({"cells": [
        {"cell_type": "code", "source": ["prepared = True"], "metadata": {"tags": ["prepared-inputs"]}},
        {"cell_type": "code", "source": ["raise AssertionError('must not train')"], "metadata": {}},
    ]}))
    result = run_notebook(notebook, prepare_only=True)
    assert result["prepared"]
    notebook.write_text(json.dumps({"cells": [
        {"cell_type": "code", "source": ["raise AssertionError('must not run')"], "metadata": {}},
    ]}))
    with pytest.raises(ValueError, match="boundary"):
        run_notebook(notebook, prepare_only=True)


def test_actual_notebook_prepares_pooled_fixture_without_model_fitting(pooled_paths):
    notebook = Path(__file__).resolve().parents[1] / "notebooks/data0_modelling_pipeline.ipynb"
    result = run_notebook(notebook, prepare_only=True, namespace={
        "PROCESSED_DIR": pooled_paths["data0"].parent.parent, "MODEL_THREADS": 1,
    })
    run_dir = result["OUTPUT_DIR"]
    assert (run_dir / "split_manifest.parquet").exists()
    assert (run_dir / "split_dataset_summary.csv").exists()
    assert not (run_dir / "cv_started.json").exists()
    assert not (run_dir / "selected_model.joblib").exists()
    for split in ["train", "val", "test"]:
        X = pd.read_parquet(run_dir / f"X_{split}.parquet")
        y = pd.read_parquet(run_dir / f"y_{split}.parquet")
        meta = pd.read_parquet(run_dir / f"metadata_{split}.parquet")
        assert set(X.index.get_level_values("dataset")) == {"data0", "data1", "data2"}
        assert X.index.equals(y.index) and X.index.equals(meta.index)
