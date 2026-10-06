"""Synthetic fixtures only; exercise all four notebooks without project-data fits."""

import ast
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

torch = pytest.importorskip("torch", reason="Install requirements.txt for advanced model tests")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from AdvancedModel.common import data as data_module
from AdvancedModel.common import ensemble as ensemble_module
from AdvancedModel.common.data import (
    KmerNormalizer, ReadBags, inner_group_split, prepare_store,
)
from AdvancedModel.common.ensemble import blend_scores, select_oof_blend
from AdvancedModel.common.evaluation import metric_report, select_threshold
from AdvancedModel.common.models import MILModel, ModelConfig, site_loss
from AdvancedModel.common.training import TrainConfig, fit_partition, run_grouped_cv
from src.feature_engineering import RAW_FEATURE_NAMES

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = sorted((ROOT / "AdvancedModel").glob("*.ipynb"))


@pytest.fixture
def synthetic_store(tmp_path):
    annotated = tmp_path / "annotated"
    annotated.mkdir()
    rng = np.random.default_rng(10)
    for dataset in ["data0", "data1", "data2"]:
        rows = []
        for group in range(4 if dataset == "data2" else 24):
            for label in [0, 1]:
                synthetic = dataset == "data2"
                n_reads = 2 + (group + label) % 5
                for read in range(n_reads):
                    values = rng.uniform(.1, 2, 9) + group * .1 + label * .2
                    row = dict(zip(RAW_FEATURE_NAMES, values))
                    row.update(transcript_id=f"mixture_{label}" if synthetic else f"TX{group}",
                               transcript_position=group + 10 if synthetic else 10 + label,
                               gene_id=f"Synthetic cc{group}" if synthetic else f"ENSG{group}",
                               sequence="AAGACCA", n_reads=n_reads, label=label, read_idx=read)
                    if synthetic:
                        row.update(geo_reference_id=f"cc{group}", geo_reference_center_0based=group + 100,
                                   reference_mapping_status="inferred_unique_ordered_alignment",
                                   label_original=float(label) * .25)
                    rows.append(row)
        pd.DataFrame(rows).to_parquet(annotated / f"{dataset}_reads.parquet", index=False)
    return prepare_store(tmp_path, seed=42, cv_folds=3, batch_size=17)


def test_cache_reuse_keeps_every_read_and_separates_dataset_observations(synthetic_store):
    store = synthetic_store
    before = (store.cache_dir / "reads.npy").stat().st_mtime_ns
    again = prepare_store(store.cache_dir.parents[2], cv_folds=3)
    assert again.cache_dir == store.cache_dir
    assert (store.cache_dir / "reads.npy").stat().st_mtime_ns == before
    assert len(store.reads) == store.sites.n_reads.sum()
    assert store.sites.groupby("group_id").split.nunique().eq(1).all()
    assert store.sites.query("split == 'train'").groupby("group_id").assessment_fold.nunique().eq(1).all()
    assert store.sites.groupby("dataset").size().to_dict() == {"data0": 48, "data1": 48, "data2": 8}
    for row_id in [0, 3, len(store.sites) - 1]:
        np.testing.assert_allclose(store.raw_bag(row_id).sum(0), store.sums[row_id], rtol=1e-6)


def test_normalizer_uses_only_fitting_sites_and_falls_back_for_unseen_kmers(synthetic_store):
    store = synthetic_store
    fit, stop = inner_group_split(store, store.ids("train"))
    expected = KmerNormalizer.fit(store, fit)
    original_sums, original_squares = store.sums, store.squares
    store.sums, store.squares = store.sums.copy(), store.squares.copy()
    excluded = np.setdiff1d(np.arange(len(store.sites)), fit)
    store.sums[excluded] += 1e8
    store.squares[excluded] += 1e16
    actual = KmerNormalizer.fit(store, fit)
    np.testing.assert_array_equal(expected.mean, actual.mean)
    np.testing.assert_array_equal(expected.scale, actual.scale)
    store.sums, store.squares = original_sums, original_squares
    assert not set(store.sites.iloc[fit].group_id) & set(store.sites.iloc[stop].group_id)
    assert np.isfinite(expected.transform(store.raw_bag(0), np.array([1023, 1023, 1023]))).all()
    with pytest.raises(ValueError, match="outer-training"):
        KmerNormalizer.fit(store, store.ids("test"))


def test_bag_sampling_is_reproducible_and_short_bags_are_masked(synthetic_store):
    store = synthetic_store
    normalizer = KmerNormalizer.fit(store, store.ids("train"))
    small = int(np.flatnonzero(store.sites.n_reads.eq(2))[0])
    bags = ReadBags(store, [small], normalizer, max_reads=8, seed=23)
    a, b = bags[0], bags[0]
    torch.testing.assert_close(a["features"], b["features"])
    assert a["mask"].sum() == 2
    assert (a["features"][~a["mask"]] == 0).all()
    long_ids = np.flatnonzero(store.sites.n_reads.eq(6))
    a = ReadBags(store, long_ids, normalizer, max_reads=3, seed=23, draw=0)
    b = ReadBags(store, long_ids[::-1], normalizer, max_reads=3, seed=23, draw=0)
    torch.testing.assert_close(a[0]["features"], b[len(b) - 1]["features"])
    c = ReadBags(store, long_ids, normalizer, max_reads=3, seed=23, draw=1)
    assert any(not torch.equal(a[i]["features"], c[i]["features"]) for i in range(len(a)))


@pytest.mark.parametrize("architecture", ["noisy_or", "gated_attention", "set_transformer"])
def test_network_is_permutation_invariant_padding_safe_and_differentiable(architecture):
    torch.manual_seed(7)
    model = MILModel(ModelConfig(architecture=architecture, hidden_dim=16, read_dim=8,
                                attention_dim=8, heads=2, dropout=0.)).eval()
    features = torch.randn(2, 6, 15)
    tokens = torch.tensor([[0, 1, 2], [3, 4, 5]])
    mask = torch.tensor([[1, 1, 0, 0, 0, 0], [1, 1, 1, 1, 1, 0]], dtype=torch.bool)
    original = model(features, tokens, mask)["logits"]
    permutation = torch.tensor([5, 2, 0, 4, 1, 3])
    reordered = model(features[:, permutation], tokens, mask[:, permutation])["logits"]
    torch.testing.assert_close(original, reordered, rtol=1e-5, atol=1e-6)
    features[~mask] = torch.nan
    padded = torch.cat([features, torch.full((2, 3, 15), float('nan'))], dim=1)
    padded_mask = torch.cat([mask, torch.zeros((2, 3), dtype=torch.bool)], dim=1)
    torch.testing.assert_close(original, model(padded, tokens, padded_mask)["logits"], rtol=1e-5, atol=1e-6)
    outputs = model(features, tokens, mask)
    loss = site_loss(outputs, {"label": torch.tensor([0., 1.]), "dataset": torch.tensor([0, 2]),
                               "fraction": torch.tensor([0., .25])}, auxiliary_weight=.1)
    loss.backward()
    gradients = [p.grad for p in model.parameters() if p.grad is not None]
    assert gradients and all(torch.isfinite(g).all() for g in gradients)
    assert any(g.abs().sum() > 0 for g in gradients)
    with pytest.raises(ValueError, match="at least one"):
        model(features, tokens, torch.zeros_like(mask))


def test_noisy_or_matches_probability_formula_and_ignores_padding():
    model = MILModel(ModelConfig(architecture="noisy_or", dropout=0.)).eval()
    torch.nn.init.zeros_(model.read_head.weight)
    torch.nn.init.constant_(model.read_head.bias, np.log(.1 / .9))
    mask = torch.tensor([[1, 1, 0], [1, 1, 1]], dtype=torch.bool)
    probabilities = model(torch.zeros(2, 3, 15), torch.zeros(2, 3, dtype=torch.long), mask)["logits"].sigmoid()
    torch.testing.assert_close(probabilities, torch.tensor([1 - .9 ** 2, 1 - .9 ** 3]), rtol=1e-5, atol=1e-6)


def test_val_test_cannot_be_used_as_fit_or_early_stop(synthetic_store):
    model = ModelConfig()
    config = TrainConfig(device="cpu", max_epochs=1, threads=1)
    with pytest.raises(ValueError, match="outer training"):
        fit_partition(synthetic_store, synthetic_store.ids("val"), model, config)
    with pytest.raises(ValueError, match="outer train"):
        fit_partition(synthetic_store, synthetic_store.ids("train"), model, config,
                      stopping_ids=synthetic_store.ids("test"))


def test_blend_rejects_misalignment_leakage_and_can_choose_single_branch():
    frame = pd.DataFrame({"label": [0, 1, 0, 1], "group_id": ["A", "A", "B", "B"],
                          "assessment_fold": [1, 1, 2, 2], "split": "train", "score": [.1, .9, .2, .8]})
    bad_model = frame.assign(score=[.9, .1, .8, .2])
    weight, _ = select_oof_blend(frame, bad_model, [0, .5, 1])
    assert weight == 1
    with pytest.raises(ValueError, match="identities"):
        select_oof_blend(frame, bad_model.iloc[::-1])
    with pytest.raises(ValueError, match="metadata mismatch"):
        select_oof_blend(frame, bad_model.assign(assessment_fold=[2, 2, 1, 1]))
    with pytest.raises(ValueError, match="outer-training"):
        select_oof_blend(frame.assign(split="test"), bad_model.assign(split="test"))
    with pytest.raises(ValueError, match="weight"):
        blend_scores([.1], [.2], 1.5)


def test_dataset_reports_include_prevalence_reference_and_shared_threshold():
    meta = pd.DataFrame({"dataset": ["data0"] * 4 + ["data2"] * 4,
                         "group_id": ["human:A"] * 4 + ["synthetic:A"] * 4,
                         "label": [0, 0, 0, 1, 0, 1, 1, 1]})
    scores = np.array([.1, .2, .4, .8, .3, .7, .8, .9])
    threshold, _ = select_threshold(meta.label, scores)
    report = metric_report(meta, scores, threshold).set_index("dataset")
    assert report.loc["data2", "prior_ap_reference"] == .75
    assert report.loc["data2", "always_positive_f1"] == pytest.approx(6 / 7)
    assert report.threshold.eq(threshold).all()


def test_notebooks_are_separate_valid_python_and_have_no_saved_outputs():
    assert len(NOTEBOOKS) == 4
    for path in NOTEBOOKS:
        notebook = json.loads(path.read_text())
        assert notebook["nbformat"] == 4 and notebook["nbformat_minor"] == 5
        assert len({cell["id"] for cell in notebook["cells"]}) == len(notebook["cells"])
        code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]
        for cell in code_cells:
            ast.parse("".join(cell["source"]))
            assert cell["execution_count"] is None and cell["outputs"] == []
        source = "\n".join("".join(cell["source"]) for cell in code_cells)
        assert "run_grouped_cv" in source and "save_evaluation" in source and "mark_complete" in source


@pytest.mark.parametrize("notebook_path", NOTEBOOKS, ids=lambda p: p.stem)
def test_each_notebook_runs_independently_end_to_end_on_tiny_fixture(synthetic_store, notebook_path, monkeypatch):
    """Execute every real code cell, scaling only data/epochs/width for this test."""
    monkeypatch.setattr(data_module, "prepare_store", lambda *args, **kwargs: synthetic_store)
    monkeypatch.setattr(ensemble_module, "H04_PARAMETERS", {
        **ensemble_module.H04_PARAMETERS, "max_iter": 3, "min_samples_leaf": 2,
    })
    namespace = {"__name__": "__main__"}
    for i, cell in enumerate(json.loads(notebook_path.read_text())["cells"]):
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        exec(compile(source, f"{notebook_path.name}:cell{i}", "exec"), namespace)
        if "NOTEBOOK_PATH =" in source:
            namespace["MODEL"] = replace(namespace["MODEL"], hidden_dim=16, read_dim=8,
                                         attention_dim=8, heads=2, dropout=0.)
            namespace["TRAIN"] = replace(namespace["TRAIN"], device="cpu", threads=1, max_epochs=1,
                                         patience=1, batch_size=32, max_reads=3, predict_draws=1)
            namespace["BOOTSTRAP_REPEATS"] = 20
        plt.close("all")
    run = namespace["RUN_DIR"]
    assert (run / "complete.json").exists()
    assert (run / "model.pt").exists()
    assert len(pd.read_csv(run / "test_metrics_by_dataset.csv")) == 4
    cv = namespace["cv"]
    assert len(cv.fold_metrics) == 3
    assert cv.oof_predictions.index.equals(synthetic_store.indexed_metadata(synthetic_store.ids("train")).index)
    for fold in [1, 2, 3]:
        inner = pd.read_parquet(run / f"fold_{fold}_inner_groups.parquet")
        fit = set(inner.loc[inner.role.eq("fit"), "group_id"])
        stop = set(inner.loc[inner.role.eq("early_stop"), "group_id"])
        assessment = set(synthetic_store.sites.loc[synthetic_store.sites.assessment_fold.eq(fold), "group_id"])
        assert fit.isdisjoint(stop) and (fit | stop).isdisjoint(assessment)
    if "ensemble" in notebook_path.name:
        assert (run / "ensemble.json").exists() and (run / "h04.joblib").exists()
        assert (run / "test_component_ranking.csv").exists()
    with pytest.raises(FileExistsError):
        run_grouped_cv(synthetic_store, namespace["MODEL"], namespace["TRAIN"], run)
