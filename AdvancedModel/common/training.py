"""Grouped CV with inner early stopping, full-fold refits and portable weights."""

from dataclasses import asdict, dataclass
from datetime import datetime
import hashlib
import json
from pathlib import Path
import random
from time import perf_counter
from uuid import uuid4
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd
import sklearn
import torch
from torch.utils.data import DataLoader

from src.data_parser import file_sha256
from .data import KmerNormalizer, ReadBags, frame_digest, inner_group_split
from .evaluation import metric_report, ranking_metrics
from .models import MILModel, ModelConfig, site_loss


@dataclass
class TrainConfig:
    seed: int = 42
    device: str = "auto"
    threads: int = 4
    batch_size: int = 128
    max_reads: int = 32
    max_epochs: int = 25
    patience: int = 4
    learning_rate: float = 1e-3
    weight_decay: float = 1e-4
    gradient_clip: float = 5.
    predict_draws: int = 5
    inner_stop_fraction: float = .15
    dataset_weights: tuple = (1., 1., 1.)
    auxiliary_weight: float = 0.
    positive_weight: str = "none"

    def __post_init__(self):
        if min(self.threads, self.batch_size, self.max_reads, self.max_epochs, self.patience, self.predict_draws) < 1:
            raise ValueError("Training counts must be positive.")
        if len(self.dataset_weights) != 3 or not np.isfinite(self.dataset_weights).all() or min(self.dataset_weights) <= 0:
            raise ValueError("Provide three positive finite dataset weights.")
        if self.positive_weight not in {"none", "ratio", "sqrt_ratio"}:
            raise ValueError("Unknown positive class weighting strategy.")
        if not 0 < self.inner_stop_fraction < 1 or self.learning_rate <= 0 or self.gradient_clip <= 0:
            raise ValueError("Invalid optimization/stopping configuration.")
        if self.auxiliary_weight < 0 or self.weight_decay < 0:
            raise ValueError("Loss/decay weights must be nonnegative.")


def resolve_device(choice="auto"):
    if choice == "auto":
        choice = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
    if choice not in {"cuda", "mps", "cpu"}:
        raise ValueError("device must be auto, cuda, mps or cpu.")
    if choice == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA requested but unavailable.")
    if choice == "mps" and not torch.backends.mps.is_available():
        raise ValueError("Apple MPS requested but unavailable.")
    return torch.device(choice)


def seed_everything(seed, threads):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.set_num_threads(threads)
    torch.use_deterministic_algorithms(True, warn_only=True)
    if torch.backends.cudnn.is_available():
        torch.backends.cudnn.benchmark = False


def create_run(store, model_config, train_config, *, name=None, notebook_path=None):
    name = name or model_config.architecture
    root = store.cache_dir.parent.parent / "runs" / name
    timestamp = datetime.now(ZoneInfo("Asia/Singapore")).strftime("%Y%m%dT%H%M%S")
    output = root / f"{timestamp}_{uuid4().hex[:8]}"
    output.mkdir(parents=True, exist_ok=False)
    provenance = {"run_name": name, "created_at": datetime.now(ZoneInfo("Asia/Singapore")).isoformat(),
                  "model": asdict(model_config), "training": asdict(train_config),
                  "device": str(resolve_device(train_config.device)),
                  "data_fingerprint": store.manifest["fingerprint"], "sources": store.manifest["sources"],
                  "split_seed": store.manifest["seed"], "cv_folds": store.manifest["cv_folds"],
                  "versions": {"torch": torch.__version__, "numpy": np.__version__, "sklearn": sklearn.__version__},
                  "code_sha256": {p.name: file_sha256(p) for p in Path(__file__).parent.glob("*.py")},
                  "selection": "Inner grouped stopping selects epochs; assessment folds only score full-fold refits.",
                  "final_epochs": "Median of training CV selected epochs; no validation/test fitting."}
    if notebook_path is not None:
        notebook = json.loads(Path(notebook_path).read_text())
        source = json.dumps([cell["source"] for cell in notebook["cells"]], sort_keys=True)
        provenance["notebook_sources_sha256"] = hashlib.sha256(source.encode()).hexdigest()
    (output / "config.json").write_text(json.dumps(provenance, indent=2))
    store.sites.to_parquet(output / "split_and_fold_manifest.parquet", index=False)
    return output


def move_batch(batch, device):
    return {key: value.to(device) for key, value in batch.items()}


def predict_bags(model, normalizer, store, ids, config, *, draws=None):
    """Deterministic repeated read subsampling, aligned exactly to ids."""
    ids = np.asarray(ids, dtype=int)
    if not len(ids):
        raise ValueError("No sites supplied for prediction.")
    device = next(model.parameters()).device
    draws = config.predict_draws if draws is None else draws
    if draws < 1:
        raise ValueError("Prediction draws must be positive.")
    scores = np.zeros(len(ids), dtype=np.float64)
    model.eval()
    with torch.inference_mode():
        for draw in range(draws):
            bags = ReadBags(store, ids, normalizer, config.max_reads, config.seed + 100_000, draw)
            loader = DataLoader(bags, batch_size=config.batch_size, shuffle=False, num_workers=0)
            cursor = 0
            for batch in loader:
                batch = move_batch(batch, device)
                result = model(batch["features"], batch["tokens"], batch["mask"])
                values = result["logits"].sigmoid().cpu().numpy()
                scores[cursor:cursor + len(values)] += values / draws
                cursor += len(values)
    if not np.isfinite(scores).all() or ((scores < 0) | (scores > 1)).any():
        raise ValueError("Model produced invalid probabilities.")
    return scores


def fit_partition(store, ids, model_config, config, *, epochs=None, stopping_ids=None, seed=None, log_path=None):
    ids = np.asarray(ids, dtype=int)
    if not len(ids) or not store.sites.iloc[ids].split.eq("train").all():
        raise ValueError("Only outer training sites may be used for fitting.")
    if store.sites.iloc[ids].label.nunique() != 2:
        raise ValueError("Fitting requires both classes.")
    if stopping_ids is not None:
        stopping_ids = np.asarray(stopping_ids, dtype=int)
        if not store.sites.iloc[stopping_ids].split.eq("train").all():
            raise ValueError("Stopping sites must come from outer train.")
        if set(store.sites.iloc[ids].group_id) & set(store.sites.iloc[stopping_ids].group_id):
            raise ValueError("Fitting/stopping groups overlap.")
    seed = config.seed if seed is None else seed
    seed_everything(seed, config.threads)
    device = resolve_device(config.device)
    normalizer = KmerNormalizer.fit(store, ids)
    model = MILModel(model_config).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.learning_rate, weight_decay=config.weight_decay)
    labels = store.sites.iloc[ids].label.to_numpy()
    ratio = float((labels == 0).sum() / (labels == 1).sum())
    positive_weight = {"none": 1., "ratio": ratio, "sqrt_ratio": np.sqrt(ratio)}[config.positive_weight]
    bags = ReadBags(store, ids, normalizer, config.max_reads, seed)
    generator = torch.Generator().manual_seed(seed)
    loader = DataLoader(bags, batch_size=config.batch_size, shuffle=True, generator=generator, num_workers=0)
    best_ap, best_epoch, stale, history = -np.inf, 1, 0, []
    for epoch in range(1, (config.max_epochs if epochs is None else epochs) + 1):
        start = perf_counter()
        bags.draw = epoch
        model.train()
        total_loss, count = 0., 0
        for batch in loader:
            batch = move_batch(batch, device)
            optimizer.zero_grad(set_to_none=True)
            result = model(batch["features"], batch["tokens"], batch["mask"])
            loss = site_loss(result, batch, config.dataset_weights, config.auxiliary_weight, positive_weight)
            if not torch.isfinite(loss):
                raise ValueError("Nonfinite loss; inspect signals and learning rate.")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), config.gradient_clip, error_if_nonfinite=True)
            optimizer.step()
            total_loss += float(loss.detach()) * len(batch["label"])
            count += len(batch["label"])
        row = {"epoch": epoch, "training_loss": total_loss / count, "fit_sites": count,
               "positive_weight": positive_weight, "seconds": perf_counter() - start}
        if stopping_ids is not None:
            scores = predict_bags(model, normalizer, store, stopping_ids, config, draws=1)
            row["inner_stop_ap"] = ranking_metrics(store.sites.iloc[stopping_ids].label, scores)["average_precision"]
            if row["inner_stop_ap"] > best_ap + 1e-6:
                best_ap, best_epoch, stale = row["inner_stop_ap"], epoch, 0
            else:
                stale += 1
        else:
            best_epoch = epoch
        history.append(row)
        if log_path is not None:
            pd.DataFrame([row]).to_csv(log_path, mode="a", header=not Path(log_path).exists(), index=False)
        suffix = f" | inner AP={row['inner_stop_ap']:.4f}" if stopping_ids is not None else ""
        print(f"Epoch {epoch:02d}: loss={row['training_loss']:.4f}{suffix}", flush=True)
        if stopping_ids is not None and stale >= config.patience:
            break
    # Inner-fit weights are used only to choose epoch count, then discarded.
    return model, normalizer, pd.DataFrame(history), best_epoch


@dataclass
class CVResult:
    fold_metrics: pd.DataFrame
    oof_predictions: pd.DataFrame
    best_epochs: list
    history: pd.DataFrame


def run_grouped_cv(store, model_config, config, run_dir):
    with (run_dir / "cv_started.json").open("x") as handle:
        json.dump({"started_at": datetime.now(ZoneInfo("Asia/Singapore")).isoformat()}, handle)
    train_ids = store.ids("train")
    oof = np.full(len(store.sites), np.nan)
    rows, epochs, histories = [], [], []
    fold_ids = store.sites.assessment_fold.to_numpy()
    for fold in range(1, store.manifest["cv_folds"] + 1):
        fit_ids, assessment_ids = train_ids[fold_ids[train_ids] != fold], train_ids[fold_ids[train_ids] == fold]
        if set(store.sites.iloc[fit_ids].group_id) & set(store.sites.iloc[assessment_ids].group_id):
            raise ValueError("Assessment groups leaked into fold fitting.")
        inner_fit, inner_stop = inner_group_split(store, fit_ids, config.seed + fold, config.inner_stop_fraction)
        membership = store.sites.iloc[np.r_[inner_fit, inner_stop]][["dataset", "transcript_id", "transcript_position", "group_id"]].copy()
        membership["role"] = ["fit"] * len(inner_fit) + ["early_stop"] * len(inner_stop)
        membership.to_parquet(run_dir / f"fold_{fold}_inner_groups.parquet", index=False)
        print(f"Fold {fold}: select epoch count using inner training groups", flush=True)
        tuning_model, _, history, selected_epoch = fit_partition(
            store, inner_fit, model_config, config, stopping_ids=inner_stop, seed=config.seed + fold,
            log_path=run_dir / f"fold_{fold}_inner_history.csv",
        )
        del tuning_model
        histories.append(history.assign(fold=fold, phase="inner_stop"))
        print(f"Fold {fold}: refit on all fitting groups for {selected_epoch} epochs", flush=True)
        model, normalizer, history, _ = fit_partition(
            store, fit_ids, model_config, config, epochs=selected_epoch, seed=config.seed + 100 + fold,
            log_path=run_dir / f"fold_{fold}_refit_history.csv",
        )
        scores = predict_bags(model, normalizer, store, assessment_ids, config)
        oof[assessment_ids] = scores
        row = {"fold": fold, "epochs": selected_epoch, "fit_sites": len(fit_ids), "assessment_sites": len(assessment_ids),
               "fit_groups": store.sites.iloc[fit_ids].group_id.nunique(),
               "assessment_groups": store.sites.iloc[assessment_ids].group_id.nunique(),
               **ranking_metrics(store.sites.iloc[assessment_ids].label, scores)}
        rows.append(row)
        epochs.append(selected_epoch)
        histories.append(history.assign(fold=fold, phase="refit"))
        pd.DataFrame(rows).to_csv(run_dir / "cv_fold_metrics.csv", index=False)
        print(f"Fold {fold}: held-out AP={row['average_precision']:.4f}", flush=True)
        del model
    if not np.isfinite(oof[train_ids]).all():
        raise ValueError("Incomplete OOF coverage.")
    predictions = store.indexed_metadata(train_ids).assign(score=oof[train_ids])
    predictions.to_parquet(run_dir / "oof_predictions.parquet")
    metric_report(store.sites.iloc[train_ids], oof[train_ids]).to_csv(run_dir / "cv_oof_metrics_by_dataset.csv", index=False)
    fold_metrics = pd.DataFrame(rows)
    fold_metrics[["average_precision", "roc_auc", "pr_auc_trapezoid"]].agg(["mean", "std"]).to_csv(run_dir / "cv_mean_sd.csv")
    combined_history = pd.concat(histories, ignore_index=True)
    combined_history.to_csv(run_dir / "training_history.csv", index=False)
    (run_dir / "cv_complete.json").write_text(json.dumps({"best_epochs": epochs, "median_epochs": int(np.median(epochs))}))
    return CVResult(fold_metrics, predictions, epochs, combined_history)


def fit_final(store, model_config, config, cv_result, run_dir):
    epochs = max(1, int(np.median(cv_result.best_epochs)))
    model, normalizer, history, _ = fit_partition(
        store, store.ids("train"), model_config, config, epochs=epochs,
        seed=config.seed + 1000, log_path=run_dir / "final_training_history.csv",
    )
    return model, normalizer, epochs


def save_mil_bundle(path, model, normalizer, config, store, *, threshold, epochs):
    torch.save({"format_version": 1, "model_config": asdict(model.config), "training_config": asdict(config),
                "state_dict": {k: v.detach().cpu() for k, v in model.state_dict().items()},
                "normalizer": normalizer.state_dict(), "threshold": float(threshold) if threshold is not None else None, "epochs": int(epochs),
                "data_fingerprint": store.manifest["fingerprint"],
                "training_site_digest": frame_digest(store.sites.iloc[store.ids("train")]),
                "feature_transform": "log dwell; per-position kmer standardization; central-minus-flank contrasts",
                "torch_version": str(torch.__version__)}, path)


def load_mil_bundle(path, device="cpu"):
    bundle = torch.load(path, map_location="cpu", weights_only=True)
    if bundle.get("format_version") != 1:
        raise ValueError("Unsupported advanced model format.")
    model = MILModel(ModelConfig(**bundle["model_config"]))
    model.load_state_dict(bundle["state_dict"])
    model.to(resolve_device(device)).eval()
    return model, KmerNormalizer.from_state_dict(bundle["normalizer"]), TrainConfig(**bundle["training_config"]), bundle


def mark_complete(run_dir, **details):
    required = ["model.pt", "evaluation.json", "test_predictions.parquet", "cv_complete.json"]
    if not all((run_dir / name).is_file() for name in required):
        raise ValueError("Cannot mark an incomplete model run as complete.")
    (run_dir / "complete.json").write_text(json.dumps({
        "completed_at": datetime.now(ZoneInfo("Asia/Singapore")).isoformat(), **details,
    }, indent=2))
