"""All-read disk cache, reproducible bags and fitting-partition-only normalization."""

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import shutil
from uuid import uuid4

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold, GroupShuffleSplit
import torch
from torch.utils.data import Dataset

from src.data_parser import file_sha256
from src.feature_engineering import RAW_FEATURE_NAMES, SITE_KEY_COLUMNS
from src.modelling_data import (
    DATASETS, OBSERVATION_KEYS, assign_group_splits, build_pooled_features,
    iter_complete_site_batches, load_site_metadata, validate_splits,
)

CACHE_VERSION = 1


def sequence_tokens(sequence):
    """Three ordered, overlapping 5-mers from the site's 7-mer; no fitted vocab."""
    bases = np.array(["ACGT".index(base) for base in sequence], dtype=np.int64)
    return np.array([bases[i:i + 5] @ (4 ** np.arange(4, -1, -1)) for i in range(3)])


def frame_digest(frame):
    return hashlib.sha256(frame.to_csv(index=False).encode()).hexdigest()


def assign_cv_folds(sites, n_splits=5, seed=42):
    """Same sorted-site GroupKFold policy as the existing pooled notebook."""
    membership = np.full(len(sites), -1, dtype=np.int16)
    train = np.flatnonzero(sites["split"].eq("train"))
    groups = sites.iloc[train]["group_id"]
    if groups.nunique() < n_splits:
        raise ValueError("Too few training groups for the requested CV folds.")
    splitter = GroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for fold, (fit, assess) in enumerate(splitter.split(train, groups=groups), 1):
        if set(groups.iloc[fit]) & set(groups.iloc[assess]):
            raise ValueError("CV group leakage.")
        for rows in (train[fit], train[assess]):
            if sites.iloc[rows]["label"].nunique() != 2:
                raise ValueError("A CV partition lacks one class; review groups without seed fishing.")
        membership[train[assess]] = fold
    return membership


@dataclass
class BagStore:
    cache_dir: Path
    sites: pd.DataFrame
    manifest: dict

    def __post_init__(self):
        self.reads = np.load(self.cache_dir / "reads.npy", mmap_mode="r")
        self.offsets = np.load(self.cache_dir / "offsets.npy", mmap_mode="r")
        self.sums = np.load(self.cache_dir / "sums.npy", mmap_mode="r")
        self.squares = np.load(self.cache_dir / "squares.npy", mmap_mode="r")
        self.tokens = np.stack(self.sites.sequence.map(sequence_tokens))

    def ids(self, split):
        return np.flatnonzero(self.sites["split"].eq(split))

    def raw_bag(self, row_id):
        return self.reads[self.offsets[row_id]:self.offsets[row_id + 1]]

    def indexed_metadata(self, ids):
        return self.sites.iloc[np.asarray(ids)].set_index(OBSERVATION_KEYS)


def prepare_store(processed_dir, *, seed=42, cv_folds=5, reference_run=None, batch_size=200_000):
    """Prepare once and reuse across all four notebooks; never fit normalization.

    Every source is fingerprinted. Partial builds never become reusable caches.
    Supplying reference_run checks exact split/fold membership and source hashes.
    """
    processed_dir = Path(processed_dir).resolve()
    paths = {name: processed_dir / "annotated" / f"{name}_reads.parquet" for name in DATASETS}
    sites = load_site_metadata(paths)
    sites["split"] = assign_group_splits(sites, seed=seed)
    validate_splits(sites)
    sites["assessment_fold"] = assign_cv_folds(sites, cv_folds, seed)
    sources = {name: {"path": str(path), "sha256": file_sha256(path)} for name, path in paths.items()}
    if reference_run is not None:
        reference_run = Path(reference_run)
        recorded = pd.read_parquet(reference_run / "split_manifest.parquet").set_index(OBSERVATION_KEYS).sort_index()
        current = sites.set_index(OBSERVATION_KEYS).sort_index()
        for column in ["split", "group_id", "label", "n_reads", "sequence"]:
            pd.testing.assert_series_equal(current[column], recorded[column], check_dtype=False)
        cv = pd.read_parquet(reference_run / "cv_fold_manifest.parquet").sort_index()
        pd.testing.assert_series_equal(current.loc[current.split.eq("train"), "assessment_fold"],
                                       cv["assessment_fold"], check_dtype=False)
        lineage = json.loads((reference_run / "cv_started.json").read_text())
        if lineage["source_parquet_sha256"] != {k: v["sha256"] for k, v in sources.items()}:
            raise ValueError("Reference run used different annotated inputs.")
    identity = {"cache_version": CACHE_VERSION, "sources": sources, "seed": seed,
                "cv_folds": cv_folds, "sites_digest": frame_digest(sites),
                "features": RAW_FEATURE_NAMES, "transform": "log dwell at 0/3/6; raw sd and mean"}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    cache_dir = processed_dir / "advanced" / "cache" / fingerprint[:20]
    if (cache_dir / "complete.json").is_file():
        manifest = json.loads((cache_dir / "complete.json").read_text())
        if manifest["fingerprint"] != fingerprint:
            raise ValueError("Cache identity mismatch.")
        for name, size in manifest["file_sizes"].items():
            if not (cache_dir / name).is_file() or (cache_dir / name).stat().st_size != size:
                raise ValueError(f"Incomplete cache file: {cache_dir / name}")
        print("Reusing all-read cache:", cache_dir)
        return BagStore(cache_dir, sites, manifest)
    if cache_dir.exists():
        raise ValueError(f"Incomplete cache directory: {cache_dir}; inspect/remove it before retrying.")
    cache_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = cache_dir.with_name(f".{cache_dir.name}.{uuid4().hex}.partial")
    staging.mkdir()
    try:
        offsets = np.r_[0, sites.n_reads.to_numpy(dtype=np.int64).cumsum()]
        np.save(staging / "offsets.npy", offsets)
        reads_map = np.lib.format.open_memmap(staging / "reads.npy", mode="w+", dtype="float32", shape=(int(offsets[-1]), 9))
        sums = np.zeros((len(sites), 9), dtype=np.float64)
        squares = np.zeros_like(sums)
        lookup = {tuple(key): i for i, key in enumerate(sites[OBSERVATION_KEYS].itertuples(index=False, name=None))}
        seen = np.zeros(len(sites), dtype=bool)
        for name, path in paths.items():
            processed = 0
            for batch in iter_complete_site_batches(path, name, batch_size):
                if batch.duplicated([*SITE_KEY_COLUMNS, "read_idx"]).any():
                    raise ValueError("Duplicate read indices within a site.")
                for key, frame in batch.groupby(SITE_KEY_COLUMNS, sort=False):
                    row_id = lookup[(name, *key)]
                    expected = sites.iloc[row_id]
                    if seen[row_id] or len(frame) != expected.n_reads:
                        raise ValueError("Duplicate or incomplete site in read cache.")
                    for column in ["sequence", "gene_id", "label", "n_reads"]:
                        if not frame[column].eq(expected[column]).all():
                            raise ValueError(f"Read/site metadata mismatch: {column}")
                    if name == "data2":
                        for column in ["geo_reference_id", "label_original", "reference_mapping_status"]:
                            if not frame[column].eq(expected[column]).all():
                                raise ValueError(f"Inconsistent synthetic annotations: {column}")
                    values = frame[RAW_FEATURE_NAMES].to_numpy(dtype=np.float64, copy=True)
                    if not np.isfinite(values).all() or (values[:, [0, 3, 6]] <= 0).any():
                        raise ValueError("Read signals must be finite and dwell times positive.")
                    values[:, [0, 3, 6]] = np.log(values[:, [0, 3, 6]])
                    reads_map[offsets[row_id]:offsets[row_id + 1]] = values
                    sums[row_id], squares[row_id] = values.sum(0), np.square(values).sum(0)
                    seen[row_id] = True
                processed += len(batch)
                print(f"Cache {name}: {processed:,} reads", flush=True)
        if not seen.all():
            raise ValueError("Some labelled sites have no cached reads.")
        reads_map.flush()
        del reads_map
        np.save(staging / "sums.npy", sums)
        np.save(staging / "squares.npy", squares)
        sites.to_parquet(staging / "sites.parquet", index=False)
        manifest = {**identity, "fingerprint": fingerprint,
                    "file_sizes": {p.name: p.stat().st_size for p in staging.iterdir()}}
        (staging / "complete.json").write_text(json.dumps(manifest, indent=2))
        staging.rename(cache_dir)
    except BaseException:
        shutil.rmtree(staging)
        raise
    return BagStore(cache_dir, sites, manifest)


@dataclass
class KmerNormalizer:
    mean: np.ndarray
    scale: np.ndarray

    @classmethod
    def fit(cls, store, ids):
        """Moment estimates from every read in ONLY these fitting sites."""
        ids = np.asarray(ids, dtype=int)
        if not len(ids) or not store.sites.iloc[ids].split.eq("train").all():
            raise ValueError("Fit normalization on nonempty outer-training site IDs only.")
        count = store.sites.iloc[ids].n_reads.to_numpy(dtype=float)
        means, scales = [], []
        for position in range(3):
            cols = slice(3 * position, 3 * position + 3)
            n = np.bincount(store.tokens[ids, position], weights=count, minlength=1024)
            total, squares = np.zeros((1024, 3)), np.zeros((1024, 3))
            np.add.at(total, store.tokens[ids, position], store.sums[ids, cols])
            np.add.at(squares, store.tokens[ids, position], store.squares[ids, cols])
            fallback_mean = total.sum(0) / count.sum()
            fallback_var = np.maximum(squares.sum(0) / count.sum() - fallback_mean ** 2, 1e-6)
            mean = np.divide(total, n[:, None], out=np.tile(fallback_mean, (1024, 1)), where=n[:, None] > 0)
            second = np.divide(squares, n[:, None], out=np.tile(fallback_var + fallback_mean ** 2, (1024, 1)), where=n[:, None] > 0)
            means.append(mean)
            scales.append(np.sqrt(np.maximum(second - mean ** 2, 1e-6)))
        return cls(np.array(means, dtype=np.float32), np.array(scales, dtype=np.float32))

    def transform(self, reads, tokens):
        # Retain flank order within a read; concatenate central-minus-flank contrasts.
        mean = self.mean[np.arange(3), tokens].reshape(9)
        scale = self.scale[np.arange(3), tokens].reshape(9)
        normalized = (np.asarray(reads, dtype=np.float32) - mean) / scale
        return np.concatenate([normalized, normalized[:, 3:6] - normalized[:, :3],
                               normalized[:, 3:6] - normalized[:, 6:9]], axis=1).astype(np.float32)

    def state_dict(self):
        return {"mean": torch.from_numpy(self.mean.copy()), "scale": torch.from_numpy(self.scale.copy())}

    @classmethod
    def from_state_dict(cls, state):
        return cls(state["mean"].cpu().numpy(), state["scale"].cpu().numpy())


class ReadBags(Dataset):
    """One label per site; padded short bags never duplicate real reads."""

    def __init__(self, store, ids, normalizer, max_reads=32, seed=42, draw=0):
        if max_reads < 1:
            raise ValueError("max_reads must be positive.")
        self.store, self.ids, self.normalizer = store, np.asarray(ids), normalizer
        self.max_reads, self.seed, self.draw = max_reads, seed, draw

    def __len__(self):
        return len(self.ids)

    def __getitem__(self, index):
        row_id = int(self.ids[index])
        reads = self.store.raw_bag(row_id)
        rng = np.random.default_rng(np.random.SeedSequence([self.seed, self.draw, row_id]))
        if len(reads) > self.max_reads:
            reads = reads[rng.choice(len(reads), self.max_reads, replace=False)]
        count = len(reads)
        features = np.zeros((self.max_reads, 15), dtype=np.float32)
        features[:count] = self.normalizer.transform(reads, self.store.tokens[row_id])
        site = self.store.sites.iloc[row_id]
        return {"features": torch.from_numpy(features),
                "tokens": torch.tensor(self.store.tokens[row_id], dtype=torch.long),
                "mask": torch.arange(self.max_reads) < count,
                "label": torch.tensor(float(site.label)),
                "dataset": torch.tensor(DATASETS.index(site.dataset)),
                "fraction": torch.tensor(float(site.label_original) if site.dataset == "data2" else 0.),
                "row_id": torch.tensor(row_id)}


def inner_group_split(store, fit_ids, seed=42, stop_fraction=0.15):
    fit_ids = np.asarray(fit_ids)
    if not store.sites.iloc[fit_ids].split.eq("train").all():
        raise ValueError("Early stopping partitions must come from outer train.")
    groups = store.sites.iloc[fit_ids].group_id
    a, b = next(GroupShuffleSplit(n_splits=1, test_size=stop_fraction, random_state=seed).split(fit_ids, groups=groups))
    if set(groups.iloc[a]) & set(groups.iloc[b]):
        raise ValueError("Early-stopping group leakage.")
    for indices in (a, b):
        if store.sites.iloc[fit_ids[indices]].label.nunique() != 2:
            raise ValueError("Inner fitting/stopping partition lacks a class; review group composition.")
    return fit_ids[a], fit_ids[b]


def load_tabular_features(store):
    """Site-local H04 features, cached once; no population statistics are fitted."""
    path = store.cache_dir / "h04_features.parquet"
    if path.exists():
        features = pd.read_parquet(path)
    else:
        paths = {name: Path(entry["path"]) for name, entry in store.manifest["sources"].items()}
        features = build_pooled_features(paths, store.sites)
        temporary = path.with_name(f".{uuid4().hex}.parquet")
        features.to_parquet(temporary)
        temporary.rename(path)
    expected = store.sites.set_index(OBSERVATION_KEYS).index
    if not features.index.equals(expected):
        raise ValueError("H04 features do not align with read bags.")
    return features
