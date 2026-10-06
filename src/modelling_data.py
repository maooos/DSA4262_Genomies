"""Prepare pooled site observations with shared human-gene / construct groups.

Dataset names distinguish observations, but never distinguish the same human
gene across datasets. Synthetic mixture transcript IDs are not gene groups.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from .feature_engineering import FEATURE_COLUMNS, RAW_FEATURE_NAMES, SITE_KEY_COLUMNS
from .inference import features_from_reads


DATASETS = ("data0", "data1", "data2")
OBSERVATION_KEYS = ["dataset", *SITE_KEY_COLUMNS]
SPLIT_NAMES = ("train", "val", "test")
RAW_MEAN_COLUMNS = [f"raw_mean__{name}" for name in RAW_FEATURE_NAMES]
MODEL_FEATURE_COLUMNS = [*FEATURE_COLUMNS, *RAW_MEAN_COLUMNS]
READ_COLUMNS = [*SITE_KEY_COLUMNS, "sequence", "n_reads", "read_idx", "gene_id", "label", *RAW_FEATURE_NAMES]
SYNTHETIC_COLUMNS = ["geo_reference_id", "geo_reference_center_0based", "reference_mapping_status", "label_original"]


def read_columns(dataset):
    return READ_COLUMNS + (SYNTHETIC_COLUMNS if dataset == "data2" else [])


def load_site_metadata(paths):
    """Read first-read metadata and validate all inputs before feature extraction."""
    if set(paths) != set(DATASETS):
        raise ValueError("Pooled training requires data0, data1 and data2 paths.")
    parts = []
    for dataset in DATASETS:
        path = Path(paths[dataset])
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing annotated {dataset}: {path}. Run "
                "python -m scripts.prepare_datasets --download from the repository root."
            )
        audit_path = path.with_suffix(".audit.json")
        if audit_path.exists():
            audit = json.loads(audit_path.read_text())
            if audit.get("max_sites_requested") is not None:
                raise ValueError(f"{dataset}: use a complete dataset, not a max-sites sample.")
        parquet = pq.ParquetFile(path)
        missing = set(read_columns(dataset)) - set(parquet.schema_arrow.names)
        if missing:
            raise ValueError(f"{dataset}: missing annotated columns {sorted(missing)}")
        columns = [c for c in read_columns(dataset) if c not in ["read_idx", *RAW_FEATURE_NAMES]]
        sites = pd.read_parquet(path, columns=columns, filters=[("read_idx", "=", 0)])
        if sites.empty or sites.duplicated(SITE_KEY_COLUMNS).any():
            raise ValueError(f"{dataset}: empty metadata or duplicate site keys.")
        if sites.isna().any().any() or sites["gene_id"].str.strip().eq("").any():
            raise ValueError(f"{dataset}: missing site metadata / gene or construct annotation.")
        if not sites["label"].isin([0, 1]).all():
            raise ValueError(f"{dataset}: binary annotated labels are required.")
        if not sites["sequence"].str.fullmatch("[ACGT]{7}").all():
            raise ValueError(f"{dataset}: expected seven-base A/C/G/T sequences.")
        if not sites["n_reads"].gt(0).all() or sites["n_reads"].sum() != parquet.metadata.num_rows:
            raise ValueError(f"{dataset}: site coverage does not match the read count.")
        if dataset == "data2":
            proportions = sites["label_original"].astype(float)
            if not proportions.between(0, 1).all() or not sites["label"].eq(proportions.gt(0)).all():
                raise ValueError("data2: label must equal (label_original > 0).")
            if sites["geo_reference_id"].nunique() != 4:
                raise ValueError("data2: expected all four annotated synthetic constructs.")
            sites["group_id"] = "synthetic:" + sites["geo_reference_id"]
            sites["group_kind"] = "synthetic_construct"
        else:
            sites["group_id"] = "human:" + sites["gene_id"]
            sites["group_kind"] = "human_gene"
        sites["dataset"] = dataset
        sites["label"] = sites["label"].astype("int8")
        sites["central_5mer"] = sites["sequence"].str[1:6]
        parts.append(sites)
    result = pd.concat(parts, ignore_index=True).sort_values(OBSERVATION_KEYS).reset_index(drop=True)
    human = result.loc[result["group_kind"].eq("human_gene")]
    if not human.groupby("transcript_id")["group_id"].nunique().eq(1).all():
        raise ValueError("A human transcript maps to inconsistent genes across data0/data1.")
    if not human.groupby(SITE_KEY_COLUMNS)["sequence"].nunique().eq(1).all():
        raise ValueError("An overlapping human site has inconsistent sequence contexts.")
    return result


def assign_group_splits(sites, ratios=(0.70, 0.15, 0.15), seed=42):
    """Split the human gene union and synthetic constructs independently.

    With four constructs and default ratios this yields 2 / 1 / 1 constructs.
    All mixtures of a construct stay together. Labels never determine the split.
    """
    ratios = np.asarray(ratios, dtype=float)
    if len(ratios) != 3 or not np.isfinite(ratios).all() or (ratios <= 0).any() or not np.isclose(ratios.sum(), 1):
        raise ValueError("Provide three positive ratios that sum to 1.")
    if sites[["group_id", "group_kind"]].isna().any().any():
        raise ValueError("Every observation needs a group ID and group kind.")
    if not sites.groupby("group_id")["group_kind"].nunique().eq(1).all():
        raise ValueError("A group cannot belong to more than one group kind.")
    mapping = {}
    for kind, table in sites.groupby("group_kind", sort=True):
        groups = np.random.default_rng(seed).permutation(sorted(table["group_id"].unique()))
        cuts = np.floor(np.cumsum(ratios[:-1]) * len(groups)).astype(int)
        parts = np.split(groups, cuts)
        if any(len(part) == 0 for part in parts):
            raise ValueError(f"Too few {kind} groups for three nonempty partitions.")
        mapping.update({group: name for name, part in zip(SPLIT_NAMES, parts) for group in part})
    return sites["group_id"].map(mapping).rename("split")


def validate_splits(sites):
    """Reject leakage and partitions unable to support per-dataset evaluation."""
    if sites["split"].isna().any() or not sites["split"].isin(SPLIT_NAMES).all():
        raise ValueError("Unassigned / unknown partitions.")
    if sites.duplicated(OBSERVATION_KEYS).any():
        raise ValueError("Duplicate dataset/site observations.")
    if not sites.groupby("group_id")["split"].nunique().eq(1).all():
        raise ValueError("Group leakage across partitions.")
    human = sites.loc[sites["group_kind"].eq("human_gene")]
    if not human.groupby("transcript_id")["split"].nunique().eq(1).all():
        raise ValueError("Human transcript leakage across datasets / partitions.")
    synthetic = sites.loc[sites["group_kind"].eq("synthetic_construct")]
    if not synthetic.empty and not synthetic.groupby("geo_reference_id")["split"].nunique().eq(1).all():
        raise ValueError("Synthetic construct leakage across mixtures / partitions.")
    expected = pd.MultiIndex.from_product([DATASETS, SPLIT_NAMES], names=["dataset", "split"])
    summary = sites.groupby(["dataset", "split"]).agg(
        sites=("label", "size"), groups=("group_id", "nunique"),
        reads=("n_reads", "sum"), positives=("label", "sum"), positive_rate=("label", "mean"),
    ).reindex(expected)
    summary["negatives"] = summary["sites"] - summary["positives"]
    if summary.isna().any().any() or not summary[["positives", "negatives"]].gt(0).all().all():
        raise ValueError("Every dataset/split must contain both classes; review group composition without seed fishing.")
    return summary


def iter_complete_site_batches(path, dataset, batch_size=200_000):
    """Keep sites intact across read batches from the original-order Parquet."""
    if batch_size < 1:
        raise ValueError("batch_size must be positive.")
    carry = None
    for batch in pq.ParquetFile(path).iter_batches(batch_size=batch_size, columns=read_columns(dataset)):
        frame = batch.to_pandas()
        if carry is not None:
            frame = pd.concat([carry, frame], ignore_index=True)
        last_site = frame[SITE_KEY_COLUMNS].eq(frame[SITE_KEY_COLUMNS].iloc[-1]).all(axis=1)
        carry = frame.loc[last_site].copy()
        complete = frame.loc[~last_site].copy()
        if not complete.empty:
            yield complete
    if carry is not None and not carry.empty:
        yield carry


def build_pooled_features(paths, sites, batch_size=200_000):
    """Aggregate each dataset separately before attaching dataset-aware keys."""
    parts = []
    for dataset in DATASETS:
        processed = 0
        for number, reads in enumerate(iter_complete_site_batches(paths[dataset], dataset, batch_size), 1):
            if reads.isna().any().any() or not reads["label"].isin([0, 1]).all():
                raise ValueError(f"{dataset}: missing read values or nonbinary labels.")
            if reads.duplicated([*SITE_KEY_COLUMNS, "read_idx"]).any():
                raise ValueError(f"{dataset}: duplicate read indices within a site.")
            repeated = ["gene_id", "label"] + (SYNTHETIC_COLUMNS if dataset == "data2" else [])
            if not reads.groupby(SITE_KEY_COLUMNS)[repeated].nunique(dropna=False).eq(1).all().all():
                raise ValueError(f"{dataset}: inconsistent annotations within a site.")
            # Shared with inference; identifiers / labels never become predictors.
            features = features_from_reads(reads)
            features["dataset"] = dataset
            parts.append(features)
            processed += len(reads)
            if number % 10 == 0:
                print(f"{dataset}: processed {processed:,} reads", flush=True)
        if processed != pq.ParquetFile(paths[dataset]).metadata.num_rows:
            raise ValueError(f"{dataset}: incomplete read processing.")
        print(f"{dataset}: complete, {processed:,} reads", flush=True)
    features = pd.concat(parts, ignore_index=True).set_index(OBSERVATION_KEYS).sort_index()
    metadata = sites.set_index(OBSERVATION_KEYS).sort_index()
    if not features.index.is_unique or not features.index.equals(metadata.index):
        raise ValueError("Feature/site identities differ or a site appears in non-contiguous batches.")
    if not features["n_reads"].eq(metadata["n_reads"]).all():
        raise ValueError("Feature coverage differs from metadata.")
    if list(features.columns) != MODEL_FEATURE_COLUMNS or not np.isfinite(features.to_numpy()).all():
        raise ValueError("Unexpected feature schema or nonfinite values.")
    return features
